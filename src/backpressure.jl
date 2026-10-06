# Client side of the commons API back-pressure contract (unify-central
# dev-docs/API_BACKPRESSURE.md), the same in the R, Python, Clojure and Julia
# libraries:
#
# - 429 and 503 are retried, up to max_retries times, when the body is a
#   retryable application/problem+json, or is not problem+json at all (e.g.
#   from a load balancer). Every other status, 502 and 504 included, is
#   returned as is. HTTP.jl's own retries are off for API calls (retry=false).
# - The wait is Retry-After when the server sends it, otherwise exponential
#   backoff with full jitter: rand() * min(max_backoff, 2^attempt) seconds.
# - Self-throttling: a RateLimit header with r=0 makes the next call wait t
#   seconds, and at most max_concurrency /query and /datoms calls run at once
#   in this process (the server's per-key cap is 4).
# - Retries are logged with @info, with the problem type and the wait.
#
# Presigned S3 downloads are not API calls and don't go through here.

const PROBLEM_PREFIX = "urn:pattern-data-commons:problem:"
const DEFAULT_MAX_TIMEOUT_MS = 120000

const _bp_lock = ReentrantLock()
const _retry_policy = Dict{Symbol,Any}(:max_retries => 5, :max_backoff => 60.0, :max_concurrency => 4)
const _slots = Ref(Base.Semaphore(4))
const _next_allowed = Ref(0.0)  # time() before which no call is sent
const _max_timeout_ms = Ref(DEFAULT_MAX_TIMEOUT_MS)

# indirection so tests don't really sleep
const _sleep = Ref{Function}(sleep)

"""
    ThrottledError

The commons API kept throttling a call past the retry limit. Fields:
`problem_type` (or `nothing`), `status`, `msg`.
"""
struct ThrottledError <: Exception
    msg::String
    problem_type::Union{Nothing,String}
    status::Int
end
Base.showerror(io::IO, e::ThrottledError) = print(io, e.msg)

"""
    set_retry_policy!(; max_retries, max_backoff, max_concurrency)

Set how this session handles API back pressure. Keywords left out keep their
current value. Defaults: `max_retries=5` retries per call, `max_backoff=60`
seconds (cap on the jittered backoff when the server sends no Retry-After),
`max_concurrency=4` concurrent /query and /datoms calls (the server's per-key
cap). Returns the previous policy.
"""
function set_retry_policy!(; max_retries=nothing, max_backoff=nothing, max_concurrency=nothing)
    lock(_bp_lock) do
        old = copy(_retry_policy)
        max_retries === nothing || (_retry_policy[:max_retries] = Int(max_retries))
        max_backoff === nothing || (_retry_policy[:max_backoff] = Float64(max_backoff))
        if max_concurrency !== nothing && Int(max_concurrency) != _retry_policy[:max_concurrency]
            max_concurrency >= 1 || throw(ArgumentError("max_concurrency must be at least 1"))
            _retry_policy[:max_concurrency] = Int(max_concurrency)
            _slots[] = Base.Semaphore(Int(max_concurrency))
        end
        old
    end
end

"""
    retry_policy()

The current API back-pressure policy (see `set_retry_policy!`).
"""
retry_policy() = copy(_retry_policy)

# The server's cap on requested query timeouts, from the last
# PDC-Query-Max-Timeout-Ms header seen (120000 until one is).
max_timeout_ms() = _max_timeout_ms[]

function _problem(resp)
    startswith(HTTP.header(resp, "Content-Type", ""), "application/problem+json") || return nothing
    try
        p = JSON3.read(String(copy(resp.body)))
        p isa JSON3.Object ? p : nothing
    catch
        nothing
    end
end

function _retry_after(resp, problem)
    h = tryparse(Float64, strip(HTTP.header(resp, "Retry-After", "")))
    h === nothing || return max(0.0, h)
    problem !== nothing && get(problem, :retry_after, nothing) isa Number && return max(0.0, Float64(problem[:retry_after]))
    nothing
end

function _note_headers(resp)
    rl = HTTP.header(resp, "RateLimit", "")
    r = match(r"\br=(\d+)", rl)
    t = match(r"\bt=(\d+)", rl)
    if r !== nothing && t !== nothing && parse(Int, r[1]) == 0
        lock(_bp_lock) do
            _next_allowed[] = max(_next_allowed[], time() + parse(Int, t[1]))
        end
    end
    mt = tryparse(Int, strip(HTTP.header(resp, "PDC-Query-Max-Timeout-Ms", "")))
    mt === nothing || (_max_timeout_ms[] = mt)
end

function _wait_for_rate_limit()
    wait = _next_allowed[] - time()
    if wait > 0
        @info "patternq: rate limit reached (RateLimit r=0); waiting $(round(wait; digits=1)) s before the next call"
        _sleep[](wait)
    end
end

# Send an API request (do_request performs it) under the back-pressure
# contract. limited: the call counts toward the per-key concurrency cap
# (/query, /datoms). Returns the final response, which may still be an error
# for raise_for; throws ThrottledError when throttling outlasts max_retries.
function send_api(do_request, what::AbstractString; limited::Bool=false)
    attempt = 0
    while true
        _wait_for_rate_limit()
        resp = limited ? Base.acquire(do_request, _slots[]) : do_request()
        _note_headers(resp)
        resp.status in (429, 503) || return resp
        problem = _problem(resp)
        problem !== nothing && get(problem, :retryable, false) !== true && return resp
        ptype = problem === nothing ? nothing : (t = get(problem, :type, nothing); t === nothing ? nothing : String(t))
        label = something(ptype, "HTTP $(resp.status)")
        max_retries = _retry_policy[:max_retries]
        if attempt >= max_retries
            detail = problem !== nothing && haskey(problem, :detail) ? String(problem[:detail]) :
                     first(strip(String(copy(resp.body))), 500)
            throw(ThrottledError("$what throttled by the commons API ($label) and still throttled after " *
                                 "$attempt retries: $detail", ptype, resp.status))
        end
        attempt += 1
        wait = _retry_after(resp, problem)
        wait === nothing && (wait = rand() * min(_retry_policy[:max_backoff], 2.0^attempt))
        @info "patternq: $what throttled ($label, HTTP $(resp.status)); retry $attempt of $max_retries in $(round(wait; digits=1)) s"
        _sleep[](wait)
    end
end
