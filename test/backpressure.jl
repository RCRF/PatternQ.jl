# API back-pressure contract (unify-central dev-docs/API_BACKPRESSURE.md),
# against a local stub server.
const HTTP = PatternQ.HTTP
const JSON3 = PatternQ.JSON3

const PROBLEM = PatternQ.PROBLEM_PREFIX

ok_resp(headers=["RateLimit" => "\"api\";r=99;t=0"]) =
    (200, ["Content-Type" => "application/json"; headers], Dict("query_result" => [[1]], "basis_t" => 7))

function throttled(; kind="rate-limited", status=429, retry_after=nothing, retryable=true)
    headers = Pair{String,String}["Content-Type" => "application/problem+json"]
    retry_after === nothing || push!(headers, "Retry-After" => string(retry_after))
    (status, headers, Dict("type" => PROBLEM * kind, "title" => "throttled", "status" => status,
                           "detail" => "$kind detail", "retryable" => retryable))
end

# Serves responses in order (the last one repeats); records paths and the peak
# number of requests in flight.
mutable struct Stub
    server::Any
    paths::Vector{String}
    inflight::Int
    peak::Int
end

function with_stub(f, responses; delay=0.0)
    queue = collect(Any, responses)
    lk = ReentrantLock()
    stub = Stub(nothing, String[], 0, 0)
    stub.server = HTTP.serve!("127.0.0.1", 0) do req
        resp = lock(lk) do
            push!(stub.paths, req.target)
            stub.inflight += 1
            stub.peak = max(stub.peak, stub.inflight)
            length(queue) > 1 ? popfirst!(queue) : queue[1]
        end
        delay > 0 && sleep(delay)
        status, headers, body = resp
        lock(() -> stub.inflight -= 1, lk)
        HTTP.Response(status, headers, body isa AbstractString ? body : JSON3.write(body))
    end
    sleeps = Float64[]
    old_policy = retry_policy()
    old_server = set_query_server("http://127.0.0.1:$(HTTP.port(stub.server))")
    set_token("test-token")
    PatternQ._sleep[] = s -> push!(sleeps, s)
    try
        f(stub, sleeps)
    finally
        close(stub.server)
        PatternQ._sleep[] = sleep
        set_query_server(old_server)
        set_token(nothing)
        set_retry_policy!(; old_policy...)
        PatternQ._next_allowed[] = 0.0
        PatternQ._max_timeout_ms[] = PatternQ.DEFAULT_MAX_TIMEOUT_MS
    end
end

const Q = Dict(":find" => ["?x"], ":where" => [["?x", ":db/ident"]])
run_query(; kw...) = PatternQ.query(Q; db="db1", cache=false, kw...)

@testset "API back pressure" begin
    @testset "429 and 503 are retried after Retry-After" begin
        with_stub([throttled(retry_after=3), throttled(kind="overloaded", status=503, retry_after=1), ok_resp()]) do s, sleeps
            res = @test_logs (:info, r"rate-limited") (:info, r"overloaded") run_query()
            @test res["basis_t"] == 7
            @test sleeps == [3.0, 1.0]
            @test length(s.paths) == 3
        end
    end
    @testset "jittered backoff without Retry-After" begin
        plain = (429, ["Content-Type" => "text/plain"], "slow down")
        with_stub([plain, plain, plain, ok_resp()]) do s, sleeps
            set_retry_policy!(max_backoff=3)
            @test_logs min_level = Base.CoreLogging.Warn run_query()
            @test length(sleeps) == 3
            @test all(0 .<= sleeps .<= min.(3, 2.0 .^ (1:3)))
        end
    end
    @testset "gives up after max_retries naming the problem" begin
        with_stub([throttled(retry_after=1)]) do s, sleeps
            set_retry_policy!(max_retries=2)
            e = try
                Base.CoreLogging.with_logger(run_query, Base.CoreLogging.NullLogger())
                nothing
            catch err
                err
            end
            @test e isa ThrottledError
            @test e.problem_type == PROBLEM * "rate-limited"
            @test occursin("rate-limited detail", e.msg) && occursin("2 retries", e.msg)
            @test length(s.paths) == 3
        end
    end
    @testset "non-retryable problems, 502, 504 and 400 are not retried" begin
        for r in [throttled(retryable=false), (502, Pair{String,String}[], "bad gateway"),
                  (504, Pair{String,String}[], "timeout"),
                  (400, ["Content-Type" => "application/json"], Dict("error" => "bad query"))]
            with_stub([r, ok_resp()]) do s, sleeps
                @test_throws r"failed \(HTTP" run_query()
                @test length(s.paths) == 1
                @test isempty(sleeps)
            end
        end
    end
    @testset "a query timeout says to narrow or page the query" begin
        with_stub([(400, ["Content-Type" => "application/json"],
                    Dict("error" => "Query canceled: timeout elapsed", "timeout" => true,
                         "type" => PROBLEM * "query-timeout"))]) do s, sleeps
            @test_throws r"Narrow the query or page it" run_query()
        end
    end
    @testset "a query-too-broad refusal is surfaced and not retried" begin
        with_stub([(400, ["Content-Type" => "application/json"],
                    Dict("error" => "Query too broad: clause [?e ?a ?v] binds neither attribute nor entity; bind the attribute", "type" => PROBLEM * "query-too-broad",
                         "reason" => "full-scan", "clause" => "[?e ?a ?v]")), ok_resp()]) do s, sleeps
            @test_throws r"bind the attribute" run_query()
            @test length(s.paths) == 1 && isempty(sleeps)
        end
    end
    @testset "RateLimit r=0 makes the next call wait" begin
        with_stub([ok_resp(["RateLimit" => "\"api\";r=0;t=1"]), ok_resp()]) do s, sleeps
            run_query()
            @test isempty(sleeps)
            @test_logs (:info, r"rate limit reached") run_query()
            @test length(sleeps) == 1 && 0 < sleeps[1] <= 1
        end
    end
    @testset "listing datasets and datoms go through back pressure" begin
        with_stub([throttled(retry_after=1), (200, ["Content-Type" => "application/json"], Dict("datasets" => [])),
                   throttled(retry_after=1),
                   (200, ["Content-Type" => "application/json"], Dict("datoms_chunk" => [], "basis_t" => 1))]) do s, sleeps
            @test_logs min_level = Base.CoreLogging.Warn begin
                @test nrow(list_datasets()) == 0
                @test nrow(datoms("eavt"; db="db1")) == 0
            end
            @test sleeps == [1.0, 1.0]
        end
    end
    @testset "concurrency is capped per process" begin
        with_stub([ok_resp()]; delay=0.15) do s, sleeps
            set_retry_policy!(max_concurrency=2)
            tasks = [Threads.@spawn run_query() for _ in 1:6]
            foreach(fetch, tasks)
            @test length(s.paths) == 6
            @test s.peak == 2
        end
    end
    @testset "timeouts above the advertised cap warn" begin
        with_stub([ok_resp(["PDC-Query-Max-Timeout-Ms" => "60000"])]) do s, sleeps
            @test_logs (:warn, r"cap of 120 s") run_query(timeout=200)
            @test PatternQ.max_timeout_ms() == 60000
            @test_logs (:warn, r"cap of 60 s") run_query(timeout=90)
        end
    end
    @testset "policy defaults" begin
        @test retry_policy() == Dict(:max_retries => 5, :max_backoff => 60.0, :max_concurrency => 4)
    end
end
