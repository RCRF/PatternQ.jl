# Session configuration: query service endpoint, API token, default database.
# Resolution order: values set in-session, then the PATTERNQ_ENDPOINT /
# PATTERNQ_API_KEY environment variables, then (endpoint only) the dev server.

const DEFAULT_ENDPOINT = "https://data-commons.rcrf-dev.org"

const _config = Dict{Symbol,Union{Nothing,String}}(:endpoint => nothing, :token => nothing, :db => nothing)

"""
    query_server()

The query service endpoint in use.
"""
function query_server()
    ep = _config[:endpoint]
    ep !== nothing && return ep
    env = get(ENV, "PATTERNQ_ENDPOINT", "")
    startswith(env, "http") ? rstrip(env, '/') : DEFAULT_ENDPOINT
end

"""
    set_query_server(url)

Set the query service endpoint for this session; returns the previous value.
"""
function set_query_server(url::Union{Nothing,AbstractString})
    old = _config[:endpoint]
    _config[:endpoint] = url === nothing ? nothing : String(rstrip(url, '/'))
    old
end

"""
    set_token(token)

Set the API token for this session, overriding `PATTERNQ_API_KEY`.
"""
function set_token(token::Union{Nothing,AbstractString})
    old = _config[:token]
    _config[:token] = token === nothing ? nothing : String(token)
    old
end

function api_token()
    t = _config[:token]
    t = (t === nothing || isempty(t)) ? get(ENV, "PATTERNQ_API_KEY", "") : t
    isempty(t) && error("No API token: set PATTERNQ_API_KEY in the environment or call set_token()")
    t
end

"""
    set_db(db)

Set the default database for this session. Every dataset is its own database,
so the database name selects the dataset; `resolve_db` maps a dataset name
(e.g. "tcga-uvm") to its current database name.
"""
function set_db(db::Union{Nothing,AbstractString})
    old = _config[:db]
    _config[:db] = db === nothing ? nothing : String(db)
    old
end

"""
    current_db()

The default database for this session, or `nothing`.
"""
current_db() = _config[:db]

function ensure_db(db)
    db = db === nothing ? _config[:db] : db
    db === nothing && error("No database given: pass db=... or call set_db()")
    String(db)
end
