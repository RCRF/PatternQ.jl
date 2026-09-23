# HTTP transport to the Pattern Data Commons query service.
#
# Queries are Dicts in the JSON form parsed by the service's datalog-json-parser:
# variables are "?x" strings, attributes/keywords ":ns/name" strings, clauses
# vectors, pull pattern maps Dicts:
#
#     Dict(":find" => ["?id"], ":where" => [["?s", ":subject/id", "?id"]])
#
# The implicit database "\$" is prepended to ":in" if absent.
#
# Endpoints (bearer API token):
#   POST /query/<db>   Accept text/plain -> presigned URL of the gzipped, S3-cached
#                      result; Accept application/json -> inline JSON, cache skipped
#   POST /datoms/<db>  -> {"datoms_chunk": [...], "basis_t": ...}
#   POST /matrix/<db>/<key> -> presigned URL of a gzipped TSV matrix
#   GET  /api-v1/list/datasets

const Query = AbstractDict

headers(accept="text/plain") = ["Authorization" => "Bearer $(api_token())", "Accept" => accept,
                                "User-Agent" => "patternq-julia", "Content-Type" => "application/json"]

function raise_for(resp, what)
    resp.status == 200 && return
    body = String(copy(resp.body))
    msg = body
    try
        parsed = JSON3.read(body)
        parsed isa JSON3.Object && haskey(parsed, :error) && (msg = string(parsed[:error]))
    catch
    end
    resp.status == 401 && (msg = "not authorized; check PATTERNQ_API_KEY / set_token(). " * msg)
    resp.status == 403 && (msg = "forbidden; your API key may not have access to this dataset. " * msg)
    error("$what failed (HTTP $(resp.status)): $msg")
end

function fetch_presigned(url)
    resp = HTTP.get(strip(url); status_exception=false, retry=true, retries=4, request_timeout=300)
    raise_for(resp, "Result download")
    b = resp.body
    length(b) >= 2 && b[1] == 0x1f && b[2] == 0x8b ? transcode(GzipDecompressor, b) : b
end

function normalize_query(q::Query)
    out = Dict{String,Any}()
    for (k, v) in q
        ks = string(k)
        out[startswith(ks, ":") ? ks : ":" * ks] = v
    end
    ins = get(out, ":in", Any[])
    ins = ins === nothing ? Any[] : collect(Any, ins)
    (isempty(ins) || ins[1] != "\$") && pushfirst!(ins, "\$")
    out[":in"] = ins
    out
end

"""
    query_body(q; args=[], timeout=30, refresh_cache=false)

The JSON request body sent to the query service.
"""
function query_body(q::Query; args=Any[], timeout=30, refresh_cache=false)
    body = Dict{String,Any}("query" => normalize_query(q), "timeout" => round(Int, timeout * 1000))
    isempty(args) || (body["args"] = collect(Any, args))
    refresh_cache && (body["refresh-cache"] = true)
    body
end

"""
    query(q; args=[], db=nothing, timeout=30, cache=true, refresh_cache=false, print_json=false)

Run a query and return the parsed response: a Dict with "query_result",
"basis_t" and "db_name". Most users want `do_query`, which returns a DataFrame.

`cache=true` uses the service's S3 result cache (presigned URL to a gzipped
cached result, computed and cached on a miss); `cache=false` returns the result
inline and skips the cache; `refresh_cache=true` recomputes and re-caches.
"""
function query(q::Query; args=Any[], db=nothing, timeout=30, cache::Bool=true,
               refresh_cache::Bool=false, print_json::Bool=false)
    db = ensure_db(db)
    body = JSON3.write(query_body(q; args=args, timeout=timeout, refresh_cache=refresh_cache))
    print_json && println(body)
    accept = cache ? "text/plain" : "application/json"
    resp = HTTP.post("$(query_server())/query/$(db)", headers(accept), body;
                     status_exception=false, request_timeout=timeout + 30)
    raise_for(resp, "Query")
    payload = strip(String(copy(resp.body)))
    raw = startswith(payload, "{") ? payload : String(fetch_presigned(payload))
    res = tonative(JSON3.read(raw))
    get(res, "error", nothing) === nothing || error("Query error: $(res["error"])")
    res["db_name"] = db
    res
end

# Pulled refs without a nested pattern come back as {":db/id": n}. Resolve those
# that are enums (have a :db/ident) with one extra query, so enum values read
# as names whichever pull pattern was used.
function resolve_enum_refs(rows, db)
    eids = Set{Any}()
    collect_eids(x) = if is_eid_map(x)
        push!(eids, x[":db/id"])
    elseif x isa AbstractDict
        foreach(collect_eids, values(x))
    elseif x isa AbstractVector
        foreach(collect_eids, x)
    end
    collect_eids(rows)
    isempty(eids) && return rows
    res = query(Dict(":find" => ["?e", "?ident"], ":in" => [["?e", "..."]],
                     ":where" => [["?e", ":db/ident", "?ident"]]);
                args=[collect(eids)], db=db)
    idents = Dict(r[1] => r[2] for r in res["query_result"])
    isempty(idents) && return rows
    replace_refs(x) = if is_eid_map(x) && haskey(idents, x[":db/id"])
        OrderedDict{String,Any}(":db/ident" => idents[x[":db/id"]])
    elseif x isa AbstractDict
        OrderedDict{String,Any}(k => replace_refs(v) for (k, v) in x)
    elseif x isa AbstractVector
        Any[replace_refs(v) for v in x]
    else
        x
    end
    replace_refs(rows)
end

"""
    do_query(q; args=[], db=nothing, exclude_ids=true, kwargs...)

Run a query and return a DataFrame. Column names come from the :find variables
("?sample-id" -> "sample_id"); pull expressions are flattened (see
`flatten_pull`). The result carries provenance (see `provenance`).
"""
function do_query(q::Query; args=Any[], db=nothing, exclude_ids::Bool=true, kwargs...)
    res = query(q; args=args, db=db, kwargs...)
    find = normalize_query(q)[":find"]
    rows = res["query_result"]
    any(is_pull, find) && (rows = resolve_enum_refs(rows, res["db_name"]))
    df = result_to_df(rows, find; exclude_ids=exclude_ids)
    with_provenance(df, res["db_name"], get(res, "basis_t", nothing))
end

"""
    across_dbs(q, dbs; args=[], kwargs...)

Run the same query against several databases and concatenate the results,
adding a `db` column. Each dataset is its own database, so cohort comparisons
run per database.
"""
function across_dbs(q::Query, dbs; args=Any[], kwargs...)
    parts = DataFrame[]
    for db in dbs
        df = do_query(q; args=args, db=db, kwargs...)
        if nrow(df) > 0
            df[!, "db"] .= db
            push!(parts, df)
        end
    end
    isempty(parts) ? DataFrame() : vcat(parts...; cols=:union)
end

"""
    datoms(index, components=[]; db=nothing, offset=0, limit=1000)

Raw datoms from an index ("eavt", "aevt", "avet", "vaet"): columns e, a, v, tx.
"""
function datoms(index::AbstractString, components=Any[]; db=nothing, offset=0, limit=1000, timeout=30)
    db = ensure_db(db)
    body = JSON3.write(Dict("index" => ":" * lstrip(index, ':'), "components" => components,
                            "offset" => offset, "limit" => limit))
    resp = HTTP.post("$(query_server())/datoms/$(db)", headers("application/json"), body;
                     status_exception=false, request_timeout=timeout + 2)
    raise_for(resp, "Datoms request")
    res = tonative(JSON3.read(resp.body))
    chunk = get(res, "datoms_chunk", Any[])
    df = DataFrame(e=Any[d[":e"] for d in chunk], a=Any[d[":a"] for d in chunk],
                   v=Any[d[":v"] for d in chunk], tx=Any[d[":tx"] for d in chunk])
    for c in names(df)
        df[!, c] = tidy_column(df[!, c])
    end
    with_provenance(df, db, get(res, "basis_t", nothing))
end

"""
    list_datasets()

Datasets available to your API key: dataset, db (current database name),
patient_count, sample_count, assays, tags.
"""
function list_datasets()
    resp = HTTP.get("$(query_server())/api-v1/list/datasets", headers("application/json");
                    status_exception=false, request_timeout=60)
    raise_for(resp, "Listing datasets")
    ds = get(tonative(JSON3.read(resp.body)), "datasets", Any[])
    getdb(d) = (x = get(d, "dataset/database", nothing); x isa AbstractDict ? get(x, "database/name", missing) : missing)
    DataFrame(dataset=[get(d, "dataset/name", missing) for d in ds],
              db=[getdb(d) for d in ds],
              patient_count=[something(get(d, "dataset/patient-count", nothing), missing) for d in ds],
              sample_count=[something(get(d, "dataset/sample-count", nothing), missing) for d in ds],
              assays=[collect(String, something(get(d, "dataset/assays", nothing), String[])) for d in ds],
              tags=[collect(String, something(get(d, "dataset/tags", nothing), String[])) for d in ds])
end

"""
    resolve_db(dataset)

Resolve a dataset name (stable, e.g. "tcga-uvm") to its current database name
(changes on re-import). A database name is passed through.
"""
function resolve_db(dataset::AbstractString)
    ds = list_datasets()
    i = findfirst(==(dataset), coalesce.(ds.dataset, ""))
    i !== nothing && !ismissing(ds.db[i]) && return ds.db[i]
    dataset in skipmissing(ds.db) && return String(dataset)
    error("Unknown dataset '$dataset'")
end

"""
    measurement_matrix(matrix_key; db=nothing)

Download a measurement matrix (e.g. single-cell counts) by its backing-file key
(see `measurement_matrices`), as a DataFrame read from the TSV.
"""
function measurement_matrix(matrix_key::AbstractString; db=nothing)
    db = ensure_db(db)
    resp = HTTP.post("$(query_server())/matrix/$(db)/$(matrix_key)", headers("text/plain"), "{}";
                     status_exception=false, request_timeout=120)
    raise_for(resp, "Matrix request")
    content = String(fetch_presigned(String(copy(resp.body))))
    lines = split(chomp(content), '\n')
    header = split(lines[1], '\t')
    rows = [split(l, '\t') for l in lines[2:end]]
    df = DataFrame()
    for (j, h) in enumerate(header)
        col = [r[j] for r in rows]
        nums = tryparse.(Float64, col)
        df[!, String(h)] = all(!isnothing, nums) ? Float64.(nums) : String.(col)
    end
    df
end
