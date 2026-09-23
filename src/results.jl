# Converting query results to DataFrames, and result provenance.

# JSON3 values -> plain Julia (Dict{String,Any}, Vector{Any}, nothing for null)
tonative(x::JSON3.Object) = OrderedDict{String,Any}(String(k) => tonative(v) for (k, v) in pairs(x))
tonative(x::JSON3.Array) = Any[tonative(v) for v in x]
tonative(x) = x

"""
    clean_name(x)

Datalog variable or keyword to a column name: `"?sample-id"` -> `"sample_id"`,
`":subject/id"` -> `"subject_id"`.
"""
clean_name(x) = replace(replace(string(x), r"^[?:]" => ""), r"[-/.]" => "_")

"""
    ident_name(x)

Keyword name without namespace: `":variant.impact/high"` -> `"high"`; other
values unchanged.
"""
ident_name(x::AbstractString) = occursin(r"^:[^/]+/", x) ? replace(x, r"^:[^/]+/" => "") : x
ident_name(x) = x

is_map(x) = x isa AbstractDict && !isempty(x)
is_ident_map(x) = x isa AbstractDict && collect(keys(x)) == [":db/ident"]
is_eid_map(x) = x isa AbstractDict && collect(keys(x)) == [":db/id"]
is_scalar(x) = !(x isa AbstractDict || x isa AbstractVector)

"""
    flatten_pull(m; exclude_ids=true)

Flatten one pulled entity into columns. Scalar attributes become columns named
after the attribute (`:subject/id` -> `subject_id`); enum refs (`{:db/ident ..}`)
become the enum name without namespace under the referring attribute's column;
nested single refs flatten recursively (colliding names are prefixed with the
referring attribute); cardinality-many values stay vectors (a one-element
vector becomes its value).
"""
function flatten_pull(m; exclude_ids::Bool=true)
    out = OrderedDict{String,Any}()
    order = String[]
    put!(k, v) = (haskey(out, k) || push!(order, k); out[k] = v)
    m isa AbstractDict || return out, order
    for (k, v) in m
        (exclude_ids && (k == ":db/id" || endswith(k, "/uid"))) && continue
        v === nothing && continue
        col = clean_name(k)
        if is_ident_map(v)
            put!(col, ident_name(v[":db/ident"]))
        elseif is_eid_map(v)
            exclude_ids || put!(col, v[":db/id"])
        elseif v isa AbstractDict
            nested, norder = flatten_pull(v; exclude_ids=exclude_ids)
            for nk in norder
                put!(haskey(out, nk) ? "$(col)_$(nk)" : nk, nested[nk])
            end
        elseif v isa AbstractVector
            vals = if all(is_ident_map, v)
                Any[ident_name(e[":db/ident"]) for e in v]
            elseif all(e -> e isa AbstractDict && length(e) == 1, v)
                Any[only(values(e)) for e in v]
            elseif all(is_scalar, v)
                Any[e for e in v]
            else
                Any[e isa AbstractDict ? first(flatten_pull(e; exclude_ids=exclude_ids)) : e for e in v]
            end
            put!(col, length(vals) == 1 ? vals[1] : vals)
        else
            put!(col, v)
        end
    end
    out, order
end

is_pull(e) = e isa AbstractVector && !isempty(e) && e[1] == "pull"

function find_column_names(find)
    names = String[]
    for e in find
        push!(names, e isa AbstractVector ?
              join([clean_name(x) for x in e if is_scalar(x)], "_") : clean_name(e))
    end
    seen = Dict{String,Int}()
    out = String[]
    for n in names
        if haskey(seen, n)
            seen[n] += 1
            push!(out, "$(n)_$(seen[n])")
        else
            seen[n] = 0
            push!(out, n)
        end
    end
    out
end

# narrow a column's element type (Any -> Float64, String, ...), nothing -> missing
function tidy_column(v::AbstractVector)
    w = Any[x === nothing ? missing : x for x in v]
    isempty(w) && return w
    identity.(w)
end

"""
    records_to_df(recs, order)

Combine records (Dicts) into a DataFrame, columns in first-seen order, filling
absent fields with `missing`.
"""
function records_to_df(recs::Vector, order::Vector{String}=String[])
    cols = copy(order)
    seen = Set(cols)
    for r in recs, k in keys(r)
        k in seen || (push!(cols, k); push!(seen, k))
    end
    df = DataFrame()
    for c in cols
        df[!, c] = tidy_column(Any[get(r, c, missing) for r in recs])
    end
    df
end

function result_to_df(rows, find; exclude_ids::Bool=true)
    names = find_column_names(find)
    pulls = map(is_pull, find)
    if !any(pulls)
        df = DataFrame()
        for (j, n) in enumerate(names)
            df[!, n] = tidy_column(Any[r[j] for r in rows])
        end
        return df
    end
    recs = Dict{String,Any}[]
    order = String[]
    for r in rows
        rec = Dict{String,Any}()
        for (j, v) in enumerate(r)
            if pulls[j]
                flat, forder = flatten_pull(v; exclude_ids=exclude_ids)
                for k in forder
                    k in order || push!(order, k)
                    rec[k] = flat[k]
                end
            else
                names[j] in order || push!(order, names[j])
                rec[names[j]] = v
            end
        end
        push!(recs, rec)
    end
    records_to_df(recs, order)
end

"""
    join_many(v; sep="; ")

Collapse a cardinality-many value (vector) to a string; scalars unchanged.
"""
join_many(v::AbstractVector; sep="; ") = isempty(v) ? missing : join(string.(v), sep)
join_many(v; sep="; ") = v

"""
    order_columns(df, first)

Known columns first, in the given order; others after. Keeps provenance.
"""
function order_columns(df::AbstractDataFrame, first::Vector{String})
    f = [c for c in first if c in names(df)]
    keep_provenance(select(df, f, Not(f)), df)
end

# -- provenance ----------------------------------------------------------------

"""
    provenance(x)

Provenance of a result returned by patternq: a NamedTuple `(db, basis_t,
timestamp)`, or `nothing`. Stored as DataFrame table metadata
("patternq_provenance"), so it follows copies and most transformations.
"""
function provenance(x::AbstractDataFrame)
    "patternq_provenance" in metadatakeys(x) ? metadata(x, "patternq_provenance") : nothing
end
provenance(x::AbstractDict) = get(x, "patternq_provenance", nothing)
provenance(x) = nothing

with_provenance(df::AbstractDataFrame, db, basis_t) =
    (metadata!(df, "patternq_provenance",
               (db=db, basis_t=basis_t, timestamp=Dates.format(now(), "yyyy-mm-dd HH:MM:SS"));
               style=:note); df)

"""
    keep_provenance(x, from)

Carry provenance over to a table derived from a patternq result (joins drop it).
"""
function keep_provenance(x::AbstractDataFrame, from)
    p = provenance(from)
    p === nothing || metadata!(x, "patternq_provenance", p; style=:note)
    x
end
keep_provenance(x, from) = x
