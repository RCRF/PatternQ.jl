# Canned queries over a dataset database.
#
# Every dataset is its own database, so queries start from samples, subjects,
# assays and measurement sets directly; there is no dataset argument. Each
# function has a *_query companion returning (query, args) as a NamedTuple, so
# the query can be inspected, modified or run elsewhere.

"""
    dq(find, where; in=nothing, args=[], with=nothing)

Build a query (Dict in the service's JSON form) and its arguments:
`(query=Dict(...), args=[...])`.
"""
function dq(find, where; in=nothing, args=Any[], with=nothing)
    q = Dict{String,Any}(":find" => collect(Any, find), ":where" => collect(Any, where))
    in === nothing || isempty(in) || (q[":in"] = collect(Any, in))
    with === nothing || (q[":with"] = collect(Any, with))
    (query=q, args=collect(Any, args))
end

runq(qa; db=nothing, kwargs...) = do_query(qa.query; args=qa.args, db=db, kwargs...)

const EnumRef = (attr) -> Dict(attr => [":db/ident"])

"""
    samples_query()
"""
samples_query() = dq([["pull", "?s", Any["*",
                        Dict(":sample/subject" => [":subject/id"]),
                        Dict(":sample/timepoint" => [":timepoint/id", ":timepoint/relative-order"]),
                        EnumRef(":sample/specimen"), EnumRef(":sample/type"), EnumRef(":sample/container"),
                        Dict(":sample/study-day" => [":study-day/id"]),
                        Dict(":sample/gdc-anatomic-site" => [":gdc-anatomic-site/name"])]]],
                     [["?s", ":sample/id"]])

"""
    samples(; db=nothing)

One row per sample: sample_id, subject_id, timepoint_id and the other sample
attributes present.
"""
samples(; db=nothing, kwargs...) = runq(samples_query(); db=db, kwargs...)

subjects_query() = dq([["pull", "?s", Any["*",
                         EnumRef(":subject/sex"), EnumRef(":subject/race"), EnumRef(":subject/ethnicity"),
                         EnumRef(":subject/smoker"), EnumRef(":subject/cause-of-death"),
                         Dict(":subject/meddra-disease" => [":meddra-disease/preferred-name"])]]],
                      [["?s", ":subject/id"]])

"""
    subjects(; db=nothing)

One row per subject: subject_id and demographic attributes present (enums as
names, e.g. subject_sex = "female").
"""
subjects(; db=nothing, kwargs...) = runq(subjects_query(); db=db, kwargs...)

dataset_summary_query() = dq(["?assay-name", "?assay-technology", "?measurement-set-name"],
                             [["?a", ":assay/name", "?assay-name"],
                              ["?a", ":assay/technology", "?t"],
                              ["?t", ":db/ident", "?assay-technology"],
                              ["?a", ":assay/measurement-sets", "?ms"],
                              ["?ms", ":measurement-set/name", "?measurement-set-name"]])

"""
    dataset_summary(; db=nothing)

One row per measurement set: assay_name, assay_technology, measurement_set_name.
"""
function dataset_summary(; db=nothing, kwargs...)
    df = runq(dataset_summary_query(); db=db, kwargs...)
    df.assay_technology = ident_name.(df.assay_technology)
    keep_provenance(sort(df, [:assay_name, :measurement_set_name]), df)
end

"""
    dataset_info(; db=nothing)

The dataset entity stored in the database: name, description, doi, url.
"""
dataset_info(; db=nothing, kwargs...) =
    runq(dq([["pull", "?d", [":dataset/name", ":dataset/description", ":dataset/doi", ":dataset/url"]]],
           [["?d", ":dataset/name"]]); db=db, kwargs...)

"""
    schema_info(; db=nothing)

Schema name and version of a database: `(name=..., version=...)`.
"""
function schema_info(; db=nothing, kwargs...)
    r = runq(dq(["?name", "?version"], [["?e", ":unify.schema/version", "?version"],
                                        ["?e", ":unify.schema/name", "?name"]]); db=db, kwargs...)
    (name=r.name[1], version=r.version[1])
end

# Datomic joins in clause order. lead_with reorders `where` to start from
# `lead`, the clauses of the subset the caller bound (genes, samples), then each
# remaining clause as soon as it shares a variable with those placed, so a query
# never starts on a whole attribute (every variant, CNV, or expression
# measurement) and filters at the end. Starting from a few samples walks only
# their measurements; a subject can own most of a dataset (a single-participant
# deep dive), so subject-only queries start from the measurement attribute.
const SAMPLE_LEAD = Any[["?s", ":sample/id", "?sample-id"], ["?m", ":measurement/sample", "?s"]]

# Without genes, a sample-subset call returns whole samples' worth of rows
# (every gene's expression or CNV call): one /query per this many samples, so
# no single query holds more on the server.
const SAMPLES_PER_QUERY = 2

function by_sample_chunks(f, samples)
    parts = [f(collect(c)) for c in Iterators.partition(collect(samples), SAMPLES_PER_QUERY)]
    rows = filter(p -> nrow(p) > 0, parts)
    keep_provenance(isempty(rows) ? parts[1] : vcat(rows...; cols=:union), parts[1])
end

# The subjects' sample ids, which subject-subset calls run by.
function subject_samples(subjects, db; kwargs...)
    r = runq(dq(["?sample-id"], [["?p", ":subject/id", "?subject-id"], ["?s", ":sample/subject", "?p"],
                                 ["?s", ":sample/id", "?sample-id"]];
                in=[["?subject-id", "..."]], args=[collect(subjects)]); db=db, kwargs...)
    nrow(r) > 0 ? sort(String.(r.sample_id)) : String[]
end

# A set's first measurement id, cached (a database name pins an immutable
# database).
const FIRST_MEASUREMENT = Dict{Tuple{String,Int},Union{Int,Nothing}}()
const FIRST_MEASUREMENT_LOCK = ReentrantLock()
function first_measurement(db, ms)
    key = (String(db), Int(ms))
    hit = lock(() -> get(FIRST_MEASUREMENT, key, missing), FIRST_MEASUREMENT_LOCK)
    hit === missing || return hit
    d = datoms("eavt", Any[ms, ":measurement-set/measurements"]; db=db, limit=1)
    v = nrow(d) > 0 ? Int(d.v[1]) : nothing
    lock(() -> (FIRST_MEASUREMENT[key] = v), FIRST_MEASUREMENT_LOCK)
    v
end

# The sets in `sets` holding `attr` measurements, judged by each set's first
# measurement (a set holds one kind of measurement): one datoms page per set
# (cached, 4 at a time), then one query over those measurements.
function sets_carrying(attr, sets, db)
    firsts = asyncmap(ms -> first_measurement(db, ms), sets; ntasks=4)
    ids = Int[f for f in firsts if f !== nothing]
    isempty(ids) && return Int[]
    r = runq(dq(["?m"], [["?m", attr]]; in=[["?m", "..."]], args=[ids]); db=db)
    hit = nrow(r) > 0 ? Set(Int.(r.m)) : Set{Int}()
    Int[ms for (ms, f) in zip(sets, firsts) if f !== nothing && f in hit]
end

const SET_WALK_CAP = 2000000  # walking more set measurements than this is never the cheaper start

# How a CNV subset query should run: (:empty, []) when no `attr` set includes
# the samples; (:sets, sets) to walk just those sets' measurements; (:samples,
# sets) to walk the samples' measurements, restricted to those sets.
# :measurement-set/includes-samples (schema 0.3.18) names the sets holding each
# sample; the server counts both walks in place, and the smaller one wins. A
# database without the enrichment walks the samples.
function cnv_plan(attr, sample_ids, db; kwargs...)
    has_includes_samples(db; kwargs...) || return :samples, Int[]
    r = runq(dq(["?s", "?ms"], [["?s", ":sample/id", "?sample-id"], ["?ms", ":measurement-set/includes-samples", "?s"]];
                in=[["?sample-id", "..."]], args=[collect(sample_ids)]); db=db, kwargs...)
    sets = nrow(r) > 0 ? sets_carrying(attr, sort(unique(Int.(r.ms))), db) : Int[]
    isempty(sets) && return :empty, Int[]
    in_sets = 0
    for ms in sets
        n = datoms_count("eavt", Any[ms, ":measurement-set/measurements"]; db=db, cap=SET_WALK_CAP)
        (n === nothing || !n[2]) && return :samples, sets
        in_sets += n[1]
    end
    walked = 0
    for e in sort(unique(Int.(r.s)))
        n = datoms_count("vaet", Any[e, ":measurement/sample"]; db=db, cap=in_sets + 1)
        walked += n === nothing ? in_sets + 1 : n[1]
        walked > in_sets && return :sets, sets
    end
    :samples, sets
end

# Run a CNV subset call by plan: run(samples, start, sets) is the query.
function run_cnv_plan(run, attr, samples, subjects, db; kwargs...)
    ids = samples !== nothing ? collect(samples) : subject_samples(subjects, db; kwargs...)
    isempty(ids) && return DataFrame()
    start, sets = cnv_plan(attr, ids, db; kwargs...)
    start === :empty && return DataFrame()
    start === :samples && length(ids) > SAMPLES_PER_QUERY &&
        return by_sample_chunks(c -> run(c, start, sets), ids)
    run(ids, start, sets)
end

function lead_with(where, lead)
    vars_of(c) = Set(x for x in c if x isa AbstractString && startswith(x, "?"))
    out = Any[lead...]
    bound = isempty(out) ? Set{String}() : union(map(vars_of, out)...)
    rest = Any[c for c in where if !(c in lead)]
    while !isempty(rest)
        i = something(findfirst(c -> !isempty(intersect(vars_of(c), bound)), rest), 1)
        c = popat!(rest, i)
        push!(out, c)
        union!(bound, vars_of(c))
    end
    out
end

# -- variants --------------------------------------------------------------------

function variants_query(; samples=nothing, genes=nothing, measurement_set=nothing)
    where = Any[["?m", ":measurement/vaf", "?vaf"], ["?m", ":measurement/variant", "?v"],
                ["?m", ":measurement/sample", "?s"], ["?s", ":sample/id", "?sample-id"],
                ["?ms", ":measurement-set/measurements", "?m"],
                ["?ms", ":measurement-set/name", "?measurement-set"]]
    ins = Any[]; args = Any[]
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    if genes !== nothing
        append!(where, [["?v", ":variant/gene", "?g"], ["?g", ":gene/hgnc-symbol", "?gene"]])
        push!(ins, ["?gene", "..."]); push!(args, collect(genes))
    end
    measurement_set === nothing || (push!(ins, "?measurement-set"); push!(args, measurement_set))
    if genes !== nothing
        where = lead_with(where, Any[["?g", ":gene/hgnc-symbol", "?gene"], ["?v", ":variant/gene", "?g"],
                                     ["?m", ":measurement/variant", "?v"]])
    elseif samples !== nothing
        where = lead_with(where, SAMPLE_LEAD)
    end
    dq(Any["?sample-id", "?measurement-set", "?vaf",
           ["pull", "?v", Any[":variant/id", ":variant/HGVSp", ":variant/HGVSc",
                              Dict(":variant/gene" => [":gene/hgnc-symbol"]), EnumRef(":variant/impact")]],
           ["pull", "?m", [":measurement/t-depth"]]],
       where; in=ins, args=args)
end

"""
    variants(; db=nothing, samples=nothing, genes=nothing, measurement_set=nothing)

Somatic variant measurements, one row per measurement: sample_id,
measurement_set, variant_id, hgnc_symbol, HGVSp, HGVSc, impact, vaf, t_depth
(where present). HGVSp/HGVSc (cardinality many) are joined with "; ".
"""
function variants(; db=nothing, samples=nothing, genes=nothing, measurement_set=nothing, kwargs...)
    if genes === nothing && samples !== nothing && length(samples) > SAMPLES_PER_QUERY
        return by_sample_chunks(c -> variants(; db=db, samples=c, measurement_set=measurement_set, kwargs...), samples)
    end
    df = runq(variants_query(samples=samples, genes=genes, measurement_set=measurement_set); db=db, kwargs...)
    rename_cols!(df, "variant_HGVSp" => "HGVSp", "variant_HGVSc" => "HGVSc", "variant_impact" => "impact",
                 "gene_hgnc_symbol" => "hgnc_symbol", "measurement_t_depth" => "t_depth")
    for c in ("HGVSp", "HGVSc")
        c in names(df) && (df[!, c] = tidy_column(Any[join_many(v) for v in df[!, c]]))
    end
    order_columns(df, ["sample_id", "measurement_set", "variant_id", "hgnc_symbol", "HGVSp", "HGVSc",
                       "impact", "vaf", "t_depth"])
end

rename_cols!(df, pairs...) = (for (a, b) in pairs; a in names(df) && rename!(df, a => b); end; df)

# -- gene expression -------------------------------------------------------------

measurement_attr(m) = ":measurement/" * replace(lstrip(string(m), ':'), r"^measurement/" => "")

function gene_expression_query(; genes=nothing, samples=nothing, measurement="tpm", measurement_set=nothing)
    where = Any[["?g", ":gene/hgnc-symbol", "?hgnc-symbol"], ["?gp", ":gene-product/gene", "?g"],
                ["?m", ":measurement/gene-product", "?gp"], ["?m", measurement_attr(measurement), "?value"],
                ["?m", ":measurement/sample", "?s"], ["?s", ":sample/id", "?sample-id"],
                ["?ms", ":measurement-set/measurements", "?m"], ["?ms", ":measurement-set/name", "?measurement-set"]]
    ins = Any[]; args = Any[]
    genes === nothing || (push!(ins, ["?hgnc-symbol", "..."]); push!(args, collect(genes)))
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    measurement_set === nothing || (push!(ins, "?measurement-set"); push!(args, measurement_set))
    genes === nothing && samples !== nothing && (where = lead_with(where, SAMPLE_LEAD))
    dq(["?sample-id", "?hgnc-symbol", "?measurement-set", "?value"], where; in=ins, args=args, with=["?m"])
end

"""
    gene_expression(; db=nothing, genes=nothing, samples=nothing, measurement="tpm", measurement_set=nothing)

Gene expression, long format: sample_id, hgnc_symbol, measurement_set, value.
`measurement` is the attribute without namespace ("tpm",
"rsem-normalized-count", ...). Datasets may have several measurement sets
carrying the attribute (bulk, pseudobulk, ...); `measurement_set` tells them apart.
"""
function gene_expression(; db=nothing, genes=nothing, samples=nothing, measurement="tpm", measurement_set=nothing,
                         kwargs...)
    if genes === nothing && samples !== nothing && length(samples) > SAMPLES_PER_QUERY
        return by_sample_chunks(c -> gene_expression(; db=db, samples=c, measurement=measurement,
                                                     measurement_set=measurement_set, kwargs...), samples)
    end
    runq(gene_expression_query(genes=genes, samples=samples, measurement=measurement,
                               measurement_set=measurement_set); db=db, kwargs...)
end
