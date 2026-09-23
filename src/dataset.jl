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
    dq(["?sample-id", "?hgnc-symbol", "?measurement-set", "?value"], where; in=ins, args=args, with=["?m"])
end

"""
    gene_expression(; db=nothing, genes=nothing, samples=nothing, measurement="tpm", measurement_set=nothing)

Gene expression, long format: sample_id, hgnc_symbol, measurement_set, value.
`measurement` is the attribute without namespace ("tpm",
"rsem-normalized-count", ...). Datasets may have several measurement sets
carrying the attribute (bulk, pseudobulk, ...); `measurement_set` tells them apart.
"""
gene_expression(; db=nothing, genes=nothing, samples=nothing, measurement="tpm", measurement_set=nothing, kwargs...) =
    runq(gene_expression_query(genes=genes, samples=samples, measurement=measurement,
                              measurement_set=measurement_set); db=db, kwargs...)
