# Measurement-set level queries: which measurement sets exist, what they
# measure, and generic measurement retrieval.

"""
    measurement_sets(; db=nothing)

assay_name, assay_technology, measurement_set_name, measurement_count.
"""
function measurement_sets(; db=nothing, kwargs...)
    df = dataset_summary(; db=db, kwargs...)
    counts = runq(dq(Any["?measurement-set-name", ["count", "?m"]],
                     [["?ms", ":measurement-set/name", "?measurement-set-name"],
                      ["?ms", ":measurement-set/measurements", "?m"]]); db=db, kwargs...)
    lookup = Dict(zip(counts.measurement_set_name, counts.count_m))
    df.measurement_count = [get(lookup, n, 0) for n in df.measurement_set_name]
    df
end

"""
    measurement_types(measurement_set; db=nothing)

Counts of every attribute across a measurement set's measurements: value
attributes (tpm, percent-of-parent, ...) and target references (gene-product,
cell-population, ...). A set can mix several kinds of measurement.
Columns: attribute, kind ("value"/"target"), count.
"""
function measurement_types(measurement_set::AbstractString; db=nothing, kwargs...)
    r = runq(dq(Any["?attribute", ["count", "?m"]],
                [["?ms", ":measurement-set/name", "?ms-name"], ["?ms", ":measurement-set/measurements", "?m"],
                 ["?m", "?a"], ["?a", ":db/ident", "?attribute"]];
                in=["?ms-name"], args=[measurement_set]); db=db, timeout=120, kwargs...)
    r = filter(row -> !(row.attribute in (":measurement/id", ":measurement/uid", ":measurement/sample")), r)
    rename!(r, "count_m" => "count")
    r.kind = [haskey(MEASUREMENT_TARGETS, a) ? "target" : "value" for a in r.attribute]
    r.attribute = [replace(a, r"^:measurement/" => "") for a in r.attribute]
    out = sort(r, [:kind, order(:count, rev=true)])
    keep_provenance(select(out, :attribute, :kind, :count), r)
end

"""
    measurement_set_attributes(measurement_set; db=nothing, measurement=nothing, n=200)

Attributes of a random sample of `n` measurements in the set (optionally only
those carrying `measurement`); cheap even for sets with millions of measurements.
"""
function measurement_set_attributes(measurement_set::AbstractString; db=nothing, measurement=nothing, n=200)
    where = Any[["?ms", ":measurement-set/name", "?ms-name"], ["?ms", ":measurement-set/measurements", "?m"]]
    measurement === nothing || push!(where, ["?m", measurement_attr(measurement)])
    s = query(Dict(":find" => [["sample", n, "?m"]], ":in" => ["?ms-name"], ":where" => where);
              args=[measurement_set], db=db, timeout=120)
    eids = reduce(vcat, [vcat(r[1]...) for r in s["query_result"]]; init=Any[])
    isempty(eids) && return String[]
    attrs = runq(dq(["?attr"], [["?m", "?a"], ["?a", ":db/ident", "?attr"]];
                    in=[["?m", "..."]], args=[eids]); db=db)
    sort(setdiff(attrs.attr, [":measurement/uid", ":measurement/id"]))
end

# Measurement target references: the entity a measurement is "of", how to
# name it, and the result variable.
const MEASUREMENT_TARGETS = Dict(
    ":measurement/gene-product" => (clauses=[["?tgp", ":gene-product/gene", "?tg"], ["?tg", ":gene/hgnc-symbol", "?hgnc-symbol"]],
                                    ref="?tgp", var="?hgnc-symbol"),
    ":measurement/variant" => (clauses=[["?tv", ":variant/id", "?variant-id"]], ref="?tv", var="?variant-id"),
    ":measurement/cnv" => (clauses=[["?tc", ":cnv/id", "?cnv-id"]], ref="?tc", var="?cnv-id"),
    ":measurement/epitope" => (clauses=[["?te", ":epitope/id", "?epitope-id"]], ref="?te", var="?epitope-id"),
    ":measurement/cell-population" => (clauses=[["?tcp", ":cell-population/name", "?cell-population"]],
                                       ref="?tcp", var="?cell-population"),
    ":measurement/tcr" => (clauses=[["?tt", ":tcr/id", "?tcr-id"]], ref="?tt", var="?tcr-id"),
    ":measurement/otu" => (clauses=[["?to", ":otu/id", "?otu-id"]], ref="?to", var="?otu-id"),
    ":measurement/sgb" => (clauses=[["?ts", ":sgb/metaphlan-id", "?sgb-id"]], ref="?ts", var="?sgb-id"),
    ":measurement/pathway" => (clauses=[["?tp", ":pathway/id", "?pathway-id"]], ref="?tp", var="?pathway-id"),
    ":measurement/metabolite-feature" => (clauses=[["?tmf", ":metabolite-feature/rt-mz-peak", "?metabolite-feature"]],
                                          ref="?tmf", var="?metabolite-feature"),
    ":measurement/nanostring-signature" => (clauses=[["?tns", ":nanostring-signature/name", "?signature"]],
                                            ref="?tns", var="?signature"),
    ":measurement/atac-peak" => (clauses=[["?tap", ":atac-peak/name", "?atac-peak"]], ref="?tap", var="?atac-peak"),
    ":measurement/single-cell" => (clauses=[["?tsc", ":single-cell/id", "?single-cell-id"]],
                                   ref="?tsc", var="?single-cell-id"))

const ENUM_MEASUREMENT_ATTRS = Set([":measurement/cnv-call", ":measurement/msi-status"])

function measurements_query(measurement; measurement_set=nothing, samples=nothing, targets=String[])
    attr = measurement_attr(measurement)
    enum = attr in ENUM_MEASUREMENT_ATTRS
    where = Any[["?m", attr, enum ? "?value-ref" : "?value"], ["?m", ":measurement/sample", "?s"],
                ["?s", ":sample/id", "?sample-id"], ["?ms", ":measurement-set/measurements", "?m"],
                ["?ms", ":measurement-set/name", "?measurement-set"]]
    enum && push!(where, ["?value-ref", ":db/ident", "?value"])
    find = Any["?sample-id", "?measurement-set"]
    for t in targets
        spec = get(MEASUREMENT_TARGETS, t, nothing)
        spec === nothing && error("Unknown measurement target $t")
        push!(where, ["?m", t, spec.ref]); append!(where, spec.clauses)
        push!(find, spec.var)
    end
    push!(find, "?value")
    ins = Any[]; args = Any[]
    measurement_set === nothing || (push!(ins, "?measurement-set"); push!(args, measurement_set))
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    dq(find, where; in=ins, args=args, with=["?m"])
end

"""
    measurements(measurement, measurement_set=nothing; db=nothing, samples=nothing, targets=nothing, wide=false, fun=mean)

Values of any measurement attribute (e.g. "tpm", "percent-of-parent",
"olink-npx", "median-channel-value") with what each measurement is of (gene,
cell population, epitope, ...). Targets are detected from a sample of the
measurements carrying the attribute unless given. Long format: sample_id,
measurement_set, target column(s), value; `wide=true` returns a samples x
targets matrix (see `to_matrix`).
"""
function measurements(measurement, measurement_set=nothing; db=nothing, samples=nothing, targets=nothing,
                      wide::Bool=false, fun=mean, kwargs...)
    db = ensure_db(db)
    if targets === nothing
        measurement_set === nothing && error("Give measurement_set (or targets) so measurement targets can be detected")
        targets = [a for a in measurement_set_attributes(measurement_set; db=db, measurement=measurement)
                   if haskey(MEASUREMENT_TARGETS, a)]
    end
    df = runq(measurements_query(measurement; measurement_set=measurement_set, samples=samples,
                                 targets=targets); db=db, kwargs...)
    eltype(df.value) <: Union{Missing,AbstractString} && (df.value = ident_name.(df.value))
    if wide
        tcols = setdiff(names(df), ["sample_id", "measurement_set", "value"])
        isempty(tcols) && (tcols = ["measurement_set"])
        return to_matrix(df, tcols; fun=fun)
    end
    df
end

"""
    sample_assays(; db=nothing)

Which samples were measured by which measurement sets: subject_id, sample_id,
timepoint_id, assay_name, assay_technology, measurement_set_name.
"""
function sample_assays(; db=nothing, kwargs...)
    db = ensure_db(db)
    sets = dataset_summary(; db=db)
    parts = DataFrame[]
    for row in eachrow(sets)
        r = runq(dq(["?sample-id"], [["?ms", ":measurement-set/name", "?ms-name"],
                                      ["?ms", ":measurement-set/measurements", "?m"],
                                      ["?m", ":measurement/sample", "?s"], ["?s", ":sample/id", "?sample-id"]];
                    in=["?ms-name"], args=[row.measurement_set_name]); db=db, kwargs...)
        nrow(r) == 0 && continue
        push!(parts, DataFrame(sample_id=r.sample_id, assay_name=row.assay_name,
                               assay_technology=row.assay_technology,
                               measurement_set_name=row.measurement_set_name))
    end
    isempty(parts) && return DataFrame()
    df = vcat(parts...)
    smp = samples(; db=db)
    keep = intersect(["sample_id", "subject_id", "timepoint_id"], names(smp))
    out = rightjoin(select(smp, keep), df; on=:sample_id)
    keep_provenance(order_columns(out, ["subject_id", "sample_id", "timepoint_id"]), smp)
end

"""
    measurement_matrices(; db=nothing)

File-backed measurement matrices: assay_name, measurement_set_name,
matrix_name, measurement_type, matrix_key (pass to `measurement_matrix`).
"""
function measurement_matrices(; db=nothing, kwargs...)
    df = runq(dq(["?assay-name", "?measurement-set-name", "?matrix-name", "?measurement-type", "?matrix-key"],
                 [["?a", ":assay/name", "?assay-name"], ["?a", ":assay/measurement-sets", "?ms"],
                  ["?ms", ":measurement-set/name", "?measurement-set-name"],
                  ["?ms", ":measurement-set/measurement-matrices", "?mm"],
                  ["?mm", ":measurement-matrix/name", "?matrix-name"],
                  ["?mm", ":measurement-matrix/measurement-type", "?mt"], ["?mt", ":db/ident", "?measurement-type"],
                  ["?mm", ":measurement-matrix/backing-file", "?matrix-key"]]); db=db, kwargs...)
    nrow(df) > 0 && (df.measurement_type = ident_name.(df.measurement_type))
    df
end

"""
    measurement_matrix_by_name(matrix_name; db=nothing)
"""
function measurement_matrix_by_name(matrix_name::AbstractString; db=nothing)
    mm = measurement_matrices(; db=db)
    i = findfirst(==(matrix_name), mm.matrix_name)
    i === nothing && error("No measurement matrix named '$matrix_name'")
    measurement_matrix(mm.matrix_key[i]; db=db)
end

"""
    isoforms(gene; db=nothing, samples=nothing)

Isoform-level expression: sample_id, transcript_id, transcript_length,
isoform_percent, effective_length.
"""
function isoforms(gene::AbstractString; db=nothing, samples=nothing, kwargs...)
    ins = Any["?hgnc"]; args = Any[gene]
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    runq(dq(["?sample-id", "?transcript-id", "?transcript-length", "?isoform-percent", "?effective-length"],
            [["?g", ":gene/hgnc-symbol", "?hgnc"], ["?gp", ":gene-product/gene", "?g"],
             ["?gp", ":gene-product/id", "?transcript-id"], ["?gp", ":gene-product/transcript-length", "?transcript-length"],
             ["?m", ":measurement/gene-product", "?gp"], ["?m", ":measurement/isoform-percent", "?isoform-percent"],
             ["?m", ":measurement/effective-transcript-length", "?effective-length"],
             ["?m", ":measurement/sample", "?s"], ["?s", ":sample/id", "?sample-id"]];
            in=ins, args=args, with=["?m"]); db=db, kwargs...)
end

# CNV data is large (segments x samples, or genes x samples): the CNV queries
# insist on a subset. Same rule in the R, Python and Clojure libraries.
function require_cnv_subset(genes, samples, subjects)
    genes === nothing && samples === nothing && subjects === nothing &&
        throw(ArgumentError("CNV queries need a subset: give genes, samples or subjects"))
end

function add_subject_filter!(where, ins, args, subjects)
    subjects === nothing && return
    append!(where, [["?s", ":sample/subject", "?p"], ["?p", ":subject/id", "?subject-id"]])
    push!(ins, ["?subject-id", "..."]); push!(args, collect(subjects))
end

"""
    cnv_segments(; db=nothing, genes=nothing, samples=nothing, subjects=nothing)

Segment-level copy number measurements (e.g. CNVkit / ASCAT calls) with the
segment coordinates; with `genes`, one row per overlapping segment and gene.
At least one of genes, samples, subjects is required.
"""
function cnv_segments(; db=nothing, genes=nothing, samples=nothing, subjects=nothing, kwargs...)
    require_cnv_subset(genes, samples, subjects)
    where = Any[["?m", ":measurement/cnv", "?c"], ["?c", ":cnv/id", "?cnv-id"], ["?m", ":measurement/sample", "?s"],
                ["?s", ":sample/id", "?sample-id"], ["?ms", ":measurement-set/measurements", "?m"],
                ["?ms", ":measurement-set/name", "?measurement-set"]]
    find = Any["?sample-id", "?measurement-set", "?cnv-id",
               ["pull", "?c", [Dict(":cnv/genomic-coordinates" => [":genomic-coordinate/contig",
                                                                  ":genomic-coordinate/start",
                                                                  ":genomic-coordinate/end"])]],
               ["pull", "?m", [":measurement/segment-mean-lrr", ":measurement/absolute-cn",
                               ":measurement/a-allele-cn", ":measurement/b-allele-cn", ":measurement/loh"]]]
    ins = Any[]; args = Any[]
    if genes !== nothing
        append!(where, [["?c", ":cnv/genes", "?g"], ["?g", ":gene/hgnc-symbol", "?hgnc-symbol"]])
        push!(find, "?hgnc-symbol"); push!(ins, ["?hgnc-symbol", "..."]); push!(args, collect(genes))
    end
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    add_subject_filter!(where, ins, args, subjects)
    df = runq(dq(find, where; in=ins, args=args); db=db, kwargs...)
    for n in names(df)
        m = match(r"^(genomic_coordinate|measurement)_(?!set$)(.*)", n)
        m === nothing || rename!(df, n => m.captures[2])
    end
    order_columns(df, ["sample_id", "measurement_set", "cnv_id", "hgnc_symbol", "contig", "start", "end"])
end

"""
    cnv_gene_calls(; db=nothing, genes=nothing, samples=nothing, subjects=nothing, measurement="cnv-call-score")

Gene-level copy number calls (e.g. GISTIC2 discrete calls): sample_id,
measurement_set, hgnc_symbol, value. At least one of genes, samples, subjects
is required.
"""
function cnv_gene_calls(; db=nothing, genes=nothing, samples=nothing, subjects=nothing,
                        measurement="cnv-call-score", kwargs...)
    require_cnv_subset(genes, samples, subjects)
    attr = measurement_attr(measurement)
    enum = attr in ENUM_MEASUREMENT_ATTRS
    # gene-first clause order: gene-level CNV sets are large (genes x samples)
    where = Any[["?g", ":gene/hgnc-symbol", "?hgnc-symbol"], ["?gp", ":gene-product/gene", "?g"],
                ["?m", ":measurement/gene-product", "?gp"], ["?m", attr, enum ? "?value-ref" : "?value"],
                ["?m", ":measurement/sample", "?s"], ["?s", ":sample/id", "?sample-id"],
                ["?ms", ":measurement-set/measurements", "?m"], ["?ms", ":measurement-set/name", "?measurement-set"]]
    enum && push!(where, ["?value-ref", ":db/ident", "?value"])
    ins = Any[]; args = Any[]
    genes === nothing || (push!(ins, ["?hgnc-symbol", "..."]); push!(args, collect(genes)))
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    add_subject_filter!(where, ins, args, subjects)
    df = runq(dq(["?sample-id", "?measurement-set", "?hgnc-symbol", "?value"], where; in=ins, args=args,
                 with=["?m"]); db=db, kwargs...)
    eltype(df.value) <: Union{Missing,AbstractString} && (df.value = ident_name.(df.value))
    df
end
