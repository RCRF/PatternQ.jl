# Measurement-set level queries: which measurement sets exist, what they
# measure, and generic measurement retrieval.

# A query that binds every measurement of a set (count, sample, [?m ?a]) makes
# the server hold the whole set in memory before any aggregate runs; on large
# expression sets that runs the service out of heap. The functions below read
# the set's measurement ids from the datoms index instead, in bounded pages,
# and query only those ids.
const PAGE = 10000  # the server's datoms page cap; each page is also an audit write

function measurement_set_eids(measurement_set, db; kwargs...)
    r = runq(dq(["?ms"], [["?ms", ":measurement-set/name", "?ms-name"]]; in=["?ms-name"],
                args=[measurement_set]); db=db, kwargs...)
    nrow(r) == 0 && error("No measurement set named '$(measurement_set)' in $(db)")
    Int.(r.ms)
end

# Up to `limit` measurement ids of set `ms` in index order, starting at
# `offset`, and whether they reach the end of the set.
function set_measurement_ids(ms, db, limit; offset=0)
    ids = Int[]
    while length(ids) < limit
        want = min(PAGE, limit - length(ids))
        page = datoms("eavt", Any[ms, ":measurement-set/measurements"]; db=db,
                      offset=offset + length(ids), limit=want)
        append!(ids, Int.(page.v))
        nrow(page) < want && return ids, true
    end
    nxt = datoms("eavt", Any[ms, ":measurement-set/measurements"]; db=db, offset=offset + length(ids), limit=1)
    ids, nrow(nxt) == 0
end

# (count, exact): the set's measurement count when it is at most `cap`, else
# (cap, false). One page for a small set. A set filling the first page gets one
# probe past the cap, then one-datom probes binary-search the 10000-wide block
# it ends in, and one page of that block finishes the count: about six round
# trips for any size, instead of a page per 10000.
function count_set_measurements(ms, db, cap)
    page(offset, limit) = nrow(datoms("eavt", Any[ms, ":measurement-set/measurements"]; db=db,
                                      offset=offset, limit=limit))
    step = 10000
    n = page(0, step)
    n < step && return n, true
    page(cap, 1) > 0 && return cap, false
    lo, hi = 1, cld(cap, step)
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        page(mid * step, 1) > 0 ? (lo = mid) : (hi = mid)
    end
    lo * step + page(lo * step, step), true
end

"""
    measurement_sets(; db=nothing, count_cap=100000)

assay_name, assay_technology, measurement_set_name, measurement_count,
measurement_count_exact. Counting never binds a whole set in one query: a set
with more than `count_cap` measurements reports `count_cap` with
measurement_count_exact = false (read it as "at least").
"""
function measurement_sets(; db=nothing, count_cap=100000, kwargs...)
    db = ensure_db(db)
    df = dataset_summary(; db=db, kwargs...)
    sets = runq(dq(["?ms", "?measurement-set-name"],
                   [["?a", ":assay/measurement-sets", "?ms"],
                    ["?ms", ":measurement-set/name", "?measurement-set-name"]]); db=db, kwargs...)
    counts = Dict{String,Tuple{Int,Bool}}()
    # one or two round trips per set, and a dataset can have dozens of sets:
    # count them concurrently (the back pressure still caps requests in flight)
    found = asyncmap(ms -> count_set_measurements(Int(ms), db, count_cap), sets.ms; ntasks=4)
    for (row, (n, exact)) in zip(eachrow(sets), found)
        prev = get(counts, row.measurement_set_name, (0, true))
        counts[row.measurement_set_name] = (prev[1] + n, prev[2] && exact)
    end
    df.measurement_count = [get(counts, n, (0, true))[1] for n in df.measurement_set_name]
    df.measurement_count_exact = [get(counts, n, (0, true))[2] for n in df.measurement_set_name]
    df
end

"""
    measurement_types(measurement_set; db=nothing, limit=20000)

Counts of every attribute across a measurement set's measurements: value
attributes (tpm, percent-of-parent, ...) and target references (gene-product,
cell-population, ...). A set can mix several kinds of measurement.
Columns: attribute, kind ("value"/"target"), count.

Counts cover the set's first `limit` measurements in index order; the result's
metadata "exact" is false when the set is larger, and the counts are then a
sample (kinds laid out after the first `limit` can be missed).
"""
function measurement_types(measurement_set::AbstractString; db=nothing, limit=20000, kwargs...)
    db = ensure_db(db)
    parts = DataFrame[]
    exact = true
    for ms in measurement_set_eids(measurement_set, db; kwargs...)
        ids, complete = set_measurement_ids(ms, db, limit)
        exact &= complete
        for chunk in Iterators.partition(ids, PAGE)
            r = runq(dq(Any["?attribute", ["count", "?m"]], [["?m", "?a"], ["?a", ":db/ident", "?attribute"]];
                        in=[["?m", "..."]], args=[collect(chunk)]); db=db, kwargs...)
            nrow(r) > 0 && push!(parts, r)
        end
    end
    if isempty(parts)
        out = DataFrame(attribute=String[], kind=String[], count=Int[])
        metadata!(out, "exact", exact; style=:note)
        return out
    end
    r = combine(groupby(vcat(parts...), :attribute), :count_m => sum => :count)
    r = filter(row -> !(row.attribute in (":measurement/id", ":measurement/uid", ":measurement/sample")), r)
    r.kind = [haskey(MEASUREMENT_TARGETS, a) ? "target" : "value" for a in r.attribute]
    r.attribute = [replace(a, r"^:measurement/" => "") for a in r.attribute]
    out = keep_provenance(select(sort(r, [:kind, order(:count, rev=true)]), :attribute, :kind, :count), parts[1])
    metadata!(out, "exact", exact; style=:note)
    out
end

"""
    measurement_set_attributes(measurement_set; db=nothing, measurement=nothing, n=200)

Attributes of the set's first `n` measurements in index order (optionally only
those carrying `measurement`, looked for among the set's first 100000
measurements, the deepest the server pages by offset). Reads bounded pages of
the datoms index, so it stays cheap for sets with millions of measurements.
"""
function measurement_set_attributes(measurement_set::AbstractString; db=nothing, measurement=nothing, n=200)
    db = ensure_db(db)
    attr = measurement === nothing ? nothing : measurement_attr(measurement)
    found = Int[]
    for ms in measurement_set_eids(measurement_set, db)
        offset = 0
        for _ in 1:10  # the server caps datoms offsets at 100000
            ids, complete = set_measurement_ids(ms, db, attr === nothing ? n : PAGE; offset=offset)
            offset += length(ids)
            if attr !== nothing && !isempty(ids)
                r = runq(dq(["?m"], [["?m", attr]]; in=[["?m", "..."]], args=[ids]); db=db)
                ids = nrow(r) > 0 ? Int.(r.m) : Int[]
            end
            append!(found, ids)
            (length(found) >= n || complete || attr === nothing) && break
        end
        length(found) >= n && break
    end
    found = first(found, n)
    isempty(found) && return String[]
    attrs = runq(dq(["?attr"], [["?m", "?a"], ["?a", ":db/ident", "?attr"]];
                    in=[["?m", "..."]], args=[found]); db=db)
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

sample_assays_query() = dq(["?sample-id", "?assay-name", "?assay-technology", "?measurement-set-name"],
                           [["?a", ":assay/name", "?assay-name"],
                            ["?a", ":assay/technology", "?t"],
                            ["?t", ":db/ident", "?assay-technology"],
                            ["?a", ":assay/measurement-sets", "?ms"],
                            ["?ms", ":measurement-set/name", "?measurement-set-name"],
                            ["?ms", ":measurement-set/includes-samples", "?s"],
                            ["?s", ":sample/id", "?sample-id"]])

has_includes_samples(db; kwargs...) =
    nrow(runq(dq(["?a"], [["?a", ":db/ident", "?ident"]]; in=["?ident"],
                 args=[":measurement-set/includes-samples"]); db=db, kwargs...)) > 0

"""
    sample_assays(; db=nothing)

Which samples were measured by which measurement sets: subject_id, sample_id,
timepoint_id, assay_name, assay_technology, measurement_set_name.

Reads `:measurement-set/includes-samples`, the enrichment that lists each set's
samples (schema 0.3.18), so nothing walks the measurements. A database from
before that enrichment falls back to walking each set's measurements, which can
time out on large expression sets.
"""
function sample_assays(; db=nothing, kwargs...)
    db = ensure_db(db)
    if has_includes_samples(db; kwargs...)
        df = runq(sample_assays_query(); db=db, kwargs...)
        nrow(df) == 0 && return DataFrame()
        df.assay_technology = ident_name.(df.assay_technology)
        smp = samples(; db=db)
        keep = intersect(["sample_id", "subject_id", "timepoint_id"], names(smp))
        out = rightjoin(select(smp, keep), df; on=:sample_id)
        return keep_provenance(order_columns(out, ["subject_id", "sample_id", "timepoint_id"]), df)
    end
    # older database: walk each set's measurement -> sample refs
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
        # Datomic joins in clause order: start from the genes and the segments
        # covering them, not from every CNV measurement in the database
        prepend!(where, [["?g", ":gene/hgnc-symbol", "?hgnc-symbol"], ["?c", ":cnv/genes", "?g"]])
        push!(find, "?hgnc-symbol"); push!(ins, ["?hgnc-symbol", "..."]); push!(args, collect(genes))
    end
    samples === nothing || (push!(ins, ["?sample-id", "..."]); push!(args, collect(samples)))
    add_subject_filter!(where, ins, args, subjects)
    genes === nothing && samples !== nothing && (where = lead_with(where, SAMPLE_LEAD))
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
    if genes === nothing
        where = lead_with(where, samples !== nothing ? SAMPLE_LEAD : Any[["?m", attr, enum ? "?value-ref" : "?value"]])
    end
    df = runq(dq(["?sample-id", "?measurement-set", "?hgnc-symbol", "?value"], where; in=ins, args=args,
                 with=["?m"]); db=db, kwargs...)
    eltype(df.value) <: Union{Missing,AbstractString} && (df.value = ident_name.(df.value))
    df
end
