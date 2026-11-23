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

# Up to `limit` measurement ids of set `ms` from `pos`, and the position after
# them, or nothing at the end of the set. A position is (offset, cursor): the
# cursor when the server pages by cursor, else the offset. A cursor page
# answered without one (an instance from before cursor paging, during a deploy
# or a rollback) is read again by offset, never repeated.
function set_page(ms, db, limit, pos)
    offset, after = pos isa Tuple ? pos : (pos, nothing)
    cursor(p) = "cursor" in metadatakeys(p) && metadata(p, "cursor")
    page = datoms("eavt", Any[ms, ":measurement-set/measurements"]; db=db,
                  offset=after === nothing ? offset : 0, limit=limit, after=after)
    if after !== nothing && !cursor(page)
        page = datoms("eavt", Any[ms, ":measurement-set/measurements"]; db=db, offset=offset, limit=limit)
    end
    ids = Int.(page.v)
    stop = offset + length(ids)
    if cursor(page)
        nxt = metadata(page, "next")
        return ids, (nxt === nothing ? nothing : (stop, nxt))
    end
    ids, (length(ids) < limit ? nothing : (stop, nothing))
end

# Up to `limit` measurement ids of set `ms` in index order, and whether they
# reach the end of the set.
function set_measurement_ids(ms, db, limit)
    ids = Int[]
    pos = 0
    while length(ids) < limit
        got, pos = set_page(ms, db, min(PAGE, limit - length(ids)), pos)
        append!(ids, got)
        pos === nothing && return ids, true
    end
    more, _ = set_page(ms, db, 1, pos)
    ids, isempty(more)
end

# (count, exact): the set's measurement count when it is at most `cap`, else
# (cap, false). Counted in place on the server when it can; otherwise one page
# for a small set, and for a set filling the first page one probe past the cap,
# then one-datom probes binary-search the 10000-wide block it ends in and one
# page of that block finishes the count: about six round trips for any size,
# instead of a page per 10000.
function count_set_measurements(ms, db, cap)
    counted = datoms_count("eavt", Any[ms, ":measurement-set/measurements"]; db=db, cap=cap)
    if counted !== nothing
        n, exact = counted
        return exact && n <= cap ? (n, true) : (cap, false)
    end
    cap = min(cap, 100000)  # the deepest offset paging goes
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
    measurement_sets(; db=nothing, count_cap=1000000)

assay_name, assay_technology, measurement_set_name, measurement_count,
measurement_count_exact. Counting never binds a whole set in one query: a set
with more than `count_cap` measurements reports `count_cap` with
measurement_count_exact = false (read it as "at least"). The server counts in
place; one without in-place counts is paged instead, capped at 100000 (the
deepest it pages by offset).
"""
function measurement_sets(; db=nothing, count_cap=1000000, kwargs...)
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
        pos = 0
        for _ in 1:10  # 100000 measurements at most, the server's offset cap
            ids, pos = set_page(ms, db, attr === nothing ? n : PAGE, pos)
            complete = pos === nothing
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
samples (schema 0.3.18), so nothing walks the measurements. A database without
it uses unify-central's stored sample index when the service has a current
copy, else datoms seeks: slow on a large cohort (about two requests per sample
and set), but no query walks the measurements.
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
    found = sample_assays_from_index(db)
    found === nothing && (found = sample_assays_by_seek(db; kwargs...))
    nrow(found) == 0 && return DataFrame()
    smp = samples(; db=db)
    keep = intersect(["sample_id", "subject_id", "timepoint_id"], names(smp))
    out = rightjoin(select(smp, keep), found; on=:sample_id)
    keep_provenance(order_columns(out, ["subject_id", "sample_id", "timepoint_id"]), smp)
end

sample_assay_frame() = DataFrame(sample_id=String[], assay_name=String[], assay_technology=String[],
                                 measurement_set_name=String[])

# sample_assays rows from unify-central's stored sample index, or nothing when
# the service has no current copy.
function sample_assays_from_index(db)
    index = sample_index(db)
    index === nothing && return nothing
    sets = get(index, "measurement_sets", Dict()); smp = get(index, "samples", Dict())
    out = sample_assay_frame()
    for (ms, keys) in get(index, "set_samples", Dict()), k in keys
        (haskey(smp, k) && haskey(sets, ms)) || continue
        m = sets[ms]
        push!(out, (string(smp[k]["id"]), string(something(get(m, "assay", nothing), missing)),
                    string(something(get(m, "assay_type", nothing), missing)), string(get(m, "name", ms))))
    end
    out
end

# The links from a measurement set to its samples: the set's ref to its
# members, and each member's ref to its sample. Matrices list their samples
# directly.
const SET_LINKS = [(":measurement-set/measurements", ":measurement/sample"),
                   (":measurement-set/single-cells", ":single-cell/sample")]
const SEEK_STEPS = 8

function seek1(index, components, db)
    d = datoms(index, components; db=db, limit=1, seek=true)
    nrow(d) == 0 ? nothing : (e=Int(d.e[1]), a=d.a[1], v=d.v[1])
end

# The highest member id of set `ms` through `link`: one offset read for a set
# within the server's offset cap, else seeks, doubling then halving the
# distance from `lo` (each seek lands on the first member at or past an id).
function last_member(ms, link, lo, db)
    n = datoms_count("eavt", Any[ms, link]; db=db, cap=100001)
    if n !== nothing && n[2]
        return Int(datoms("eavt", Any[ms, link]; db=db, offset=n[1] - 1, limit=1).v[1])
    end
    function at(x)
        d = seek1("eavt", Any[ms, link, x], db)
        d !== nothing && d.a == link && d.e == ms ? Int(d.v) : nothing
    end
    low, step = lo, 1
    high = lo
    while true
        v = at(lo + step)
        if v === nothing
            high = lo + step
            break
        end
        low, step = v, step * 2
    end
    while high - low > 1
        mid = (low + high) ÷ 2
        v = at(mid)
        v === nothing ? (high = mid) : (low = v)
    end
    low
end

# One page of entity `e`'s `link` values from `pos` (an offset or a cursor), and
# the next position, or nothing at the end.
function link_page(e, link, db, pos)
    after = pos isa AbstractString ? pos : nothing
    page = datoms("eavt", Any[e, link]; db=db, offset=after === nothing ? pos : 0, limit=PAGE, after=after)
    vals = Int.(page.v)
    if "cursor" in metadatakeys(page) && metadata(page, "cursor")
        return vals, metadata(page, "next")
    end
    vals, (length(vals) < PAGE ? nothing : pos + length(vals))
end

# sample_assays rows for a database without :measurement-set/includes-samples
# or a stored sample index, from datoms seeks only: no query binds a set's or a
# sample's measurements. For each set its first and last member ids; then for
# each sample and set, a seek into the sample's members (VAET) at the set's
# first id, and one lookup of the hit's sets. Sets' ids interleave, so a hit in
# another set (which it also records) moves the seek past it, up to SEEK_STEPS
# times; past that, one query settles that sample alone. About two requests per
# sample and set, each also an audit write on the service: slow on a large
# cohort, but it never walks the data.
function sample_assays_by_seek(db; kwargs...)
    meta = runq(dq(["?ms", "?assay-name", "?assay-technology", "?measurement-set-name"],
                   [["?a", ":assay/name", "?assay-name"], ["?a", ":assay/technology", "?t"],
                    ["?t", ":db/ident", "?assay-technology"], ["?a", ":assay/measurement-sets", "?ms"],
                    ["?ms", ":measurement-set/name", "?measurement-set-name"]]); db=db, kwargs...)
    ids = runq(dq(["?s", "?sample-id"], [["?s", ":sample/id", "?sample-id"]]); db=db, kwargs...)
    (nrow(meta) == 0 || nrow(ids) == 0) && return sample_assay_frame()
    sample_id = Dict(Int(e) => String(x) for (e, x) in zip(ids.s, ids.sample_id))
    members = Dict(Int(ms) => Set{Int}() for ms in meta.ms)
    ranges = NamedTuple[]
    for ms in Int.(meta.ms)
        for (set_link, sample_link) in SET_LINKS
            first_page = datoms("eavt", Any[ms, set_link]; db=db, limit=1)
            nrow(first_page) == 0 && continue
            lo = Int(first_page.v[1])
            push!(ranges, (ms=ms, set_link=set_link, sample_link=sample_link, lo=lo, hi=last_member(ms, set_link, lo, db)))
        end
        for mm in Int.(datoms("eavt", Any[ms, ":measurement-set/measurement-matrices"]; db=db, limit=10000).v)
            pos = 0
            while pos !== nothing
                got, pos = link_page(mm, ":measurement-matrix/samples", db, pos)
                union!(members[ms], got)
            end
        end
    end
    owners(m, set_link) = Int.(datoms("vaet", Any[m, set_link]; db=db, limit=100).e)
    function check(s)
        hits = Set{Int}()
        for r in ranges
            r.ms in hits && continue
            pos, settled = r.lo, false
            for _ in 1:SEEK_STEPS
                d = seek1("vaet", Any[s, r.sample_link, pos], db)
                if d === nothing || d.a != r.sample_link || Int(d.v) != s || d.e > r.hi
                    settled = true
                    break
                end
                of = owners(d.e, r.set_link)
                union!(hits, of)
                if r.ms in of
                    settled = true
                    break
                end
                pos = d.e + 1
            end
            if !settled
                q = runq(dq(["?m"], [["?m", r.sample_link, "?s"], ["?ms", r.set_link, "?m"]];
                            in=["?s", "?ms"], args=[s, r.ms]); db=db, kwargs...)
                nrow(q) > 0 && push!(hits, r.ms)
            end
        end
        s, hits
    end
    for (s, hits) in asyncmap(check, collect(keys(sample_id)); ntasks=4), ms in hits
        haskey(members, ms) && push!(members[ms], s)
    end
    out = sample_assay_frame()
    for row in eachrow(meta), s in members[Int(row.ms)]
        haskey(sample_id, s) && push!(out, (sample_id[s], row.assay_name, ident_name(row.assay_technology),
                                            row.measurement_set_name))
    end
    out
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
    db = ensure_db(db)
    genes === nothing || return cnv_segments_run(db, genes, samples, subjects, nothing, Int[]; kwargs...)
    run_cnv_plan((ids, start, sets) -> cnv_segments_run(db, nothing, ids, subjects, start, sets; kwargs...),
                 ":measurement/cnv", samples, subjects, db; kwargs...)
end

function cnv_segments_run(db, genes, samples, subjects, start, sets; kwargs...)
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
    isempty(sets) || (push!(ins, ["?ms", "..."]); push!(args, sets))
    start === :samples && (where = lead_with(where, SAMPLE_LEAD))
    start === :sets && (where = lead_with(where, Any[["?ms", ":measurement-set/measurements", "?m"]]))
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
    db = ensure_db(db)
    genes === nothing || return cnv_gene_calls_run(db, genes, samples, subjects, measurement, nothing, Int[]; kwargs...)
    run_cnv_plan((ids, start, sets) -> cnv_gene_calls_run(db, nothing, ids, subjects, measurement, start, sets; kwargs...),
                 measurement_attr(measurement), samples, subjects, db; kwargs...)
end

function cnv_gene_calls_run(db, genes, samples, subjects, measurement, start, sets; kwargs...)
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
    isempty(sets) || (push!(ins, ["?ms", "..."]); push!(args, sets))
    start === :samples && (where = lead_with(where, SAMPLE_LEAD))
    start === :sets && (where = lead_with(where, Any[["?ms", ":measurement-set/measurements", "?m"]]))
    df = runq(dq(["?sample-id", "?measurement-set", "?hgnc-symbol", "?value"], where; in=ins, args=args,
                 with=["?m"]); db=db, kwargs...)
    eltype(df.value) <: Union{Missing,AbstractString} && (df.value = ident_name.(df.value))
    df
end
