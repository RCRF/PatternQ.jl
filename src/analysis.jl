# Broad-scale descriptive expression analysis, backported from the Clojure
# variant-forensics reports (unify-central/analysis): a sample vs a reference
# cohort (z-scores, percentile ranks, top genes), geneset views against cohort
# distributions, two-sample change (log fold change, MA plots), ssGSEA, top
# varying genes, expression distances / nearest samples.

const GENESET_DIR = joinpath(@__DIR__, "..", "data", "genesets")

"""
    genesets()
    geneset(name)

Gene sets shipped with patternq (hallmark apoptosis / DNA repair / hypoxia /
inflammatory, antibody therapy targets, housekeeping, germline multi-cancer,
melanoma phenotypes, neural crest, adult kidney): a Dict name => symbols.
"""
function genesets()
    out = Dict{String,Vector{String}}()
    for f in readdir(GENESET_DIR)
        endswith(f, ".txt") || continue
        g = filter(!isempty, strip.(readlines(joinpath(GENESET_DIR, f))))
        out[replace(f, r"\.txt$" => "")] = String.(g)
    end
    out
end

function geneset(name::AbstractString)
    gs = genesets()
    haskey(gs, name) || error("Unknown gene set '$name'")
    gs[name]
end

"""
    sample_expression(sample; db=nothing, measurement="tpm", measurement_set=nothing, genes=nothing)

Expression of one sample as a Dict hgnc_symbol => value (values of several
gene products of the same gene summed).
"""
function sample_expression(sample::AbstractString; db=nothing, measurement="tpm", measurement_set=nothing,
                           genes=nothing)
    gx = gene_expression(; db=db, genes=genes, samples=[sample], measurement=measurement,
                         measurement_set=measurement_set)
    out = Dict{String,Float64}()
    for r in eachrow(gx)
        out[r.hgnc_symbol] = get(out, r.hgnc_symbol, 0.0) + r.value
    end
    out
end

"""
    percentile_rank(values, x)

Mid-rank percentile of x within values: 100 * (count below + half the ties) / n.
"""
function percentile_rank(values, x::Real)
    v = collect(skipmissing(values))
    isempty(v) && return NaN
    100 * (count(<(x), v) + 0.5 * count(==(x), v)) / length(v)
end

"""
    compare_to_cohort(sample; db=nothing, cohort_db, measurement="tpm", cohort_measurement=measurement,
                      genes=nothing, log=true, fill_missing=true, anchor_gene="GAPDH",
                      sd_floor=0.25, batch_size=2000, min_cohort=5)

Place a sample's expression within a reference cohort's distribution, gene by
gene: z-score and percentile on log2(1 + x). The cohort is streamed in gene
batches. With `fill_missing`, genes without a stored value in a cohort sample
count as 0 (imports often omit zeros); cohort samples are those with a value
for `anchor_gene`. `sd_floor` keeps near-constant genes from producing
enormous z-scores. Columns: hgnc_symbol, value, z, percentile, cohort_n,
cohort_observed (cohort samples with a stored value; low coverage of a
usually-expressed gene points to an import problem), cohort_mean, cohort_sd,
cohort_median.
"""
function compare_to_cohort(sample::AbstractString; db=nothing, cohort_db, measurement="tpm",
                           cohort_measurement=measurement, genes=nothing, measurement_set=nothing,
                           cohort_measurement_set=nothing, log::Bool=true, fill_missing::Bool=true,
                           anchor_gene="GAPDH", sd_floor=0.25, batch_size=2000, min_cohort=5)
    x = sample_expression(sample; db=db, measurement=measurement, measurement_set=measurement_set, genes=genes)
    isempty(x) && error("No $measurement expression for sample $sample")
    tr(v) = log ? log2(1 + max(v, 0)) : v
    ncohort = 0
    if fill_missing
        anchor = gene_expression(; db=cohort_db, genes=[anchor_gene], measurement=cohort_measurement,
                                 measurement_set=cohort_measurement_set)
        ncohort = length(unique(anchor.sample_id))
        ncohort == 0 && error("No $cohort_measurement values for $anchor_gene in $cohort_db; pass anchor_gene or fill_missing=false")
    end
    syms = collect(keys(x))
    rows = NamedTuple[]
    for b in Iterators.partition(syms, batch_size)
        cg = gene_expression(; db=cohort_db, genes=collect(b), measurement=cohort_measurement,
                             measurement_set=cohort_measurement_set)
        nrow(cg) == 0 && continue
        agg = combine(groupby(cg, [:sample_id, :hgnc_symbol]), :value => sum => :value)
        for g in groupby(agg, :hgnc_symbol)
            vals = tr.(g.value)
            observed = length(vals)
            fill_missing && (vals = vcat(vals, fill(tr(0.0), max(0, ncohort - observed))))
            gene = g.hgnc_symbol[1]
            push!(rows, (hgnc_symbol=gene, value=x[gene], cohort_n=length(vals), cohort_observed=observed,
                         cohort_mean=mean(vals), cohort_sd=length(vals) > 1 ? std(vals) : NaN,
                         cohort_median=median(vals), percentile=percentile_rank(vals, tr(x[gene]))))
        end
    end
    out = DataFrame(rows)
    out = filter(r -> r.cohort_n >= min_cohort, out)
    sdv = [max(isnan(s) ? 0.0 : s, sd_floor) for s in out.cohort_sd]
    out.z = (tr.(out.value) .- out.cohort_mean) ./ sdv
    out = select(out, :hgnc_symbol, :value, :z, :percentile, :cohort_n, :cohort_observed, :cohort_mean,
                 :cohort_sd, :cohort_median)
    sort!(out, order(:z, by=abs, rev=true))
    metadata!(out, "patternq_comparison", (sample=sample, db=ensure_db(db), cohort_db=cohort_db,
                                           measurement=measurement, cohort_measurement=cohort_measurement,
                                           log=log, cohort_size=ncohort); style=:note)
    out
end

"""
    top_by_zscore(comparison; n=25, direction=:both, min_value=1, min_observed=0.5)

Top genes by z-score (`:both` = largest |z|, `:up`, `:down`), ignoring genes
whose sample value and cohort median are both below `min_value`, and genes
observed in less than `min_observed` of the cohort.
"""
function top_by_zscore(cmp::AbstractDataFrame; n=25, direction::Symbol=:both, min_value=1, min_observed=0.5)
    c = filter(r -> !ismissing(r.z) && !isnan(r.z), cmp)
    c = filter(r -> r.cohort_observed >= min_observed * r.cohort_n, c)
    c = filter(r -> r.value >= min_value || (2^r.cohort_median - 1) >= min_value, c)
    c = direction == :both ? sort(c, order(:z, by=abs, rev=true)) :
        direction == :up ? filter(r -> r.z > 0, sort(c, :z, rev=true)) :
        filter(r -> r.z < 0, sort(c, :z))
    keep_provenance(first(c, min(n, nrow(c))), cmp)
end

"""
    log_fold_change(a, b; pseudocount=1)

a, b: Dicts gene => value (missing genes count as 0). Columns: hgnc_symbol,
value_a, value_b, lfc = log2((b + pc) / (a + pc)), avg_log10 = mean of
log10(1 + a) and log10(1 + b) (the MA plot's x, as in the Clojure reports).
"""
function log_fold_change(a::AbstractDict, b::AbstractDict; pseudocount=1)
    genes = unique(vcat(collect(keys(a)), collect(keys(b))))
    va = [get(a, g, 0.0) for g in genes]; vb = [get(b, g, 0.0) for g in genes]
    DataFrame(hgnc_symbol=genes, value_a=va, value_b=vb,
              lfc=log2.((vb .+ pseudocount) ./ (va .+ pseudocount)),
              avg_log10=(log10.(1 .+ va) .+ log10.(1 .+ vb)) ./ 2)
end

"""
    compare_samples(sample_a, sample_b; db=nothing, db_b=db, measurement="tpm", min_avg=0.5, pseudocount=1)

Expression change from sample_a (reference) to sample_b: `log_fold_change`
output with avg_log10 >= min_avg, ordered by lfc (largest first).
"""
function compare_samples(sample_a::AbstractString, sample_b::AbstractString; db=nothing, db_b=db,
                         measurement="tpm", min_avg=0.5, pseudocount=1)
    a = sample_expression(sample_a; db=db, measurement=measurement)
    b = sample_expression(sample_b; db=db_b, measurement=measurement)
    out = filter(r -> r.avg_log10 >= min_avg, log_fold_change(a, b; pseudocount=pseudocount))
    sort!(out, :lfc, rev=true)
    metadata!(out, "patternq_comparison", (sample_a=sample_a, sample_b=sample_b, measurement=measurement); style=:note)
    out
end

"""
    ssgsea_score(x, gene_set; alpha=0.25)

Single-sample GSEA (Barbie et al. 2009): genes ranked by expression; sums the
difference between the weighted (rank^alpha, highest expression = highest
rank) running distribution of the set and that of the other genes. Same
definition as the R, Python and Clojure libraries. `NaN` when no gene (or every
gene) is in the set.
"""
function ssgsea_score(x::AbstractDict, gene_set; alpha=0.25)
    g = [k for (k, v) in x if !ismissing(v)]
    v = [Float64(x[k]) for k in g]
    o = sortperm(v, rev=true, alg=MergeSort)   # stable: ties keep first-seen order
    n = length(g)
    gs = Set(gene_set)
    hits = [g[i] in gs for i in o]
    nh = count(hits)
    (nh == 0 || nh == n) && return NaN
    asc = sortperm(v, alg=MergeSort)        # R's rank(ties.method = "first")
    rank = zeros(Int, n); rank[asc] = 1:n
    r = rank[o]                             # highest expression -> rank n
    w = [hits[i] ? r[i]^alpha : 0.0 for i in 1:n]
    p_hit = cumsum(w) ./ sum(w)
    p_miss = cumsum(.!hits) ./ (n - nh)
    sum(p_hit .- p_miss)
end

"""
    ssgsea(m, gene_sets=genesets(); alpha=0.25)

ssGSEA scores for every sample (row of a LabeledMatrix) and gene set: a
LabeledMatrix samples x gene sets.
"""
function ssgsea(m::LabeledMatrix, gene_sets::AbstractDict=genesets(); alpha=0.25)
    names_ = sort(collect(keys(gene_sets)))
    out = Matrix{Union{Missing,Float64}}(missing, length(m.rows), length(names_))
    for i in eachindex(m.rows)
        row = Dict(m.cols[j] => m.data[i, j] for j in eachindex(m.cols) if !ismissing(m.data[i, j]))
        for (j, s) in enumerate(names_)
            out[i, j] = ssgsea_score(row, gene_sets[s]; alpha=alpha)
        end
    end
    LabeledMatrix(out, m.rows, names_)
end

"""
    top_varying_genes(m; n=500)

Variance of log2(1 + x) across samples (rows): hgnc_symbol, mean, variance.
"""
function top_varying_genes(m::LabeledMatrix; n=500)
    stats = [(let v = [log2(1 + max(x, 0)) for x in skipmissing(m.data[:, j])]
                  (hgnc_symbol=m.cols[j], mean=isempty(v) ? NaN : mean(v),
                   variance=length(v) > 1 ? var(v) : NaN)
              end) for j in eachindex(m.cols)]
    df = sort(DataFrame(stats), :variance, rev=true, lt=(a, b) -> isnan(b) ? !isnan(a) : a < b)
    first(df, min(n, nrow(df)))
end

"""
    expression_distance(a, b; method=:cosine)

Distance between two profiles (Dicts gene => value); genes missing from b
count as 0. `:cosine` (1 - cosine similarity) or `:euclidean`.
"""
function expression_distance(a::AbstractDict, b::AbstractDict; method::Symbol=:cosine)
    a = Dict(k => v for (k, v) in a if !ismissing(v)); b = Dict(k => v for (k, v) in b if !ismissing(v))
    ks = collect(keys(a))
    va = [a[k] for k in ks]; vb = [get(b, k, 0.0) for k in ks]
    method == :euclidean && return sqrt(sum((va .- vb) .^ 2))
    ma = sqrt(sum(va .^ 2)); mb = sqrt(sum(abs2, values(b)))
    (ma == 0 || mb == 0) && return 1.0
    1 - sum(va .* vb) / (ma * mb)
end

"""
    nearest_samples(x, m; method=:cosine, n=10)

Samples (rows of a LabeledMatrix) closest to the profile x. Raw-expression
nearest neighbours are sensitive to batch, vendor and pipeline effects:
compare within a consistently processed cohort.
"""
function nearest_samples(x::AbstractDict, m::LabeledMatrix; method::Symbol=:cosine, n=10)
    d = [expression_distance(x, Dict(m.cols[j] => coalesce(m.data[i, j], 0.0) for j in eachindex(m.cols));
                             method=method) for i in eachindex(m.rows)]
    o = sortperm(d)
    first(DataFrame(sample_id=m.rows[o], distance=d[o]), min(n, length(o)))
end

"""
    examine_geneset(genes, samples; db=nothing, cohort_dbs, measurement="tpm",
                    cohort_measurement=measurement, type=:violin, title=nothing)

Samples against reference cohort distributions for a gene set (or a
`genesets()` name): fetches expression and draws `plot_vs_cohort`.
`cohort_dbs`: a database name, a vector of names, or a Dict label => db.
"""
function examine_geneset(genes, samples; db=nothing, cohort_dbs, measurement="tpm",
                         cohort_measurement=measurement, type::Symbol=:violin, title=nothing)
    if genes isa AbstractString && haskey(genesets(), genes)
        title === nothing && (title = genes)
        genes = geneset(genes)
    end
    genes = genes isa AbstractString ? [genes] : collect(genes)
    samples = samples isa AbstractString ? [samples] : collect(samples)
    se = gene_expression(; db=db, genes=genes, samples=samples, measurement=measurement)
    cohorts = cohort_dbs isa AbstractDict ? collect(cohort_dbs) :
              [c => c for c in (cohort_dbs isa AbstractString ? [cohort_dbs] : cohort_dbs)]
    parts = DataFrame[]
    for (label, cdb) in cohorts
        d = gene_expression(; db=cdb, genes=genes, measurement=cohort_measurement)
        nrow(d) > 0 && (d.cohort .= label; push!(parts, d))
    end
    ce = isempty(parts) ? DataFrame() : vcat(parts...)
    plot_vs_cohort(se, ce; type=type, title=something(title, "Samples vs cohort expression"),
                   xlab="$measurement (log scale)")
end
