# Reshaping measurement tables and joining context (samples, subjects,
# variants, ...), ported from wick.

"""
    LabeledMatrix(data, rows, cols)

A numeric matrix with row and column names (rows = samples, cols = targets by
default), as returned by `to_matrix`. `missing` where no value.
"""
struct LabeledMatrix
    data::Matrix{Union{Missing,Float64}}
    rows::Vector{String}
    cols::Vector{String}
end
Base.size(m::LabeledMatrix, args...) = size(m.data, args...)
Base.transpose(m::LabeledMatrix) = LabeledMatrix(permutedims(m.data), m.cols, m.rows)
Base.adjoint(m::LabeledMatrix) = transpose(m)
function Base.show(io::IO, m::LabeledMatrix)
    println(io, "LabeledMatrix $(length(m.rows)) x $(length(m.cols))")
    k = min(10, length(m.rows))
    show(io, DataFrame(hcat(m.rows[1:k], m.data[1:k, :]), vcat(["row"], m.cols)); truncate=12)
    length(m.rows) > k && print(io, "\n  ⋮ $(length(m.rows) - k) more rows")
end
Base.getindex(m::LabeledMatrix, r::AbstractString, c::AbstractString) =
    m.data[findfirst(==(r), m.rows), findfirst(==(c), m.cols)]
export LabeledMatrix

"""
    to_matrix(df, col; row="sample_id", value="value", fun=mean)

Long measurements to a LabeledMatrix (rows x cols); several `col` columns are
pasted with "|"; duplicate cells are aggregated with `fun`.
"""
function to_matrix(df::AbstractDataFrame, col; row="sample_id", value="value", fun=mean)
    cols = col isa AbstractVector ? String.(col) : [String(col)]
    ckey = length(cols) == 1 ? string.(df[!, cols[1]]) : [join(string.(collect(r)), "|") for r in eachrow(df[!, cols])]
    rkey = string.(df[!, row])
    vals = df[!, value]
    cells = Dict{Tuple{String,String},Vector{Float64}}()
    for (r, c, v) in zip(rkey, ckey, vals)
        ismissing(v) && continue
        push!(get!(cells, (r, c), Float64[]), Float64(v))
    end
    rn = sort(unique(rkey)); cn = sort(unique(ckey))
    ri = Dict(r => i for (i, r) in enumerate(rn)); ci = Dict(c => i for (i, c) in enumerate(cn))
    m = Matrix{Union{Missing,Float64}}(missing, length(rn), length(cn))
    for ((r, c), v) in cells
        m[ri[r], ci[c]] = fun(v)
    end
    LabeledMatrix(m, rn, cn)
end

"""
    to_long(m; row_name="sample_id", col_name="target", value_name="value")

LabeledMatrix to long format, without missing cells.
"""
function to_long(m::LabeledMatrix; row_name="sample_id", col_name="target", value_name="value")
    r = String[]; c = String[]; v = Float64[]
    for j in eachindex(m.cols), i in eachindex(m.rows)
        ismissing(m.data[i, j]) && continue
        push!(r, m.rows[i]); push!(c, m.cols[j]); push!(v, m.data[i, j])
    end
    DataFrame(row_name => r, col_name => c, value_name => v)
end

"""
    split_by_measurement_set(df; col=nothing, wide=true, fun=mean)

One element per measurement set (wick's group_by_assay_meas_set), optionally
cast to a LabeledMatrix.
"""
function split_by_measurement_set(df::AbstractDataFrame; col=nothing, wide::Bool=true, fun=mean)
    col === nothing && (col = setdiff(names(df), ["sample_id", "measurement_set", "value"]))
    out = Dict{String,Any}()
    for g in groupby(df, :measurement_set)
        out[g.measurement_set[1]] = wide ? to_matrix(DataFrame(g), col; fun=fun) : DataFrame(g)
    end
    out
end

"""
    select_targets(m; include=nothing, exclude=nothing)
"""
function select_targets(m::LabeledMatrix; include=nothing, exclude=nothing)
    keep = trues(length(m.cols))
    include === nothing || (keep .&= in(Set(include)).(m.cols))
    exclude === nothing || (keep .&= .!in(Set(exclude)).(m.cols))
    LabeledMatrix(m.data[:, keep], m.rows, m.cols[keep])
end

ensure_long(tab, col_name="target") = tab isa LabeledMatrix ? to_long(tab; col_name=col_name) : tab

"""
    add_sample_context(tab; db=nothing, include_subjects=true, include_outcomes=false)
    add_subject_context(tab; db=nothing, include_outcomes=false)
    add_variant_context(tab; db=nothing)
    add_cnv_context(tab; db=nothing)

Join sample attributes (optionally subject attributes and outcomes) by
sample_id; subject attributes by subject_id; variant annotations by variant_id;
CNV segments by cnv_id. Provenance of `tab` is kept.
"""
function add_sample_context(tab; db=nothing, include_subjects::Bool=true, include_outcomes::Bool=false)
    tab = ensure_long(tab)
    smp = samples(; db=db)
    smp = select(smp, vcat(["sample_id"], setdiff(names(smp), names(tab))))
    out = leftjoin(tab, smp; on=:sample_id, matchmissing=:notequal)
    if include_subjects && "subject_id" in names(out)
        out = add_subject_context(out; db=db, include_outcomes=include_outcomes)
    end
    keep_provenance(out, tab)
end

function add_subject_context(tab; db=nothing, include_outcomes::Bool=false)
    "subject_id" in names(tab) || error("tab has no subject_id column")
    sub = subjects(; db=db)
    sub = select(sub, vcat(["subject_id"], setdiff(names(sub), names(tab))))
    out = leftjoin(tab, sub; on=:subject_id, matchmissing=:notequal)
    if include_outcomes
        oc = subject_outcomes(; db=db)
        oc = select(oc, vcat(["subject_id"], setdiff(names(oc), names(out))))
        out = leftjoin(out, oc; on=:subject_id, matchmissing=:notequal)
    end
    keep_provenance(out, tab)
end

function add_variant_context(tab; db=nothing)
    tab = ensure_long(tab, "variant_id")
    "variant_id" in names(tab) || error("tab has no variant_id column")
    va = variant_annotations(; db=db, variant_ids=unique(tab.variant_id))
    va = select(va, vcat(["variant_id"], setdiff(names(va), names(tab))))
    keep_provenance(leftjoin(tab, va; on=:variant_id), tab)
end

function add_cnv_context(tab; db=nothing)
    tab = ensure_long(tab, "cnv_id")
    "cnv_id" in names(tab) || error("tab has no cnv_id column")
    cn = cnvs(; db=db)
    cn = select(cn, vcat(["cnv_id"], setdiff(names(cn), names(tab))))
    keep_provenance(leftjoin(tab, cn; on=:cnv_id), tab)
end

function taxonomy_levels(taxa)
    lv = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    prefix = any(startswith("otu_"), names(taxa)) ? "otu_" : "sgb_"
    [prefix * l for l in lv if prefix * l in names(taxa)]
end

"""
    deduplicate_taxonomy(taxa; na_value="Unclassified", sep="_")

Make taxon names unique at each level by prefixing ancestors where the same
name occurs under different parents (OTU `otu_*` or SGB `sgb_*` columns).
"""
function deduplicate_taxonomy(taxa::AbstractDataFrame; na_value="Unclassified", sep="_")
    levels = taxonomy_levels(taxa)
    taxa = copy(taxa); out = copy(taxa)
    for (i, cur) in enumerate(levels)
        taxa[!, cur] = coalesce.(taxa[!, cur], na_value)
        out[!, cur] = coalesce.(out[!, cur], na_value)
        uniq = unique(taxa[!, levels[1:i]])
        counts = Dict{Any,Int}()
        for v in uniq[!, cur]; counts[v] = get(counts, v, 0) + 1; end
        dup = Set(k for (k, n) in counts if n > 1)
        for r in 1:nrow(out)
            out[r, cur] in dup && (out[r, cur] = join([string(taxa[r, l]) for l in levels[1:i]], sep))
        end
    end
    out
end

"""
    aggregate_taxa(tab, taxa, id_col; normalize=true, wide=false, na_value="Unclassified")

Sum measurements to each taxonomic level (optionally per-sample proportions):
a Dict level => long DataFrame (or LabeledMatrix with `wide=true`).
"""
function aggregate_taxa(tab::AbstractDataFrame, taxa::AbstractDataFrame, id_col; normalize::Bool=true,
                        wide::Bool=false, na_value="Unclassified")
    taxa = deduplicate_taxonomy(taxa; na_value=na_value)
    tab = innerjoin(dropmissing(tab, :value), taxa; on=Symbol(id_col))
    out = Dict{String,Any}()
    for cur in taxonomy_levels(taxa)
        m = combine(groupby(tab, [:sample_id, Symbol(cur)]), :value => sum => :value)
        rename!(m, cur => "taxon")
        if normalize
            m = transform(groupby(m, :sample_id), :value => (v -> v ./ sum(v)) => :value)
        end
        if wide
            lm = to_matrix(m, "taxon"; fun=sum)
            lm.data[ismissing.(lm.data)] .= 0.0
            out[cur] = lm
        else
            out[cur] = m
        end
    end
    out
end
