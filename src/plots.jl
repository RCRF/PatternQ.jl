# plotly plots, hand-rolled on PlotlyBase (mirrors the R and Python plots).
# Each returns a PlotlyBase.Plot. The theme (theme/plotly-theme.json, shared by
# all patternq libraries) uses a colorblind-validated categorical order: colors
# follow entities in fixed order and are never cycled; series past 8 fold into
# "Other".

const THEME_PATH = joinpath(@__DIR__, "..", "data", "plotly-theme.json")
const _theme = Ref{Any}(nothing)
const OTHER_COLOR = "#8a8983"

"""
    plot_theme()

Colors and fonts used by patternq plots (shared with R, Python and Clojure).
"""
function plot_theme()
    _theme[] === nothing && (_theme[] = tonative(JSON3.read(read(THEME_PATH, String))))
    _theme[]
end

series_colors(n) = (c = String.(plot_theme()["categorical"]); vcat(c, fill(OTHER_COLOR, max(0, n - length(c))))[1:n])

function fold_other(x; maxn=8)
    vals = collect(skipmissing(x))
    counts = Dict{Any,Int}()
    for v in vals; counts[v] = get(counts, v, 0) + 1; end
    length(counts) <= maxn && return x
    keep = Set(first.(sort(collect(counts), by=p -> -p[2]))[1:maxn-1])
    [ismissing(v) ? missing : (v in keep ? string(v) : "Other") for v in x]
end

sequential_scale() = (s = plot_theme()["sequential"]; [[(i - 1) / (length(s) - 1), s[i]] for i in eachindex(s)])
diverging_scale() = (d = plot_theme()["diverging"]; [[0, d["low"]], [0.5, d["mid"]], [1, d["high"]]])

function alpha_color(hex, a)
    h = lstrip(hex, '#')
    "rgba($(parse(Int, h[1:2], base=16)),$(parse(Int, h[3:4], base=16)),$(parse(Int, h[5:6], base=16)),$(round(a, digits=2)))"
end

function theme_axis(extra)
    th = plot_theme()
    base = Dict{Symbol,Any}(:gridcolor => th["grid"], :zerolinecolor => th["grid"], :linecolor => th["grid"],
                            :tickfont => attr(color=th["text_secondary"]))
    merge(base, Dict{Symbol,Any}(pairs(extra)))
end

function pq_layout(title; xaxis=(;), yaxis=(;), kwargs...)
    th = plot_theme()
    Layout(; title=attr(text=something(title, ""), x=0, xanchor="left"),
           font=attr(family=th["font"], color=th["text_primary"]),
           colorway=th["categorical"], paper_bgcolor=th["surface"], plot_bgcolor=th["surface"],
           hoverlabel=attr(font=attr(family=th["font"])),
           xaxis=attr(; theme_axis(xaxis)...), yaxis=attr(; theme_axis(yaxis)...), kwargs...)
end

axis_title(t) = attr(text=t)

"""
    plot_vaf_histogram(variants; samples=nothing, bin_size=0.05, title="VAF histogram")

Overlaid VAF histograms per sample, from `variants` output.
"""
function plot_vaf_histogram(v::AbstractDataFrame; samples=nothing, bin_size=0.05, title="VAF histogram")
    samples === nothing || (v = filter(r -> r.sample_id in samples, v))
    ids = sort(unique(v.sample_id))
    tr = [histogram(x=v.vaf[v.sample_id .== id], name=id, opacity=0.7,
                    xbins=attr(start=0, size=bin_size, var"end"=1),
                    marker=attr(color=c, line=attr(color=plot_theme()["surface"], width=1)))
          for (id, c) in zip(ids, series_colors(length(ids)))]
    Plot(tr, pq_layout(title; barmode="overlay", showlegend=length(ids) > 1,
                       xaxis=(title=axis_title("VAF"), range=[0, 1]), yaxis=(title=axis_title("variants"),)))
end

"""
    plot_gene_expression(expr; title="Gene expression", ylab="value", log=false)

Grouped bars of expression per gene and sample, from `gene_expression` output.
"""
function plot_gene_expression(expr::AbstractDataFrame; title="Gene expression", ylab="value", log::Bool=false)
    agg = combine(groupby(expr, [:sample_id, :hgnc_symbol]), :value => sum => :value)
    ids = sort(unique(agg.sample_id))
    tr = [(a = agg[agg.sample_id .== id, :]; bar(x=a.hgnc_symbol, y=a.value, name=id, marker=attr(color=c)))
          for (id, c) in zip(ids, series_colors(length(ids)))]
    Plot(tr, pq_layout(title; barmode="group", bargap=0.2, bargroupgap=0.05, showlegend=length(ids) > 1,
                       yaxis=(title=axis_title(ylab), type=log ? "log" : "linear"), xaxis=(title=axis_title(""),)))
end

"""
    plot_sample_overview(sample_assays; title=...)

Heatmap of sample counts per subject (rows) and measurement set (columns);
subjects without samples in a set are blank.
"""
function plot_sample_overview(sa::AbstractDataFrame; title="Samples per subject and measurement set")
    u = unique(select(sa, :subject_id, :sample_id, :measurement_set_name))
    u.n .= 1.0
    m = to_matrix(u, "measurement_set_name"; row="subject_id", value="n", fun=sum)
    z = [ismissing(x) ? nothing : x for x in m.data]
    Plot(heatmap(x=m.cols, y=m.rows, z=z, zmin=0, colorscale=sequential_scale(), xgap=1, ygap=1,
                 colorbar=attr(title=attr(text="samples")),
                 hovertemplate="%{y}<br>%{x}<br>%{z} samples<extra></extra>"),
         pq_layout(title; xaxis=(title=axis_title(""), tickangle=-30, automargin=true),
                   yaxis=(title=axis_title("subject"), autorange="reversed", automargin=true)))
end

function timepoint_order(tab, timepoints)
    timepoints !== nothing && return [t for t in timepoints if t in tab.timepoint_id]
    if "timepoint_relative_order" in names(tab)
        s = sort(unique(select(tab, :timepoint_id, :timepoint_relative_order)), :timepoint_relative_order)
        return unique(s.timepoint_id)
    end
    sort(unique(tab.timepoint_id))
end

"""
    plot_by_timepoint(tab; group=nothing, lines=false, value="value", timepoints=nothing,
                      levels=nothing, title=nothing, ylab="value")

Box plots of values per timepoint (ordered by relative order or `timepoints`),
optionally split by a group column (e.g. "bor", "status_1y"); `lines=true`
draws per-subject trajectories, colored by group when given.
"""
function plot_by_timepoint(tab::AbstractDataFrame; group=nothing, lines::Bool=false, value="value",
                           timepoints=nothing, levels=nothing, title=nothing, ylab="value")
    tab = filter(r -> !ismissing(r[value]) && !ismissing(r.timepoint_id), tab)
    timepoints === nothing || (tab = filter(r -> r.timepoint_id in timepoints, tab))
    ord = timepoint_order(tab, timepoints)
    rank = Dict(t => i for (i, t) in enumerate(ord))
    traces = GenericTrace[]
    trajectory(s, color, legendgroup, label) = begin
        s = sort(s, order(:timepoint_id, by=t -> rank[t]))
        nrow(s) < 2 && return
        push!(traces, scatter(x=s.timepoint_id, y=s[!, value], mode="lines", legendgroup=legendgroup,
                              line=attr(color=color, width=1), showlegend=false, hoverinfo="text",
                              text=fill(label, nrow(s))))
    end
    if group === nothing
        if lines && "subject_id" in names(tab)
            for sid in unique(tab.subject_id)
                trajectory(tab[tab.subject_id .== sid, :], "rgba(82,81,78,0.25)", "", sid)
            end
        end
        c = series_colors(1)[1]
        push!(traces, box(x=tab.timepoint_id, y=tab[!, value], name=ylab, marker=attr(color=c),
                          line=attr(color=c), boxpoints="all", jitter=0.3, pointpos=0, showlegend=false))
        lay = pq_layout(title)
    else
        g = fold_other(tab[!, group])
        lv = levels === nothing ? sort(unique(skipmissing(g))) : [l for l in levels if l in skipmissing(g)]
        cols = series_colors(length(lv))
        if lines && "subject_id" in names(tab)
            for (l, c) in zip(lv, cols)
                sub = tab[coalesce.(g .== l, false), :]
                for sid in unique(sub.subject_id)
                    trajectory(sub[sub.subject_id .== sid, :], alpha_color(c, 0.35), l, "$sid<br>$l")
                end
            end
        end
        for (l, c) in zip(lv, cols)
            s = tab[coalesce.(g .== l, false), :]
            push!(traces, box(x=s.timepoint_id, y=s[!, value], name=l, legendgroup=l,
                              marker=attr(color=c, size=4), line=attr(color=c), fillcolor=alpha_color(c, 0.15),
                              boxpoints=lines ? false : "all", jitter=0.3, pointpos=0))
        end
        lay = pq_layout(title; boxmode=lines ? "overlay" : "group")
    end
    lay[:xaxis] = attr(; theme_axis((title=axis_title("timepoint"), type="category",
                                     categoryorder="array", categoryarray=ord))...)
    lay[:yaxis] = attr(; theme_axis((title=axis_title(ylab),))...)
    Plot(traces, lay)
end

"""
    plot_by_group(tab, group; violin=false, title=nothing, ylab="value")

Box (or violin) plot of values per group, e.g. a measurement by best overall
response.
"""
function plot_by_group(tab::AbstractDataFrame, group; violin::Bool=false, title=nothing, ylab="value")
    tab = filter(r -> !ismissing(r.value) && !ismissing(r[group]), tab)
    g = fold_other(string.(tab[!, group]))
    lv = sort(unique(g))
    traces = GenericTrace[]
    for (l, c) in zip(lv, series_colors(length(lv)))
        y = tab.value[g .== l]
        push!(traces, violin ?
            PlotlyBase.violin(x=fill(l, length(y)), y=y, name=l, line=attr(color=c, width=1.5),
                              fillcolor=alpha_color(c, 0.25),
                              box=attr(visible=true, fillcolor=plot_theme()["surface"], line=attr(color=c), width=0.15),
                              meanline=attr(visible=false), points="all", jitter=0.4, pointpos=0,
                              marker=attr(color=c, size=5, opacity=0.8)) :
            box(x=fill(l, length(y)), y=y, name=l, marker=attr(color=c), line=attr(color=c),
                boxpoints="all", jitter=0.3, pointpos=0))
    end
    Plot(traces, pq_layout(title; showlegend=false, xaxis=(title=axis_title(string(group)),),
                           yaxis=(title=axis_title(ylab),)))
end

format_p(p) = isnan(p) ? "NA" : p < 1e-4 ? "< 1e-04" : string(round(p, sigdigits=2))

"""
    plot_survival(tab; time="os", event="os_event", group=nothing, title="Survival", xlab=time,
                  pvalue=true, levels=nothing)

Kaplan-Meier curves (hand-rolled), optionally by group, with censoring ticks
and the log-rank p-value when there are 2+ groups.
"""
function plot_survival(tab::AbstractDataFrame; time="os", event="os_event", group=nothing, title="Survival",
                       xlab=time, pvalue::Bool=true, levels=nothing)
    tab = filter(r -> !ismissing(r[time]) && !ismissing(r[event]), tab)
    g = group === nothing ? fill("all", nrow(tab)) : fold_other(tab[!, group])
    keep = .!ismissing.(g)
    tab = tab[keep, :]; g = String.(g[keep])
    lv = levels === nothing ? sort(unique(g)) : [l for l in levels if l in g]
    traces = GenericTrace[]
    for (l, c) in zip(lv, series_colors(length(lv)))
        sel = g .== l
        km = kaplan_meier(tab[sel, time], tab[sel, event])
        nm = "$l (n=$(count(sel)))"
        push!(traces, scatter(x=vcat(0, km.time), y=vcat(1, km.surv), mode="lines", name=nm,
                              line=attr(shape="hv", color=c, width=2),
                              hovertemplate="$nm<br>t=%{x:.1f}<br>S=%{y:.2f}<extra></extra>"))
        cens = km[km.n_censor .> 0, :]
        nrow(cens) > 0 && push!(traces, scatter(x=cens.time, y=cens.surv, mode="markers", showlegend=false,
                                                marker=attr(symbol="line-ns-open", size=9, color=c), hoverinfo="skip"))
    end
    ann = Any[]
    if pvalue && length(lv) > 1
        lr = logrank_test(tab[!, time], tab[!, event], g)
        push!(ann, attr(xref="paper", yref="paper", x=0.02, y=0.04, xanchor="left", showarrow=false,
                        text="log-rank p = $(format_p(lr.p))", font=attr(color=plot_theme()["text_secondary"])))
    end
    Plot(traces, pq_layout(title; showlegend=length(lv) > 1, annotations=ann,
                           xaxis=(title=axis_title(xlab), rangemode="tozero"),
                           yaxis=(title=axis_title("survival probability"), range=[0, 1.02])))
end

# -- clustering for heatmaps (average linkage, as R's hclust(method = "average"))

function hclust_order(x::AbstractMatrix{<:Real})
    n = size(x, 1)
    n < 3 && return collect(1:n)
    d = [sqrt(sum((x[i, :] .- x[j, :]) .^ 2)) for i in 1:n, j in 1:n]
    clusters = Dict(i => [i] for i in 1:n)
    sizes = Dict(i => 1 for i in 1:n)
    dist = Dict{Tuple{Int,Int},Float64}()
    for i in 1:n, j in i+1:n; dist[(i, j)] = d[i, j]; end
    next = n + 1
    while length(clusters) > 1
        _, (a, b) = findmin(dist)
        # merged cluster keeps a's members first
        clusters[next] = vcat(clusters[a], clusters[b])
        sizes[next] = sizes[a] + sizes[b]
        for k in keys(clusters)
            (k == a || k == b || k == next) && continue
            dak = dist[minmax(a, k)]; dbk = dist[minmax(b, k)]
            dist[minmax(k, next)] = (sizes[a] * dak + sizes[b] * dbk) / (sizes[a] + sizes[b])
        end
        for k in collect(keys(dist))
            (a in k || b in k) && delete!(dist, k)
        end
        delete!(clusters, a); delete!(clusters, b)
        next += 1
    end
    only(values(clusters))
end

function scale_rows(m::Matrix{Union{Missing,Float64}})
    out = copy(m)
    for i in 1:size(m, 1)
        v = collect(skipmissing(m[i, :]))
        mu = isempty(v) ? 0.0 : mean(v)
        s = length(v) > 1 ? std(v) : 1.0
        (isnan(s) || s == 0) && (s = 1.0)
        out[i, :] = (m[i, :] .- mu) ./ s
    end
    out
end

"""
    plot_heatmap(m::LabeledMatrix; scale=:none, cluster_rows=true, cluster_cols=true, title=nothing,
                 zlab=nothing, col_groups=nothing)

Clustered heatmap (average linkage). With `scale=:row`/`:column` values are
z-scored and drawn on the diverging scale; otherwise sequential.
`col_groups` (Dict column name => group, e.g. survival status) adds an
annotation strip above the heatmap.
"""
function plot_heatmap(m::LabeledMatrix; scale::Symbol=:none, cluster_rows::Bool=true, cluster_cols::Bool=true,
                      title=nothing, zlab=nothing, col_groups=nothing)
    (isempty(m.rows) || isempty(m.cols)) && error("plot_heatmap: empty matrix")
    data = copy(m.data)
    scale == :row && (data = scale_rows(data))
    scale == :column && (data = permutedims(scale_rows(permutedims(data))))
    filled = [ismissing(x) || !isfinite(x) ? 0.0 : Float64(x) for x in data]
    ro = cluster_rows ? hclust_order(filled) : collect(eachindex(m.rows))
    co = cluster_cols ? hclust_order(permutedims(filled)) : collect(eachindex(m.cols))
    data = data[ro, co]; rows = m.rows[ro]; cols = m.cols[co]
    div = scale != :none
    lim = div ? maximum(abs, skipmissing(data)) : nothing
    z = [ismissing(x) ? nothing : x for x in data]
    hm = heatmap(x=cols, y=rows, z=z, colorscale=div ? diverging_scale() : sequential_scale(),
                 zmin=div ? -lim : nothing, zmax=div ? lim : nothing,
                 colorbar=attr(title=attr(text=something(zlab, div ? "z-score" : "value"))),
                 hovertemplate="%{y}<br>%{x}<br>%{z:.3g}<extra></extra>")
    col_groups === nothing && return Plot(hm, pq_layout(title;
        xaxis=(title=axis_title(""), tickangle=-45, type="category", automargin=true),
        yaxis=(title=axis_title(""), type="category", autorange="reversed", automargin=true)))
    grp = [get(col_groups, c, missing) for c in cols]
    lv = sort(unique(skipmissing(grp)))
    colors = series_colors(length(lv))
    idx = [ismissing(x) ? nothing : findfirst(==(x), lv) for x in grp]
    cs = length(colors) == 1 ? [[0, colors[1]], [1, colors[1]]] :
         [[(i - 1) / (length(colors) - 1), colors[i]] for i in eachindex(colors)]
    # 1 x n matrices (not vectors of vectors), serialized the same way as the main heatmap's z
    strip = heatmap(x=cols, y=["group"], z=permutedims(idx), showscale=false, xgap=1, colorscale=cs, zmin=1,
                    zmax=max(1, length(lv)), yaxis="y2", text=permutedims([coalesce(x, "") for x in grp]),
                    hovertemplate="%{x}<br>%{text}<extra></extra>")
    hm.fields[:yaxis] = "y"
    ann = [attr(xref="paper", yref="paper", x=1, y=1.02 + 0.045 * (length(lv) - i), xanchor="right",
                yanchor="bottom", showarrow=false, text="■ $(lv[i])", font=attr(color=colors[i], size=12))
           for i in eachindex(lv)]
    lay = pq_layout(title; annotations=ann, margin=attr(t=60 + 18 * length(lv)),
                    xaxis=(title=axis_title(""), tickangle=-45, type="category", automargin=true),
                    yaxis=(title=axis_title(""), type="category", autorange="reversed", automargin=true,
                           domain=[0, 0.94]))
    lay[:yaxis2] = attr(; theme_axis((domain=[0.955, 1.0], showticklabels=false, type="category"))...)
    Plot([strip, hm], lay)
end
plot_heatmap(m::AbstractMatrix; rows, cols, kwargs...) =
    plot_heatmap(LabeledMatrix(Matrix{Union{Missing,Float64}}(m), rows, cols); kwargs...)

"""
    plot_mutation_landscape(variants; n_genes=25, genes=nothing, title="Mutation landscape")

Genes (most frequently mutated first) by samples; cells show the most severe
impact (modifier < low < moderate < high), or simply mutated / not mutated
when impact isn't annotated.
"""
function plot_mutation_landscape(v::AbstractDataFrame; n_genes=25, genes=nothing, title="Mutation landscape")
    v = filter(r -> !ismissing(r.hgnc_symbol), v)
    levels = ["modifier", "low", "moderate", "high"]
    has_impact = "impact" in names(v) && any(!ismissing, v.impact)
    sev = has_impact ? [ismissing(i) ? 1.0 : Float64(something(findfirst(==(i), levels), 1)) for i in v.impact] :
          ones(nrow(v))
    w = DataFrame(sample_id=v.sample_id, hgnc_symbol=v.hgnc_symbol, severity=sev)
    freq = combine(groupby(w, :hgnc_symbol), :sample_id => (x -> length(unique(x))) => :n)
    fq = Dict(zip(freq.hgnc_symbol, freq.n))
    genes === nothing && (genes = first(sort(freq, [order(:n, rev=true), :hgnc_symbol]).hgnc_symbol, n_genes))
    w = filter(r -> r.hgnc_symbol in genes, w)
    m = to_matrix(w, "sample_id"; row="hgnc_symbol", value="severity", fun=maximum)
    ri = [findfirst(==(g), m.rows) for g in genes if g in m.rows]
    data = m.data[ri, :]; rows = m.rows[ri]
    pres = .!ismissing.(data)
    co = sortperm([Tuple(.!pres[:, j]) for j in 1:size(data, 2)])
    data = data[:, co]; cols = m.cols[co]
    s = plot_theme()["sequential"]
    lab = ["$g ($(fq[g]))" for g in rows]
    z = [ismissing(x) ? nothing : x for x in data]
    hm = if has_impact
        heatmap(x=cols, y=lab, z=z, text=[ismissing(x) ? "" : levels[Int(x)] for x in data],
                colorscale=[[0, s[2]], [0.33, s[3]], [0.34, s[4]], [0.66, s[5]], [0.67, s[6]], [1, s[7]]],
                zmin=1, zmax=4, xgap=1, ygap=1,
                colorbar=attr(title=attr(text="impact"), tickvals=1:4, ticktext=levels),
                hovertemplate="%{y}<br>%{x}<br>%{text}<extra></extra>")
    else
        heatmap(x=cols, y=lab, z=z, colorscale=[[0, s[5]], [1, s[5]]], showscale=false, xgap=1, ygap=1,
                hovertemplate="%{y}<br>%{x}<br>mutated<extra></extra>")
    end
    Plot(hm, pq_layout(title; xaxis=(title=axis_title("samples ($(length(cols)))"),
                                     showticklabels=length(cols) <= 40, type="category"),
                       yaxis=(title=axis_title(""), autorange="reversed", type="category", automargin=true)))
end

# -- analysis plots ------------------------------------------------------------

"""
    plot_zscores(comparison; n=30, min_value=1, min_observed=0.5, title=nothing)

Horizontal diverging bars of the top genes by |z| against the reference cohort.
"""
function plot_zscores(cmp::AbstractDataFrame; n=30, min_value=1, min_observed=0.5, title=nothing)
    top = sort(top_by_zscore(cmp; n=n, min_value=min_value, min_observed=min_observed), :z)
    d = plot_theme()["diverging"]
    if title === nothing && "patternq_comparison" in metadatakeys(cmp)
        info = metadata(cmp, "patternq_comparison")
        title = "$(info.sample) vs $(info.cohort_db): top genes by z-score"
    end
    txt = ["value $(round(v, sigdigits=3)) · cohort median $(round(2^m - 1, sigdigits=3)) · $(round(p, digits=1)) pct"
           for (v, m, p) in zip(top.value, top.cohort_median, top.percentile)]
    Plot(bar(x=top.z, y=top.hgnc_symbol, orientation="h", text=txt, textposition="none",
             marker=attr(color=[z >= 0 ? d["high"] : d["low"] for z in top.z]),
             hovertemplate="%{y}<br>z = %{x:.2f}<br>%{text}<extra></extra>"),
         pq_layout(title; showlegend=false, bargap=0.25,
                   xaxis=(title=axis_title("z-score vs cohort (log2(1+x))"), zeroline=true),
                   yaxis=(title=axis_title(""), type="category", categoryorder="array",
                          categoryarray=top.hgnc_symbol, automargin=true)))
end

"""
    plot_vs_cohort(sample_expr, cohort_expr; type=:violin, log=true, floor=0.01,
                   title="Samples vs cohort expression", xlab="expression")

One row per gene: the cohort's distribution (violin or box, one per `cohort`
column value) with each sample's value marked. On a log scale, log10 values are
drawn on a linear axis (so densities are estimated on the log scale) with
power-of-ten tick labels.
"""
function plot_vs_cohort(se::AbstractDataFrame, ce::AbstractDataFrame; type::Symbol=:violin, log::Bool=true,
                        floor=0.01, title="Samples vs cohort expression", xlab="expression")
    ce = "cohort" in names(ce) ? ce : transform(ce, [] => (() -> "cohort") => :cohort)
    ce = combine(groupby(ce, [:cohort, :sample_id, :hgnc_symbol]), :value => sum => :value)
    se = combine(groupby(se, [:sample_id, :hgnc_symbol]), :value => sum => :value)
    genes = unique(vcat(se.hgnc_symbol, ce.hgnc_symbol))
    fl(v) = log ? log10.(max.(v, floor)) : v
    th = plot_theme()
    cohorts = sort(unique(ce.cohort))
    traces = GenericTrace[]
    for (c, col) in zip(cohorts, series_colors(length(cohorts)))
        s = ce[ce.cohort .== c, :]
        push!(traces, type == :violin ?
            PlotlyBase.violin(orientation="h", x=fl(s.value), y=s.hgnc_symbol, name=c, legendgroup=c,
                              line=attr(color=col, width=1), fillcolor=alpha_color(col, 0.25), points=false,
                              spanmode="hard", scalemode="width", width=0.8,
                              box=attr(visible=true, fillcolor=th["surface"], line=attr(color=col), width=0.2),
                              meanline=attr(visible=false), hoverinfo="y+name") :
            box(orientation="h", x=fl(s.value), y=s.hgnc_symbol, name=c, legendgroup=c,
                marker=attr(color=col, size=3), line=attr(color=col), fillcolor=alpha_color(col, 0.2),
                boxpoints="outliers"))
    end
    symbols = ["diamond", "square", "circle", "triangle-up", "x", "star"]
    for (j, sid) in enumerate(sort(unique(se.sample_id)))
        s = se[se.sample_id .== sid, :]
        push!(traces, scatter(mode="markers", x=fl(s.value), y=s.hgnc_symbol, name=sid, customdata=s.value,
                              marker=attr(symbol=symbols[mod1(j, length(symbols))], size=11,
                                          color=th["text_primary"], line=attr(color=th["surface"], width=1.5)),
                              hovertemplate="$sid<br>%{y}: %{customdata:.3g}<extra></extra>"))
    end
    xaxis = Dict{Symbol,Any}(:title => axis_title(xlab))
    if log
        vals = vcat(fl(ce.value), fl(se.value))
        vals = vals[isfinite.(vals)]
        ticks = Base.floor(Int, minimum(vals)):ceil(Int, maximum(vals))
        xaxis[:tickvals] = collect(ticks)
        xaxis[:ticktext] = [t >= 3 ? replace(string(round(Int, 10.0^t)), r"(?<=\d)(?=(\d{3})+$)" => ",") :
                            (t >= 0 ? string(round(Int, 10.0^t)) : string(10.0^t)) for t in ticks]
    end
    lay = pq_layout(title; boxmode="group", violinmode="group", height=160 + 42 * length(genes),
                    yaxis=(title=axis_title(""), type="category", categoryorder="array",
                           categoryarray=reverse(genes), automargin=true))
    lay[:xaxis] = attr(; theme_axis(NamedTuple(xaxis))...)
    Plot(traces, lay)
end

"""
    plot_ma(change; lfc_threshold=2.5, highlight=nothing, label_top=8, label_min_avg=1.5, title=nothing)

MA plot: average expression (avg_log10) vs log2 fold change; points past the
threshold colored (up red, down blue), highlight genes and the top up/down
genes labelled (greedy thinning so labels don't pile up).
"""
function plot_ma(ch::AbstractDataFrame; lfc_threshold=2.5, highlight=nothing, label_top=8, label_min_avg=1.5,
                 title=nothing)
    d = plot_theme()["diverging"]; th = plot_theme()
    if title === nothing && "patternq_comparison" in metadatakeys(ch)
        info = metadata(ch, "patternq_comparison")
        haskey(info, :sample_a) && (title = "Expression change: $(info.sample_a) → $(info.sample_b)")
    end
    cls = [l >= lfc_threshold ? "up" : l <= -lfc_threshold ? "down" : "unchanged" for l in ch.lfc]
    spec = Dict("unchanged" => (alpha_color(OTHER_COLOR, 0.35), "within threshold", 5),
                "down" => (d["low"], "down (lfc ≤ -$lfc_threshold)", 7),
                "up" => (d["high"], "up (lfc ≥ $lfc_threshold)", 7))
    traces = GenericTrace[]
    for k in ("unchanged", "down", "up")
        s = ch[cls .== k, :]
        nrow(s) == 0 && continue
        col, nm, sz = spec[k]
        push!(traces, scattergl(mode="markers", x=s.avg_log10, y=s.lfc, text=s.hgnc_symbol, name=nm,
                                marker=attr(color=col, size=sz), customdata=hcat(s.value_a, s.value_b),
                                hovertemplate="%{text}<br>lfc %{y:.2f}<br>a %{customdata[0]:.3g} → b %{customdata[1]:.3g}<extra></extra>"))
    end
    lab = ch[ch.avg_log10 .>= label_min_avg, :]
    lab = vcat(first(sort(lab, :lfc, rev=true), label_top), first(sort(lab, :lfc), label_top))
    lab = lab[abs.(lab.lfc) .>= lfc_threshold, :]
    highlight === nothing || (lab = unique(vcat(ch[in.(ch.hgnc_symbol, Ref(Set(highlight))), :], lab)))
    xr = maximum(ch.avg_log10) - minimum(ch.avg_log10); yr = maximum(ch.lfc) - minimum(ch.lfc)
    kept = Int[]
    for i in 1:nrow(lab)
        any(j -> abs(lab.avg_log10[j] - lab.avg_log10[i]) < 0.06xr && abs(lab.lfc[j] - lab.lfc[i]) < 0.05yr, kept) ||
            push!(kept, i)
    end
    lab = lab[kept, :]
    ann = [attr(x=r.avg_log10, y=r.lfc, text=r.hgnc_symbol, showarrow=true, arrowhead=0, arrowwidth=1,
                arrowcolor=th["text_secondary"], ax=18, ay=r.lfc > 0 ? -16 : 16,
                font=attr(size=11, color=th["text_primary"])) for r in eachrow(lab)]
    guide(y) = attr(type="line", xref="paper", x0=0, x1=1, y0=y, y1=y,
                    line=attr(color=th["text_secondary"], width=1, dash=y == 0 ? "solid" : "dot"))
    Plot(traces, pq_layout(title; annotations=ann, hovermode="closest",
                           shapes=[guide(0), guide(lfc_threshold), guide(-lfc_threshold)],
                           xaxis=(title=axis_title("average expression, log10(1 + x)"),),
                           yaxis=(title=axis_title("log2 fold change"),)))
end

"""
    plot_fold_change(change; genes=nothing, title="Change in gene expression")

Fold change bars for selected genes (default: top 15 up and 15 down).
"""
function plot_fold_change(ch::AbstractDataFrame; genes=nothing, title="Change in gene expression")
    if genes === nothing
        s = sort(ch, :lfc)
        genes = unique(vcat(first(s.hgnc_symbol, 15), last(s.hgnc_symbol, 15)))
    end
    s = sort(ch[in.(ch.hgnc_symbol, Ref(Set(genes))), :], :lfc)
    d = plot_theme()["diverging"]
    Plot(bar(x=s.lfc, y=s.hgnc_symbol, orientation="h", marker=attr(color=[l >= 0 ? d["high"] : d["low"] for l in s.lfc]),
             hovertemplate="%{y}: %{x:.2f}<extra></extra>"),
         pq_layout(title; showlegend=false, height=140 + 22 * nrow(s), xaxis=(title=axis_title("log2 fold change"),),
                   yaxis=(title=axis_title(""), type="category", categoryorder="array", categoryarray=s.hgnc_symbol,
                          automargin=true)))
end
