# Survival helpers for outcome-association analyses (PRINCE-style: OS
# stratified at the median of a baseline biomarker, log-rank p-values,
# landmark survival status).

"""
    kaplan_meier(time, event)

Kaplan-Meier estimate: DataFrame of time, n_risk, n_event, n_censor, surv.
"""
function kaplan_meier(time, event)
    t = Float64.(collect(time)); e = Bool.(collect(event))
    o = sortperm(t); t = t[o]; e = e[o]
    ut = unique(t)
    n_risk = [count(>=(u), t) for u in ut]
    n_event = [count(i -> t[i] == u && e[i], eachindex(t)) for u in ut]
    n_censor = [count(i -> t[i] == u && !e[i], eachindex(t)) for u in ut]
    surv = cumprod(1 .- n_event ./ n_risk)
    DataFrame(time=ut, n_risk=n_risk, n_event=n_event, n_censor=n_censor, surv=surv)
end

# upper tail of the chi-square distribution: Q(k/2, x/2) (regularized gamma)
function chisq_pvalue(x::Real, df::Integer)
    x <= 0 && return 1.0
    a = df / 2; z = x / 2
    lg = loggamma_(a)
    if z < a + 1
        # series for P, then Q = 1 - P
        term = 1 / a; s = term; n = a
        for _ in 1:10_000
            n += 1; term *= z / n; s += term
            abs(term) < abs(s) * 1e-15 && break
        end
        return max(0.0, 1 - s * exp(-z + a * log(z) - lg))
    end
    # continued fraction for Q (Lentz)
    b = z + 1 - a; c = 1 / 1e-300; d = 1 / b; h = d
    for i in 1:10_000
        an = -i * (i - a); b += 2
        d = an * d + b; abs(d) < 1e-300 && (d = 1e-300)
        c = b + an / c; abs(c) < 1e-300 && (c = 1e-300)
        d = 1 / d; del = d * c; h *= del
        abs(del - 1) < 1e-15 && break
    end
    exp(-z + a * log(z) - lg) * h
end

# log gamma for positive half-integers / reals (Lanczos)
function loggamma_(x::Real)
    g = 7
    coef = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
            -176.61502916214059, 12.507343278686905, -0.13857109526572012,
            9.9843695780195716e-6, 1.5056327351493116e-7)
    x < 0.5 && return log(pi / abs(sin(pi * x))) - loggamma_(1 - x)
    x -= 1
    s = coef[1]
    for i in 1:g+1
        s += coef[i+1] / (x + i)
    end
    t = x + g + 0.5
    0.5 * log(2pi) + (x + 0.5) * log(t) - t + log(s)
end

"""
    logrank_test(time, event, group)

Mantel-Haenszel log-rank test for a difference in survival between two or more
groups (no Distributions.jl dependency). Matches R's survival::survdiff.
Returns a NamedTuple: chisq, df, p, observed, expected, n (Dicts by group).
"""
function logrank_test(time, event, group)
    keep = [!(ismissing(a) || ismissing(b) || ismissing(c)) for (a, b, c) in zip(time, event, group)]
    t = Float64.(collect(time)[keep]); e = Bool.(collect(event)[keep]); g = string.(collect(group)[keep])
    lv = sort(unique(g)); k = length(lv)
    k < 2 && return (chisq=NaN, df=0, p=NaN, observed=Dict{String,Float64}(), expected=Dict{String,Float64}(),
                     n=Dict{String,Int}())
    O = zeros(k); E = zeros(k); V = zeros(k, k)
    for tt in sort(unique(t[e]))
        at_risk = [count(i -> t[i] >= tt && g[i] == l, eachindex(t)) for l in lv]
        dg = [count(i -> t[i] == tt && e[i] && g[i] == l, eachindex(t)) for l in lv]
        n = sum(at_risk); d = sum(dg)
        n < 1 && continue
        O .+= dg
        E .+= d .* at_risk ./ n
        if n > 1
            f = d * (n - d) / (n^2 * (n - 1))
            V .+= f .* (Diagonal(at_risk .* n) .- at_risk * at_risk')
        end
    end
    idx = 1:k-1
    diff = (O .- E)[idx]
    chisq = diff' * (V[idx, idx] \ diff)
    (chisq=chisq, df=k - 1, p=chisq_pvalue(chisq, k - 1),
     observed=Dict(zip(lv, O)), expected=Dict(zip(lv, E)),
     n=Dict(l => count(==(l), g) for l in lv))
end

"""
    median_split(x; labels=("low", "high"))

"high" at or above the median, "low" below, `missing` where x is missing.
"""
function median_split(x; labels=("low", "high"))
    vals = collect(skipmissing(x))
    isempty(vals) && return Vector{Union{Missing,String}}(missing, length(x))
    m = median(vals)
    [ismissing(v) ? missing : (v >= m ? labels[2] : labels[1]) for v in x]
end

"""
    survival_status(outcomes; time="os", event="os_event", at=12)

Landmark status: "alive at 1 year" when time >= at, "died within 1 year" when
the event happened before `at`, `missing` when censored before `at` (unknown).
"""
function survival_status(outcomes; time="os", event="os_event", at=12,
                         labels=(alive="alive at 1 year", died="died within 1 year"))
    [ismissing(t) || t === nothing ? missing :
     (t >= at ? labels.alive : (!ismissing(e) && e === true ? labels.died : missing))
     for (t, e) in zip(outcomes[!, time], outcomes[!, event])]
end

"""
    survival_by_median(values, outcomes; time="os", event="os_event")

Split subjects at the median of `values` (a Dict subject_id => value, e.g. a
baseline cell-population frequency) and test the difference in survival.
Returns `outcomes` restricted to those subjects plus `value` and `group`
("low"/"high"); the log-rank test is in the table metadata ("logrank").
"""
function survival_by_median(values::AbstractDict, outcomes::AbstractDataFrame; time="os", event="os_event")
    vals = Dict(k => v for (k, v) in values if !ismissing(v))
    tab = filter(r -> haskey(vals, r.subject_id), outcomes)
    tab.value = [vals[s] for s in tab.subject_id]
    tab.group = median_split(tab.value)
    lr = logrank_test(tab[!, time], tab[!, event], tab.group)
    metadata!(tab, "logrank", lr; style=:note)
    tab
end

"""
    change_from_baseline(tab; baseline="C1D1", by=nothing, method=:log2_ratio, pseudocount=0)

Per subject (and target), the change of `value` relative to the subject's value
at the baseline timepoint: `:log2_ratio` (frequencies), `:difference` (values
already on a log scale, like Olink NPX) or `:ratio`. Adds baseline_value and
change; subjects without a baseline value are dropped.
"""
function change_from_baseline(tab::AbstractDataFrame; baseline="C1D1", by=nothing,
                              method::Symbol=:log2_ratio, pseudocount=0)
    by === nothing && (by = intersect(["cell_population", "epitope_id", "hgnc_symbol", "signature",
                                       "measurement_set", "target"], names(tab)))
    keycols = vcat(["subject_id"], by)
    key(r) = Tuple(r[c] for c in keycols)
    base = Dict{Any,Vector{Float64}}()
    for r in eachrow(tab)
        (r.timepoint_id == baseline && !ismissing(r.value)) || continue
        push!(get!(base, key(r), Float64[]), r.value)
    end
    out = filter(r -> haskey(base, key(r)), tab)
    out.baseline_value = [mean(base[key(r)]) for r in eachrow(out)]
    ch = if method == :log2_ratio
        log2.((out.value .+ pseudocount) ./ (out.baseline_value .+ pseudocount))
    elseif method == :difference
        out.value .- out.baseline_value
    elseif method == :ratio
        (out.value .+ pseudocount) ./ (out.baseline_value .+ pseudocount)
    else
        throw(ArgumentError("method must be :log2_ratio, :difference or :ratio"))
    end
    out.change = [ismissing(x) || !isfinite(x) ? missing : x for x in ch]
    keep_provenance(out, tab)
end
