# Timepoints, clinical observations, adverse events, clinical interventions,
# and the cell population / TCR / microbiome entities of measurement sets.

"""
    timepoints(; db=nothing)

timepoint_id, timepoint_relative_order, timepoint_type, offset, cycle/day where
present, ordered by relative order.
"""
function timepoints(; db=nothing, kwargs...)
    df = runq(dq([["pull", "?t", Any["*", EnumRef(":timepoint/type")]]], [["?t", ":timepoint/id"]]); db=db, kwargs...)
    "timepoint_relative_order" in names(df) && (df = keep_provenance(sort(df, :timepoint_relative_order), df))
    order_columns(df, ["timepoint_id", "timepoint_relative_order", "timepoint_type"])
end

"""
    clinical_observation_sets(; db=nothing)
"""
clinical_observation_sets(; db=nothing, kwargs...) =
    runq(dq([["pull", "?c", [":clinical-observation-set/name", ":clinical-observation-set/description"]]],
            [["?c", ":clinical-observation-set/name"]]); db=db, kwargs...)

"""
    clinical_observations(; db=nothing, obs_type=nothing, set_name=nothing, subjects=nothing)

With `obs_type` (e.g. "os", "pfs", "bor", "recist"): subject_id, timepoint_id
and a column named after the type (enums as names). With `set_name`: one row
per observation in the set with all its attributes.
"""
function clinical_observations(; db=nothing, obs_type=nothing, set_name=nothing, subjects=nothing, kwargs...)
    obs_type === nothing && set_name === nothing && error("Give obs_type or set_name")
    if obs_type !== nothing
        t = replace(lstrip(string(obs_type), ':'), r"^clinical-observation/" => "")
        attr = ":clinical-observation/" * t
        ins = Any[]; args = Any[]
        subjects === nothing || (push!(ins, ["?subject-id", "..."]); push!(args, collect(subjects)))
        df = runq(dq(Any["?subject-id", ["pull", "?o", Any[Dict(":clinical-observation/timepoint" => [":timepoint/id"]),
                                                          Dict(":clinical-observation/study-day" => [":study-day/id"]),
                                                          attr]]],
                     [["?o", attr], ["?o", ":clinical-observation/subject", "?p"], ["?p", ":subject/id", "?subject-id"]];
                     in=ins, args=args); db=db, kwargs...)
        rename_cols!(df, clean_name(attr) => clean_name(t))
        return df
    end
    df = runq(dq([["pull", "?o", Any["*",
                    Dict(":clinical-observation/subject" => [":subject/id"]),
                    Dict(":clinical-observation/timepoint" => [":timepoint/id"]),
                    Dict(":clinical-observation/study-day" => [":study-day/id"]),
                    EnumRef(":clinical-observation/recist"), EnumRef(":clinical-observation/bor"),
                    EnumRef(":clinical-observation/pfs-reason"), EnumRef(":clinical-observation/os-reason"),
                    EnumRef(":clinical-observation/disease-stage"),
                    Dict(":clinical-observation/metastasis-gdc-anatomic-sites" => [":gdc-anatomic-site/name"])]]],
                 [["?cos", ":clinical-observation-set/name", "?set-name"],
                  ["?cos", ":clinical-observation-set/clinical-observations", "?o"]];
                 in=["?set-name"], args=[set_name]); db=db, kwargs...)
    for n in names(df)
        startswith(n, "clinical_observation_") && rename!(df, n => replace(n, r"^clinical_observation_" => ""))
    end
    df
end

"""
    subject_outcomes(; db=nothing)

One row per subject with bor (from :clinical-observation/bor, or derived from
RECIST observations when absent: CR > PR > SD > PD), pfs, pfs_event, os,
os_event where present. Errors if a subject has several values of one outcome.
"""
function subject_outcomes(; db=nothing, kwargs...)
    db = ensure_db(db)
    function get1(type)
        r = try
            clinical_observations(; db=db, obs_type=type, kwargs...)
        catch
            nothing
        end
        (r === nothing || nrow(r) == 0) && return nothing
        col = clean_name(type)
        col in names(r) || return nothing
        r = dropmissing(select(r, "subject_id", col), col)
        length(unique(r.subject_id)) < nrow(r) && error("More than one $type value for some subjects")
        r
    end
    bor = get1("bor")
    if bor === nothing
        recist = try clinical_observations(; db=db, obs_type="recist", kwargs...) catch; nothing end
        if recist !== nothing && nrow(recist) > 0
            rank = Dict("CR" => 1, "PR" => 2, "SD" => 3, "PD" => 4)
            g = combine(groupby(recist, :subject_id), :recist => (x -> begin
                v = [r for r in skipmissing(x) if haskey(rank, r)]
                isempty(v) ? "Unknown" : v[argmin([rank[r] for r in v])]
            end) => :bor)
            bor = g
        end
    end
    ids = subjects(; db=db)
    out = select(ids, :subject_id)
    for p in (bor, get1("pfs"), get1("pfs-event"), get1("os"), get1("os-event"))
        p === nothing && continue
        out = leftjoin(out, p; on=:subject_id)
    end
    keep_provenance(sort(out, :subject_id), ids)
end

"""
    adverse_events(; db=nothing, set_name=nothing)
"""
function adverse_events(; db=nothing, set_name=nothing, kwargs...)
    where = Any[["?o", ":adverse-event/subject"]]; ins = Any[]; args = Any[]
    if set_name !== nothing
        where = Any[["?cos", ":clinical-observation-set/name", "?set-name"],
                    ["?cos", ":clinical-observation-set/adverse-events", "?o"]]
        ins = Any["?set-name"]; args = Any[set_name]
    end
    runq(dq([["pull", "?o", Any["*", Dict(":adverse-event/subject" => [":subject/id"]),
                                Dict(":adverse-event/timepoint" => [":timepoint/id"]),
                                Dict(":adverse-event/meddra-adverse-event" => [":meddra-disease/preferred-name"]),
                                EnumRef(":adverse-event/ctcae-grade"), EnumRef(":adverse-event/ae-causality"),
                                Dict(":adverse-event/study-day" => [":study-day/id"])]]],
            where; in=ins, args=args); db=db, kwargs...)
end

"""
    clinical_interventions(; db=nothing, subjects=nothing)

One row per intervention (treatments, surgeries, biopsies, ...).
"""
function clinical_interventions(; db=nothing, subjects=nothing, kwargs...)
    ins = Any[]; args = Any[]
    subjects === nothing || (push!(ins, ["?subject-id", "..."]); push!(args, collect(subjects)))
    df = runq(dq([["pull", "?ci", Any["*",
                    Dict(":clinical-intervention/subject" => [":subject/id"]),
                    Dict(":clinical-intervention/timepoint" => [":timepoint/id", ":timepoint/relative-order"]),
                    Dict(":clinical-intervention/treatment-regimen" => Any[":treatment-regimen/name",
                        Dict(":treatment-regimen/drug-regimens" => Any[Dict(":drug-regimen/drug" => [":drug/preferred-name"]),
                                                                       ":drug-regimen/freetext-drug"])]),
                    EnumRef(":clinical-intervention/surgery-type"),
                    EnumRef(":clinical-intervention/cancer-medication-category"),
                    EnumRef(":clinical-intervention/radiation-therapy-category"),
                    EnumRef(":clinical-intervention/biospecimen-type"),
                    EnumRef(":clinical-intervention/biospecimen-collection"),
                    Dict(":clinical-intervention/biospecimen-derived-samples" => [":sample/id"])]]],
                 [["?ci", ":clinical-intervention/subject", "?p"], ["?p", ":subject/id", "?subject-id"]];
                 in=ins, args=args); db=db, kwargs...)
    for n in names(df)
        startswith(n, "clinical_intervention_") && rename!(df, n => replace(n, r"^clinical_intervention_" => ""))
    end
    df
end

function ms_entities(ref_attr, pattern, measurement_set, db; kwargs...)
    runq(dq([["pull", "?e", pattern]], [["?ms", ":measurement-set/name", "?ms-name"], ["?ms", ref_attr, "?e"]];
            in=["?ms-name"], args=[measurement_set]); db=db, kwargs...)
end

"""
    cell_populations(measurement_set; db=nothing)
    tcrs(measurement_set; db=nothing)
    otus(measurement_set; db=nothing)
    sgbs(measurement_set; db=nothing)

Entities attached to a measurement set.
"""
cell_populations(ms::AbstractString; db=nothing, kwargs...) =
    ms_entities(":measurement-set/cell-populations",
                Any["*", Dict(":cell-population/cell-type" => [":cell-type/co-name"]),
                    Dict(":cell-population/positive-markers" => [":epitope/id"]),
                    Dict(":cell-population/negative-markers" => [":epitope/id"]),
                    Dict(":cell-population/parent" => [":cell-population/name"])], ms, db; kwargs...)
tcrs(ms::AbstractString; db=nothing, kwargs...) = ms_entities(":measurement-set/tcrs", ["*"], ms, db; kwargs...)
otus(ms::AbstractString; db=nothing, kwargs...) = ms_entities(":measurement-set/otus", ["*"], ms, db; kwargs...)
sgbs(ms::AbstractString; db=nothing, kwargs...) = ms_entities(":measurement-set/sgbs", ["*"], ms, db; kwargs...)
