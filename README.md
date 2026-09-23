# PatternQ.jl

Query and analysis tools for the Pattern Data Commons: Datalog queries over the
query service, canned queries for every kind of data in a dataset database,
survival and expression analysis, and plotly plots. Part of the patternq family
([patternq](https://github.com/RCRF/patternq) for Python,
[patternq-r](https://github.com/RCRF/patternq-r),
[patternq-clj](https://github.com/RCRF/patternq-clj) and
[PatternQ.jl](https://github.com/RCRF/PatternQ.jl)), which share one function
catalog, the same result columns and the same plots.

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/RCRF/PatternQ.jl")
```

Julia 1.10 or newer. Dependencies: HTTP, JSON3, CodecZlib, DataFrames,
OrderedCollections and PlotlyBase.

## Configure

```sh
export PATTERNQ_ENDPOINT=https://data-commons.rcrf-dev.org
export PATTERNQ_API_KEY=...   # user settings page of the Pattern Data Commons dashboard
```

or in a session: `set_query_server(url)`, `set_token(token)`.

## Quick start

Every dataset is its own database. Resolve a dataset name to its current
database and set it as the session default, or pass `db=` to any function.

```julia
using PatternQ, DataFrames

list_datasets()                          # dataset, db, counts, assays, tags
set_db(resolve_db("prince-2022"))

dataset_summary()                        # assays and measurement sets
measurement_types("PICI CyTOF Immune Profiling")
subjects(); samples(); timepoints()
outcomes = subject_outcomes()            # BOR / PFS / OS per subject

cy = measurements("percent-of-parent", "PICI CyTOF Immune Profiling")
cy = add_sample_context(cy; include_outcomes=true)
plot_by_timepoint(cy[cy.cell_population .== "HLA-DR+ Non-Naive CD8 T Cells", :]; group="bor")

plot_survival(outcomes; group="bor")
plot_mutation_landscape(variants(); n_genes=20)

# CNV queries always take a subset: genes, samples or subjects
cnv_segments(db=resolve_db("H37004"), genes=["TP53", "BAP1"])
```

### Survival by a baseline biomarker

```julia
base = cy[(cy.cell_population .== "HLA-DR+ Non-Naive CD8 T Cells") .& (coalesce.(cy.timepoint_id, "") .== "C1D1"), :]
values = Dict(k.subject_id => sum(g.value) / nrow(g) for (k, g) in pairs(groupby(base, :subject_id)))
km = survival_by_median(values, outcomes)
metadata(km, "logrank").p
plot_survival(km; group="group", levels=["low", "high"])
```

### A sample against a reference cohort

```julia
h = resolve_db("H37001"); uvm = resolve_db("tcga-uvm")
cmp = compare_to_cohort("H37001-003"; db=h, cohort_db=uvm)
top_by_zscore(cmp; n=25)
plot_zscores(cmp)
examine_geneset(["BAP1", "GNAQ", "PRAME", "PMEL"], ["H37001-003"]; db=h, cohort_dbs=uvm)

ch = compare_samples("H37001-003", "H37001-001"; db=h)   # two-sample change
plot_ma(ch); plot_fold_change(ch)
```

Only compare measurements that are comparable across the two databases (same
units and normalization; TPM is the usual common ground), and check
`cohort_observed`: a usually-expressed gene with few stored cohort values
points to an import or annotation problem in the cohort.

### Writing queries

Queries are the JSON form of Datomic Datalog, written as Julia data:

```julia
q = Dict(":find" => ["?sample-id", "?vaf"],
         ":in" => [["?hgnc", "..."]],
         ":where" => [["?g", ":gene/hgnc-symbol", "?hgnc"], ["?v", ":variant/gene", "?g"],
                      ["?m", ":measurement/variant", "?v"], ["?m", ":measurement/vaf", "?vaf"],
                      ["?m", ":measurement/sample", "?s"], ["?s", ":sample/id", "?sample-id"]])
do_query(q; args=[["TP53", "KRAS"]])
```

Results are DataFrames whose columns come from the Datalog variables and
attributes (`?sample-id`, `:sample/id` -> `sample_id`); enum values come back as
names. Every result carries `provenance(df)`: database, basis t and query time
(DataFrame metadata; patternq's context joins keep it, use `keep_provenance`
after your own joins). Every canned query has a `*_query` companion returning
`(query, args)`.

`cache=true` (the default) uses the service's S3 result cache; `cache=false`
returns results inline and skips it; `refresh_cache=true` recomputes.

Wide data (samples x targets) is a `LabeledMatrix` (`data`, `rows`, `cols`),
from `to_matrix` or `measurements(...; wide=true)`.

Plots are `PlotlyBase.Plot` objects; they display as HTML in Quarto, Jupyter,
Pluto and VS Code.

## Tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'   # live tests run when PATTERNQ_API_KEY is set
```

## Examples

`examples/` has a tutorial and a PRINCE trial demo as Quarto notebooks (Julia
engine); the same notebooks exist for every language in the family.

## Contributing

patternq is one library in four languages: [patternq](https://github.com/RCRF/patternq)
(Python), [patternq-r](https://github.com/RCRF/patternq-r) (R),
[patternq-clj](https://github.com/RCRF/patternq-clj) (Clojure) and
[PatternQ.jl](https://github.com/RCRF/PatternQ.jl) (Julia). All four are generated
from a common source and published here, so this repository does not accept pull
requests.

Please report bugs and feature requests as
[issues](https://github.com/RCRF/PatternQ.jl/issues). Code is welcome in an issue: a
minimal example (the call or query, the dataset, what you expected and what you
got), or a proposed change as a snippet, is the most useful way to suggest one.
