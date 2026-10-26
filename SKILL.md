---
name: patternq-julia
description: Query the Pattern Data Commons (PDC, RCRF datasets) from Julia with PatternQ.jl. Use when asked about Pattern Data Commons datasets, subjects, samples, measurements, variants, copy number, gene expression, or clinical outcomes, or when Julia code uses PatternQ.
---

# Pattern Data Commons from Julia (PatternQ.jl)

PatternQ queries the Pattern Data Commons query service over HTTPS. You act on
behalf of the user whose API key is configured. The key can see only that
user's datasets, the data is de-identified, and every call is read-only.

## Setup

```julia
using Pkg
Pkg.add(url="https://github.com/RCRF/PatternQ.jl")
using PatternQ, DataFrames
```

```sh
export PATTERNQ_ENDPOINT=https://data-commons.rcrf-dev.org   # the default
export PATTERNQ_API_KEY=...   # the user's key, from the dashboard's user settings
```

The user supplies the key. Read it from the environment and never print it,
log it, write it into scripts or notebooks, or commit it.

## Every dataset is its own database

Resolve the dataset name to its current database first. Database names
change when a dataset is re-imported, so don't hard-code them.

```julia
list_datasets()                          # dataset, db, patient_count, sample_count, assays, tags
db = resolve_db("tcga-uvm")              # dataset name -> current database name
set_db(db)                               # session default, or pass db= to any call
```

To cover several datasets, loop over `list_datasets()` and pass `db=` each
time. Add a `dataset` column yourself before `vcat(...; cols=:union)`.

## The calls you need most

| task | call |
|---|---|
| what's in a dataset | `dataset_summary()`, `measurement_sets()`, `measurement_types(set_name)` |
| subjects (demographics, disease) | `subjects()` |
| samples (subject, timepoint, type) | `samples()`, `timepoints()`, `sample_assays()` |
| measurement values | `measurements("tpm", set_name; samples=..., targets=...)` |
| somatic variants | `variants(; genes=[...], samples=[...])` |
| gene expression | `gene_expression(; genes=[...], measurement="tpm")` |
| copy number | `cnv_gene_calls(; genes=[...])`, `cnv_segments(; genes=[...])` (genes, samples, or subjects required) |
| clinical | `subject_outcomes()` (BOR, PFS, OS), `clinical_observation_sets()`, `clinical_observations(; obs_type="os")` |
| sample or subject context | `add_sample_context(df)`, `add_subject_context(df)` |
| large matrices | `measurement_matrices()`, `measurement_matrix_by_name(name)` |

Results are DataFrames with snake_case columns, and enum values come back as
names. `provenance(df)` gives the database, basis t, and query time (use
`keep_provenance` after your own joins). Cite these when you report numbers.

## Raw queries are the fallback

Use the catalog first. When nothing in it fits, write the JSON form of
Datomic Datalog as Julia data:

```julia
do_query(Dict(":find" => ["?id", "?sex"],
              ":where" => [["?s", ":subject/id", "?id"],
                           ["?s", ":subject/sex", "?x"],
                           ["?x", ":db/ident", "?sex"]]))
```

Every catalog function has a `*_query` companion that returns `(query, args)`,
a good starting point to adapt. `schema_info()` gives the schema version.
Bind an attribute or entity in every clause: the server refuses full scans
(`query-too-broad`).

## Rate limits

Each API key gets about 600 requests a minute and 4 queries in flight.
PatternQ handles HTTP 429 and 503 itself: it honors Retry-After, retries up
to 5 times, and caps concurrent queries at 4 across tasks. Don't add your own
retry loop. `ThrottledError` means the limit outlasted the retries, so slow
down. A server-side query timeout isn't retried; narrow or page the query
instead (the server caps timeouts at 120 s). Prefer a few broad catalog
calls over many tiny ones.

## More for agents

- Agent guide: https://data-commons.rcrf-dev.org/docs/agents (also
  https://data-commons.rcrf-dev.org/llms.txt).
- Without a client library, the commons also serves MCP at
  https://data-commons.rcrf-dev.org/mcp (Streamable HTTP, `Authorization:
  Bearer <API key>`), with `list_datasets`, `describe_dataset`, and `query`
  (EDN Datalog, at most 1,000 rows, no matrices). For larger results,
  matrices, and the analysis functions, use this library.
