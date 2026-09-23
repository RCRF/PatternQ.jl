"""
    PatternQ

Query and analysis tools for the Pattern Data Commons: Datalog queries over the
query service, canned queries for every kind of data in a dataset database,
survival and expression analysis, and plotly plots.

Configure access with the `PATTERNQ_ENDPOINT` / `PATTERNQ_API_KEY` environment
variables, or `set_query_server` / `set_token`. Every dataset is its own
database: select one with `set_db` or pass `db=` to any query function.

Part of the patternq family (Python, R, Clojure, Julia) sharing one function
catalog, the same result columns and the same plots.
"""
module PatternQ

using CodecZlib
using DataFrames
using Dates
using HTTP
using JSON3
using LinearAlgebra
using OrderedCollections
using PlotlyBase
using Statistics

include("config.jl")
include("results.jl")
include("transport.jl")
include("dataset.jl")
include("measurements.jl")
include("clinical.jl")
include("reference.jl")
include("context.jl")
include("survival.jl")
include("analysis.jl")
include("plots.jl")

# config / transport
export query_server, set_query_server, set_token, set_db, current_db,
       query, query_body, do_query, across_dbs, datoms, list_datasets, resolve_db,
       measurement_matrix, measurement_matrix_by_name
# results
export provenance, keep_provenance, flatten_pull, clean_name, ident_name, join_many
# dataset / measurements / clinical / reference
export samples, samples_query, subjects, subjects_query, dataset_summary, dataset_summary_query,
       dataset_info, schema_info, measurement_sets, measurement_types, measurement_set_attributes,
       sample_assays, measurement_matrices, measurements, measurements_query, variants, variants_query,
       gene_expression, gene_expression_query, isoforms, cnv_segments, cnv_gene_calls,
       timepoints, clinical_observation_sets, clinical_observations, subject_outcomes,
       adverse_events, clinical_interventions, cell_populations, tcrs, otus, sgbs,
       gene_symbols, genes, gene_products, gene_coordinates, variant_annotations, cnvs,
       gdc_anatomic_sites, proteins, epitopes, cell_types, meddra_diseases, drugs, db_idents,
       map_gene_symbols
# context
export to_matrix, to_long, split_by_measurement_set, select_targets, add_sample_context,
       add_subject_context, add_variant_context, add_cnv_context, deduplicate_taxonomy, aggregate_taxa
# survival
export kaplan_meier, logrank_test, median_split, survival_status, survival_by_median,
       change_from_baseline
# analysis
export genesets, geneset, sample_expression, percentile_rank, compare_to_cohort, top_by_zscore,
       log_fold_change, compare_samples, ssgsea_score, ssgsea, top_varying_genes,
       expression_distance, nearest_samples, examine_geneset
# plots
export plot_theme, plot_vaf_histogram, plot_gene_expression, plot_sample_overview, plot_by_timepoint,
       plot_by_group, plot_survival, plot_heatmap, plot_mutation_landscape, plot_zscores,
       plot_vs_cohort, plot_ma, plot_fold_change

end # module
