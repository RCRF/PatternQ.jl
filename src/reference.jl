# Reference data queries. Reference entities (genes, gene products, proteins,
# epitopes, cell types, ...) are included in each dataset database.

"""
    gene_symbols(; db=nothing)

Vector of HGNC symbols.
"""
gene_symbols(; db=nothing, kwargs...) =
    runq(dq(["?hgnc-symbol"], [["_", ":gene/hgnc-symbol", "?hgnc-symbol"]]); db=db, kwargs...).hgnc_symbol

"""
    genes(; db=nothing)

One row per gene: gene_hgnc_symbol, ids, previous/alias symbols (vectors).
"""
genes(; db=nothing, kwargs...) =
    runq(dq([["pull", "?g", Any[":gene/hgnc-symbol", ":gene/hgnc-id", ":gene/hgnc-name", ":gene/ensembl-id",
                               ":gene/previous-hgnc-symbols", ":gene/alias-hgnc-symbols",
                               EnumRef(":gene/hgnc-locus-group")]]],
            [["?g", ":gene/hgnc-symbol"]]); db=db, kwargs...)

"""
    gene_products(; db=nothing)
"""
gene_products(; db=nothing, kwargs...) =
    runq(dq(["?gene-product-id", "?hgnc-symbol"],
            [["?gp", ":gene-product/id", "?gene-product-id"], ["?gp", ":gene-product/gene", "?g"],
             ["?g", ":gene/hgnc-symbol", "?hgnc-symbol"]]); db=db, kwargs...)

"""
    gene_coordinates(; db=nothing, genes=nothing)

hgnc_symbol, assembly, contig, strand, start, end.
"""
function gene_coordinates(; db=nothing, genes=nothing, kwargs...)
    ins = Any[]; args = Any[]
    genes === nothing || (push!(ins, ["?hgnc-symbol", "..."]); push!(args, collect(genes)))
    df = runq(dq(["?hgnc-symbol", "?assembly", "?contig", "?strand", "?start", "?end"],
                 [["?g", ":gene/hgnc-symbol", "?hgnc-symbol"], ["?g", ":gene/genomic-coordinates", "?gc"],
                  ["?gc", ":genomic-coordinate/assembly", "?a"], ["?a", ":db/ident", "?assembly"],
                  ["?gc", ":genomic-coordinate/contig", "?contig"], ["?gc", ":genomic-coordinate/strand", "?strand"],
                  ["?gc", ":genomic-coordinate/start", "?start"], ["?gc", ":genomic-coordinate/end", "?end"]];
                 in=ins, args=args); db=db, kwargs...)
    nrow(df) > 0 && (df.assembly = ident_name.(df.assembly))
    df
end

"""
    variant_annotations(; db=nothing, variant_ids=nothing, genes=nothing)

Variant reference entities: id, gene, HGVS, impact, consequences, classification.
"""
function variant_annotations(; db=nothing, variant_ids=nothing, genes=nothing, kwargs...)
    where = Any[["?v", ":variant/id", "?variant-id"]]; ins = Any[]; args = Any[]
    variant_ids === nothing || (push!(ins, ["?variant-id", "..."]); push!(args, collect(variant_ids)))
    if genes !== nothing
        append!(where, [["?v", ":variant/gene", "?g"], ["?g", ":gene/hgnc-symbol", "?gene"]])
        push!(ins, ["?gene", "..."]); push!(args, collect(genes))
    end
    df = runq(dq([["pull", "?v", Any[":variant/id", ":variant/HGVSp", ":variant/HGVSc", ":variant/ref-allele",
                                   ":variant/alt-allele", ":variant/coordinate-string", ":variant/dbSNP",
                                   ":variant/max-af", ":variant/external-ids",
                                   Dict(":variant/gene" => [":gene/hgnc-symbol"]), EnumRef(":variant/impact"),
                                   EnumRef(":variant/classification"), EnumRef(":variant/type"),
                                   Dict(":variant/so-consequences" => [":so-sequence-feature/name"])]]],
                 where; in=ins, args=args); db=db, kwargs...)
    for c in ("variant_HGVSp", "variant_HGVSc", "variant_so_consequences", "variant_external_ids")
        c in names(df) && (df[!, c] = tidy_column(Any[join_many(v) for v in df[!, c]]))
    end
    for n in names(df)
        startswith(n, "variant_") && rename!(df, n => replace(n, r"^variant_" => ""))
    end
    rename_cols!(df, "gene_hgnc_symbol" => "hgnc_symbol", "id" => "variant_id")
    order_columns(df, ["variant_id", "hgnc_symbol", "HGVSp", "HGVSc", "impact"])
end

"""
    cnvs(; db=nothing)

CNV reference entities: cnv_id, coordinates, genes (vector).
"""
cnvs(; db=nothing, kwargs...) =
    runq(dq([["pull", "?c", Any[":cnv/id",
                               Dict(":cnv/genomic-coordinates" => [":genomic-coordinate/contig",
                                                                  ":genomic-coordinate/start",
                                                                  ":genomic-coordinate/end"]),
                               Dict(":cnv/genes" => [":gene/hgnc-symbol"])]]], [["?c", ":cnv/id"]]); db=db, kwargs...)

simple_names(attr; db=nothing, kwargs...) =
    runq(dq(["?name"], [["_", attr, "?name"]]); db=db, kwargs...).name

"Names of reference entities present in a database."
gdc_anatomic_sites(; db=nothing, kwargs...) = simple_names(":gdc-anatomic-site/name"; db=db, kwargs...)
epitopes(; db=nothing, kwargs...) = simple_names(":epitope/id"; db=db, kwargs...)
cell_types(; db=nothing, kwargs...) = simple_names(":cell-type/co-name"; db=db, kwargs...)
meddra_diseases(; db=nothing, kwargs...) = simple_names(":meddra-disease/preferred-name"; db=db, kwargs...)
drugs(; db=nothing, kwargs...) = simple_names(":drug/preferred-name"; db=db, kwargs...)
proteins(; db=nothing, kwargs...) =
    runq(dq([["pull", "?p", Any[":protein/preferred-name", ":protein/uniprot-name", ":protein/uniprot-accessions",
                               Dict(":protein/gene" => [":gene/hgnc-symbol"])]]],
            [["?p", ":protein/preferred-name"]]); db=db, kwargs...)

"""
    db_idents(; db=nothing)
"""
db_idents(; db=nothing, kwargs...) =
    runq(dq(["?db-id", "?db-ident"], [["?db-id", ":db/ident", "?db-ident"]]); db=db, kwargs...)

"""
    map_gene_symbols(symbols; all_genes=nothing, db=nothing, warn=true)

Resolve previous and alias symbols to current HGNC symbols, case
insensitively: a Dict input symbol => HGNC symbol (`missing` if unmapped).
"""
function map_gene_symbols(symbols; all_genes=nothing, db=nothing, warn::Bool=true)
    all_genes === nothing && (all_genes = genes(; db=db))
    lookup = Dict(uppercase(s) => s for s in all_genes.gene_hgnc_symbol)
    for col in ("gene_previous_hgnc_symbols", "gene_alias_hgnc_symbols")
        col in names(all_genes) || continue
        for (sym, v) in zip(all_genes.gene_hgnc_symbol, all_genes[!, col])
            ismissing(v) && continue
            for a in (v isa AbstractVector ? v : [v])
                k = uppercase(string(a))
                haskey(lookup, k) || (lookup[k] = sym)   # current symbols take precedence
            end
        end
    end
    out = Dict(s => get(lookup, uppercase(s), missing) for s in symbols)
    nmiss = count(ismissing, values(out))
    warn && nmiss > 0 && @warn "$nmiss of $(length(out)) gene symbols could not be mapped"
    out
end
