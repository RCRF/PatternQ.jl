using PatternQ, DataFrames, Statistics, Test

@testset "offline" begin
    @testset "results" begin
        m = Dict(":sample/id" => "S1", ":db/id" => 1, ":sample/subject" => Dict(":subject/id" => "P1"),
                 ":sample/specimen" => Dict(":db/ident" => ":sample.specimen/ffpe"))
        flat, _ = flatten_pull(m)
        @test Dict(flat) == Dict("sample_id" => "S1", "subject_id" => "P1", "sample_specimen" => "ffpe")
        @test PatternQ.find_column_names(["?sample-id", ["count", "?s"]]) == ["sample_id", "count_s"]
        @test clean_name(":subject/id") == "subject_id"
        @test ident_name(":variant.impact/high") == "high"
        q = PatternQ.normalize_query(Dict("find" => ["?id"], "where" => []))
        @test q[":in"] == ["\$"]
        @test join_many(["a", "b"]) == "a; b"
    end
    @testset "shared cross-language reference values" begin
        x = Dict("G$i" => Float64(101 - i) for i in 1:100)
        @test ssgsea_score(x, ["G$i" for i in 1:10]) ≈ 50.0216062883 atol = 1e-8
        @test ssgsea_score(x, ["G3", "G17", "G50", "G88"]) ≈ 17.3034501639 atol = 1e-8
        @test isnan(ssgsea_score(x, ["nope"]))
        lr = logrank_test(1:10, [true, true, false, true, true, true, false, true, false, true],
                          ["a", "b", "a", "a", "b", "a", "b", "b", "a", "b"])
        @test lr.chisq ≈ 0.2201257417 atol = 1e-8
        @test lr.df == 1
        @test PatternQ.chisq_pvalue(3.841458820694124, 1) ≈ 0.05 atol = 1e-10
        @test PatternQ.chisq_pvalue(5.991464547107979, 2) ≈ 0.05 atol = 1e-10
    end
    @testset "survival helpers" begin
        oc = DataFrame(os=[20, 5, 5, missing], os_event=[false, true, false, true])
        @test isequal(survival_status(oc), ["alive at 1 year", "died within 1 year", missing, missing])
        @test isequal(median_split([1, 2, 3, 4, missing]), ["low", "low", "high", "high", missing])
        km = kaplan_meier([1, 2, 3, 4], [true, false, true, true])
        @test km.surv ≈ [0.75, 0.75, 0.375, 0.0]
        tab = DataFrame(subject_id=["p1", "p1", "p2", "p2", "p3"], timepoint_id=["C1D1", "C2D1", "C1D1", "C2D1", "C2D1"],
                        cell_population="T", value=[1.0, 4, 2, 2, 9])
        ch = change_from_baseline(tab)
        @test nrow(ch) == 4
        @test ch.change[(ch.subject_id .== "p1") .& (ch.timepoint_id .== "C2D1")] == [2.0]
        d = change_from_baseline(tab; method=:difference)
        @test d.change[(d.subject_id .== "p1") .& (d.timepoint_id .== "C2D1")] == [3.0]
        oc2 = DataFrame(subject_id=["s$i" for i in 1:8], os=[1, 2, 3, 4, 10, 11, 12, 13], os_event=true)
        km2 = survival_by_median(Dict("s$i" => Float64(i) for i in 1:8), oc2)
        @test sort(unique(km2.group)) == ["high", "low"]
        @test metadata(km2, "logrank").p < 0.05
    end
    @testset "analysis helpers" begin
        @test percentile_rank(1:4, 2.5) == 50
        @test percentile_rank([1, 2, 2, 3], 2) == 50
        lfc = log_fold_change(Dict("G1" => 0.0, "G2" => 3.0, "G3" => 7.0), Dict("G1" => 1.0, "G2" => 3.0, "G4" => 15.0))
        get1(g) = only(lfc.lfc[lfc.hgnc_symbol .== g])
        @test get1("G1") == 1 && get1("G2") == 0 && get1("G3") == -3 && get1("G4") == 4
        @test expression_distance(Dict("A" => 1.0, "B" => 0.0), Dict("A" => 0.0, "B" => 1.0)) ≈ 1
        @test expression_distance(Dict("A" => 1.0, "B" => 0.0), Dict("A" => 0.0, "B" => 1.0); method=:euclidean) ≈ sqrt(2)
        gs = genesets()
        @test "GAPDH" in gs["housekeeping"] && haskey(gs, "hallmark_hypoxia")
        m = to_matrix(DataFrame(sample_id=["a", "a", "b"], g=["x", "y", "x"], value=[1.0, 2, 3]), "g")
        @test m["a", "y"] == 2.0 && ismissing(m["b", "y"])
        @test nrow(to_long(m)) == 3
        @test PatternQ.hclust_order([0.0 0; 10 10; 0.1 0; 10 10.2]) == [1, 3, 2, 4]
    end
    @testset "CNV subset rule" begin
        @test_throws ArgumentError cnv_segments(db="any")
        @test_throws ArgumentError cnv_gene_calls(db="any")
    end
    @testset "plots build offline" begin
        v = DataFrame(sample_id=["a", "a", "b"], vaf=[0.1, 0.5, 0.3])
        @test length(plot_vaf_histogram(v).data) == 2
        oc = DataFrame(os=[1, 2, 3, 4], os_event=[true, true, false, true], g=["x", "y", "x", "y"])
        @test plot_survival(oc; group="g") isa PatternQ.PlotlyBase.Plot
        m = to_matrix(DataFrame(sample_id=["a", "a", "b", "b", "c", "c"], g=["x", "y", "x", "y", "x", "y"],
                                value=[1.0, 2, 3, 1, 2, 5]), "g")
        p = plot_heatmap(transpose(m); scale=:row, col_groups=Dict("a" => "g1", "b" => "g2", "c" => "g1"))
        @test length(p.data) == 2
    end
end

if !isempty(get(ENV, "PATTERNQ_API_KEY", ""))
    dbs = Dict{String,String}()
    tdb(name) = get!(() -> resolve_db(name), dbs, name)
    @testset "live" begin
        @testset "H37001" begin
            db = tdb("H37001")
            @test nrow(samples(db=db)) == 6
            @test subjects(db=db).subject_id == ["H37001"]
            v = variants(db=db)
            @test nrow(v) > 100 && "BAP1" in v.hgnc_symbol
            @test names(v)[1:4] == ["sample_id", "measurement_set", "variant_id", "hgnc_symbol"]
            @test provenance(v).db == db
            gx = gene_expression(db=db, genes=["BAP1", "GNAQ"], measurement="rsem-normalized-count")
            @test Set(gx.hgnc_symbol) == Set(["BAP1", "GNAQ"])
        end
        @testset "tcga-uvm" begin
            db = tdb("tcga-uvm")
            @test nrow(subjects(db=db)) == 80
            @test Set(dataset_summary(db=db).assay_technology) == Set(["WES", "RNA-seq"])
            n = do_query(Dict(":find" => [["count", "?s"]], ":where" => [["?s", ":sample/id"]]); db=db, cache=false)
            @test n.count_s == [80]
            d = datoms("aevt", [":subject/id"]; db=db, limit=3)
            @test nrow(d) == 3 && d.a[1] == ":subject/id"
        end
        @testset "prince: measurements, clinical, survival" begin
            db = tdb("prince-2022")
            mt = measurement_types("PICI CyTOF Immune Profiling"; db=db)
            @test issubset(["percent-of-parent", "median-channel-value", "cell-population", "epitope"], mt.attribute)
            pp = measurements("percent-of-parent", "PICI CyTOF Immune Profiling"; db=db)
            @test names(pp) == ["sample_id", "measurement_set", "cell_population", "value"]
            @test nrow(pp) == only(mt.count[mt.attribute .== "percent-of-parent"])
            @test size(measurements("olink-npx", "PICI Olink Proteomics"; db=db, wide=true)) == (262, 30)
            oc = subject_outcomes(db=db)
            @test issubset(["subject_id", "bor", "os", "os_event", "pfs"], names(oc))
            @test nrow(oc) == nrow(subjects(db=db))
            x50 = add_sample_context(measurements("percent-of-parent", "PICI X50 Immune Profiling"; db=db);
                                     db=db, include_subjects=false)
            @test provenance(x50).db == db
            base(pop) = (b = x50[(x50.cell_population .== pop) .& (coalesce.(x50.timepoint_id, "") .== "C1D1"), :];
                         Dict(k.subject_id => mean(g.value) for (k, g) in pairs(groupby(b, :subject_id))))
            # cross-language check: same p as the R, Python and Clojure libraries
            km = survival_by_median(base("PD-1+Tbet+ non-Naive CD4+ T cells (% of CD4 not naive T cells)"), oc)
            @test metadata(km, "logrank").p ≈ 0.01152 atol = 5e-5
            km = survival_by_median(base("CD38+ effector memory CD8+ T cells (% effector memory CD8 T cells)"), oc)
            @test metadata(km, "logrank").p ≈ 0.7055 atol = 5e-4
            @test only(values(map_gene_symbols(["p53"]; db=db, warn=false))) == "TP53"
        end
        @testset "H37004: CNV, isoforms, interventions" begin
            db = tdb("H37004")
            cs = cnv_segments(db=db, genes=["TP53"])
            @test all(==("TP53"), cs.hgnc_symbol)
            @test issubset(["measurement_set", "contig", "start", "end", "segment_mean_lrr"], names(cs))
            @test nrow(cnv_segments(db=db, subjects=["H37004"])) > nrow(cs)
            @test "isoform_percent" in names(isoforms("TP53"; db=db))
            @test nrow(clinical_interventions(db=db)) > 0
            @test nrow(measurement_matrices(db=db)) == 2
        end
        @testset "painter: gene-level CNV calls" begin
            db = tdb("painter-2025-angiosarc")
            calls = cnv_gene_calls(db=db, genes=["MYC", "CDKN2A"])
            @test all(in(-2:2), calls.value)
            @test length(unique(calls.sample_id)) < nrow(samples(db=db))
        end
        @testset "cohort comparison and two-sample change" begin
            h = tdb("H37001"); uvm = tdb("tcga-uvm")
            cmp = compare_to_cohort("H37001-003"; db=h, cohort_db=uvm, genes=["MLANA", "PMEL", "TYR", "GAPDH", "BAP1", "GZMB"])
            @test all(==(80), cmp.cohort_n)
            # identical to R / Clojure
            @test only(cmp.z[cmp.hgnc_symbol .== "PMEL"]) ≈ -9.743 atol = 1e-3
            ch = compare_samples("H37001-003", "H37001-001"; db=h)
            @test all(>=(0.5), ch.avg_log10)
            @test ch.hgnc_symbol[1] == "TYRP1"
            @test plot_ma(ch) isa PatternQ.PlotlyBase.Plot
        end
    end
else
    @info "PATTERNQ_API_KEY not set: skipping live tests"
end
