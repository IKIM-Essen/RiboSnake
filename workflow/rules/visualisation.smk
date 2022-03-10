rule visualise_samples:
    input:
        "results/{date}/out/demux-paired-end.qza"
    output:
        "results/{date}/visual/paired-seqs.qzv"
    params:
        date = get_date()
    shell:
        "qiime demux summarize "
            "--i-data {input} "
            "--o-visualization {output}"

rule visualise_table:
    input:
        "results/{date}/out/table-cluster-lengthfilter.qza"
    output:
        "results/{date}/visual/table-cluster-lengthfilter.qzv"
    shell:
        "qiime feature-table summarize " 
            "--i-table {input} "
            "--o-visualization {output}"

rule visualise_fastq:
    input:
        "results/{date}/out/demux-paired-end.qza"
    output:
        "results/{date}/visual/fastq_stats.qzv"
    shell:
        "qiime vsearch fastq-stats "
            "--i-sequences {input} "
            "--o-visualization {output}"

rule demux_stats:
    input:
        "results/{date}/out/demux-joined-filter-stats.qza"
    output:
        "results/{date}/visual/demux-joined-filter-stats.qzv"
    shell:
        "qiime metadata tabulate " 
            "--m-input-file {input} "
            "--o-visualization {output} "
        

rule taxa_heatmap:
    input:
        "results/{date}/out/taxa_collapsed.qza"
    output:
        "results/{date}/visual/heatmap.qzv"
    params:
        metadata = "swab-site",
        cluster = "features"
    shell:
        "qiime feature-table heatmap "
            "--i-table {input} "
            "--m-sample-metadata-file config/pep/sample.tsv "
            "--m-sample-metadata-column {params.metadata} "
            "--p-cluster {params.cluster} "
            "--o-visualization {output}"

rule taxa_barplot:
    input:
        table = "results/{date}/out/table-taxa-filtered.qza",#"results/{date}/out/table-cluster-filtered.qza",
        taxonomy = "results/{date}/out/taxonomy.qza"
    output:
        "results/{date}/visual/taxa-bar-plots.qzv"
    shell:
        "qiime taxa barplot "
            "--i-table {input.table} "
            "--i-taxonomy {input.taxonomy} "
            "--m-metadata-file config/pep/sample.tsv "
            "--o-visualization {output}"

rule tabulate_taxa:
    input:
        "results/{date}/out/taxa_collapsed.qza"
    output:
        "results/{date}/visual/taxonomy.qzv"
    shell:
        "qiime metadata tabulate "
            "--m-input-file {input} "
            "--o-visualization {output}"

rule newick_tree:
    input:
        "results/{date}/visual/rooted-tree.qza"
    output:
        "results/{date}/visual"
    shell:
        "qiime tools export "
            "--input-path {input} "
            "--output-path {output}"

rule alpha:
    input:
        faith = "results/{date}/out/core-metrics-results/faith_pd_vector.qza",
        evenness = "results/{date}/out/core-metrics-results/evenness_vector.qza"
    output:
        faith = "results/{date}/visual/faith-pd-group-significance.qzv",
        evenness = "results/{date}/visual/evenness-group-significance.qzv"
    shell:
        "qiime diversity alpha-group-significance "
            "--i-alpha-diversity {input.faith} "
            "--m-metadata-file config/pep/sample.tsv "
            "--o-visualization {output.faith} \n"
        "qiime diversity alpha-group-significance "
            "--i-alpha-diversity {input.evenness} "
            "--m-metadata-file config/pep/sample.tsv "
            "--o-visualization {output.evenness}"

rule beta:
    input:
        "results/{date}/out/core-metrics-results/unweighted_unifrac_distance_matrix.qza"
    output:
        "results/{date}/visual/unweighted-unifrac-body-site-significance.qzv"
    params:
        metadata = "swab-site"
    shell:
        "qiime diversity beta-group-significance "
            "--i-distance-matrix {input} "
            "--m-metadata-file config/pep/sample.tsv "
            "--m-metadata-column {params.metadata} "
            "--o-visualization {output} "
            "--p-pairwise"

rule emperor:
    input:
        unifrac_pcoa = "results/{date}/out/core-metrics-results/unweighted_unifrac_pcoa_results.qza",
        bray_curtis_pcoa =  "results/{date}/out/core-metrics-results/bray_curtis_pcoa_results.qza"
    output:
        unifrac = "results/{date}/visual/unweighted-unifrac-emperor-days-since-experiment-start.qzv",
        bray_curtis = "results/{date}/visual/bray-curtis-emperor-days-since-experiment-start.qzv"
    shell:
        "qiime emperor plot "
            "--i-pcoa {input.unifrac_pcoa} "
            "--m-metadata-file config/pep/sample.tsv "
            "--p-custom-axes year "
            "--o-visualization {output.unifrac} \n"

        "qiime emperor plot "
            "--i-pcoa {input.bray_curtis_pcoa} "
            "--m-metadata-file config/pep/sample.tsv "
            "--p-custom-axes year "
            "--o-visualization {output.bray_curtis}"

rule rarefaction:
    input:
        table = "results/{date}/out/table-taxa-filtered.qza",
        phylogeny = "results/{date}/visual/rooted-tree.qza"
    output:
        alpha = "results/{date}/visual/alpha-rarefaction.qzv",
        beta = "results/{date}/visual/beta-rarefaction.qzv"
    shell:
        "qiime diversity alpha-rarefaction "
            "--i-table {input.table} "
            "--i-phylogeny {input.phylogeny} "
            "--p-max-depth 500 "
            "--m-metadata-file config/pep/sample.tsv "
            "--o-visualization {output.alpha} \n"
        "qiime diversity beta-rarefaction "
            "--i-table {input.table} "
            "--i-phylogeny {input.phylogeny} "
            "--p-metric euclidean "
            "--p-clustering-method nj "
            "--m-metadata-file config/pep/sample.tsv "
            "--p-sampling-depth 100 "
            "--o-visualization {output.beta}"

rule gneiss:
    input:
        table = "results/{date}/out/taxa_collapsed.qza"
    output:
        tree = "results/{date}/out/hirarchy_gneiss.qza",
        heatmap_gneiss = "results/{date}/visual/heatmap_gneiss.qzv"
    shell:
        "qiime gneiss correlation-clustering "
            "--i-table {input} "
            "--o-clustering {output.tree} \n"

        "qiime gneiss dendrogram-heatmap "
            "--i-table {input} "
            "--i-tree {output.tree} "
            "--m-metadata-file config/pep/sample.tsv "
            "--m-metadata-column subject "
            "--p-color-map seismic "
            "--o-visualization {output.heatmap_gneiss}"

rule rename_taxonomy:
    input:
        "results/{date}/out/taxonomy_biom/"
    output:
        "results/{date}/out/taxonomy_biom.tsv"
    script:
        "../scripts/rename_taxonomy.py"

rule add_biom_files:
    input:
        biom = "results/{date}/out/binary_biom/",
        taxonomy = "results/{date}/out/taxonomy_biom.tsv"
    output:
        biom = "results/{date}/out/table.w-taxa.biom",
        txt = "results/{date}/out/table.from_biom_w_taxonomy.txt"
    shell:
        """
            biom add-metadata -i {input.biom}/feature-table.biom -o {output.biom} --observation-metadata-fp {input.taxonomy}

            biom convert -i {output.biom} -o {output.txt} --to-tsv --header-key taxonomy
        """

rule binary_heatmap:
    input:
        "results/{date}/out/table.from_biom_w_taxonomy.txt"
    output:
        report(
            "results/{date}/visual/heatmap_binary.png",
            category = "1. Heatmap",
            subcategory = "Presence/absence heatmap",
        ),
    script:
        "../scripts/binaryheatmap.py"
