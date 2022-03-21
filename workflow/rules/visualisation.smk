rule visualise_samples:
    input:
        "results/{date}/out/demux-paired-end.qza",
    output:
        "results/{date}/visual/paired-seqs.qzv",
    log:
        "logs/{date}/visualisation/visualise-samples.log",
    conda:
        "../envs/qiime-metadata-visual.yaml"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output}"


rule visualise_table:
    input:
        "results/{date}/out/table-cluster-lengthfilter.qza",
    output:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
    log:
        "logs/{date}/visualisation/visualise-table.log",
    conda:
        "../envs/qiime-metadata-visual.yaml"
    shell:
        "qiime feature-table summarize "
        "--i-table {input} "
        "--o-visualization {output}"


rule visualise_fastq:
    input:
        "results/{date}/out/demux-paired-end.qza",
    output:
        "results/{date}/visual/fastq_stats.qzv",
    log:
        "logs/{date}/visualisation/visualise-fastq.log",
    conda:
        "envs/qiime-vsearch.yaml"
    shell:
        "qiime vsearch fastq-stats "
        "--i-sequences {input} "
        "--o-visualization {output}"


rule demux_stats:
    input:
        "results/{date}/out/demux-joined-filter-stats.qza",
    output:
        "results/{date}/visual/demux-joined-filter-stats.qzv",
    log:
        "logs/{date}/visualisation/demux-stats.log",
    conda:
        "../envs/qiime-metadata-visual.yaml"
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input} "
        "--o-visualization {output} "


rule taxa_heatmap:
    input:
        "results/{date}/out/taxa_collapsed.qza",
    output:
        "results/{date}/visual/heatmap.qzv",
    params:
        metadata=config["metadata-parameters"]["taxa-heatmap-column"],  #"extract-group-no",#"swab-site", config["metadata-parameters"]["taxa-heatmap-column"]
        cluster="features",
    log:
        "logs/{date}/visualisation/taxa-heatmap.log",
    conda:
        "../envs/qiime-taxonomy.yaml"
    shell:
        "qiime feature-table heatmap "
        "--i-table {input} "
        "--m-sample-metadata-file config/pep/sample.tsv "
        "--m-sample-metadata-column {params.metadata} "
        "--p-cluster {params.cluster} "
        "--o-visualization {output}"


rule taxa_barplot:
    input:
        table="results/{date}/out/table-taxa-filtered.qza",  #"results/{date}/out/table-cluster-filtered.qza",
        taxonomy="results/{date}/out/taxonomy.qza",
    output:
        "results/{date}/visual/taxa-bar-plots.qzv",
    log:
        "logs/{date}/visualisation/taxa-barplot.log",
    conda:
        "../envs/qiime-taxonomy.yaml"
    shell:
        "qiime taxa barplot "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--m-metadata-file config/pep/sample.tsv "
        "--o-visualization {output}"


rule tabulate_taxa:
    input:
        "results/{date}/out/taxa_collapsed.qza",
    output:
        "results/{date}/visual/taxonomy.qzv",
    log:
        "logs/{date}/visualisation/tabulate-taxa.log",
    conda:
        "../envs/qiime-metadata-visual.yaml"
    shell:
        "qiime metadata tabulate "
        "--m-input-file {input} "
        "--o-visualization {output}"


rule newick_tree:
    input:
        "results/{date}/visual/rooted-tree.qza",
    output:
        "results/{date}/visual",
    log:
        "logs/{date}/visualisation/newick-tree.log",
    conda:
        "../envs/qiime-qiime-exports.yaml"
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output}"


rule alpha:
    input:
        faith="results/{date}/out/core-metrics-results/faith_pd_vector.qza",
        evenness="results/{date}/out/core-metrics-results/evenness_vector.qza",
    output:
        faith="results/{date}/visual/faith-pd-group-significance.qzv",
        evenness="results/{date}/visual/evenness-group-significance.qzv",
    log:
        "logs/{date}/visualisation/alpha-rarefaction.log",
    conda:
        "../envs/qiime-diversity.yaml"
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
        "results/{date}/out/core-metrics-results/unweighted_unifrac_distance_matrix.qza",
    output:
        "results/{date}/visual/unweighted-unifrac-body-site-significance.qzv",
    params:
        metadata=config["metadata-parameters"]["beta-metadata-column"],  #"extract-group-no"#"swab-site", config["metadata-parameters"]["beta-metadata-column"]
    log:
        "logs/{date}/visualisation/beta-rarefaction.log",
    conda:
        "../envs/qiime-diversity.yaml"
    shell:
        "qiime diversity beta-group-significance "
        "--i-distance-matrix {input} "
        "--m-metadata-file config/pep/sample.tsv "
        "--m-metadata-column {params.metadata} "
        "--o-visualization {output} "
        "--p-pairwise"


rule emperor:
    input:
        unifrac_pcoa="results/{date}/out/core-metrics-results/unweighted_unifrac_pcoa_results.qza",
        bray_curtis_pcoa=(
            "results/{date}/out/core-metrics-results/bray_curtis_pcoa_results.qza"
        ),
    output:
        unifrac="results/{date}/visual/unweighted-unifrac-emperor-days-since-experiment-start.qzv",
        bray_curtis=(
            "results/{date}/visual/bray-curtis-emperor-days-since-experiment-start.qzv"
        ),
    log:
        "logs/{date}/visualisation/emperor.log",
    conda:
        "../envs/qiime-plots.yaml"
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
        table="results/{date}/out/table-taxa-filtered.qza",
        phylogeny="results/{date}/visual/rooted-tree.qza",
    output:
        alpha="results/{date}/visual/alpha-rarefaction.qzv",
        beta="results/{date}/visual/beta-rarefaction.qzv",
    params:
        max_depth=config["metadata-parameters"]["rarefaction-max-depth"],
    log:
        "logs/{date}/visualisation/rarefaction.log",
    conda:
        "../envs/qiime-diversity.yaml"
    shell:
        "qiime diversity alpha-rarefaction "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-max-depth {params.max_depth} "
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
        #400,500


rule gneiss:
    input:
        table="results/{date}/out/taxa_collapsed.qza",
    output:
        tree="results/{date}/out/hirarchy_gneiss.qza",
        heatmap_gneiss="results/{date}/visual/heatmap_gneiss.qzv",
    params:
        metadata=config["metadata-parameters"]["gneiss-metadata-column"],
    log:
        "logs/{date}/visualisation/gneiss.log",
    conda:
        "../envs/qiime-plots.yaml"
    shell:
        "qiime gneiss correlation-clustering "
        "--i-table {input} "
        "--o-clustering {output.tree} \n"

        "qiime gneiss dendrogram-heatmap "
        "--i-table {input} "
        "--i-tree {output.tree} "
        "--m-metadata-file config/pep/sample.tsv "
        "--m-metadata-column {params.metadata} "
        "--p-color-map seismic "
        "--o-visualization {output.heatmap_gneiss}"
        #subject


rule rename_taxonomy:
    input:
        "results/{date}/out/taxonomy_biom/",
    output:
        "results/{date}/out/taxonomy_biom.tsv",
    log:
        "logs/{date}/visualisation/rename-taxonomy.log",
    conda:
        "../envs/qiime-taxonomy.yaml"
    script:
        "../scripts/rename_taxonomy.py"


rule add_biom_files:
    input:
        biom="results/{date}/out/binary_biom/",
        taxonomy="results/{date}/out/taxonomy_biom.tsv",
    output:
        biom="results/{date}/out/table.w-taxa.biom",
        txt="results/{date}/out/table.from_biom_w_taxonomy.txt",
    log:
        "logs/{date}/visualisation/add-biom-files.log",
    conda:
        "../envs/biom.yaml"
    shell:
        """
        biom add-metadata -i {input.biom}/feature-table.biom -o {output.biom} --observation-metadata-fp {input.taxonomy}

        biom convert -i {output.biom} -o {output.txt} --to-tsv --header-key taxonomy
        """


rule binary_heatmap:
    input:
        "results/{date}/out/table.from_biom_w_taxonomy.txt",
    output:
        report(
            "results/{date}/visual/heatmap_binary.png",
            category="1. Heatmap",
            subcategory="Presence/absence heatmap",
        ),
    log:
        "logs/{date}/visualisation/binary-heatmap.log",
    conda:
        "../envs/binary-heatmap.yaml"
    script:
        "../scripts/binaryheatmap.py"
