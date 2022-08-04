rule visualise_samples:
    input:
        "results/{date}/out/demux-paired-end.qza",
    output:
        "results/{date}/visual/paired-seqs.qzv",
    log:
        "logs/{date}/visualisation/visualise-samples.log",
    conda:
        "../envs/qiime-only-env.yaml"
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
        "../envs/qiime-only-env.yaml"
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
        "../envs/qiime-only-env.yaml"
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
        "../envs/qiime-only-env.yaml"
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
        cluster=config["metadata-parameters"]["cluster"],
    log:
        "logs/{date}/visualisation/taxa-heatmap.log",
    conda:
        "../envs/qiime-only-env.yaml"
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
        "../envs/qiime-only-env.yaml"
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
        "../envs/qiime-only-env.yaml"
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
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output}"


rule alpha_significance:
    input:
        direc="results/{date}/core-metrics-results",
    output:
        faith="results/{date}/visual/faith-pd-group-significance.qzv",
        evenness="results/{date}/visual/evenness-group-significance.qzv",
    params:
        faith="results/{date}/core-metrics-results/faith_pd_vector.qza",
        evenness="results/{date}/core-metrics-results/evenness_vector.qza",
    log:
        "logs/{date}/visualisation/alpha-significance.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha-group-significance "
        "--i-alpha-diversity {params.faith} "
        "--m-metadata-file config/pep/sample.tsv "
        "--o-visualization {output.faith} \n"
        "qiime diversity alpha-group-significance "
        "--i-alpha-diversity {params.evenness} "
        "--m-metadata-file config/pep/sample.tsv "
        "--o-visualization {output.evenness}"


rule beta_significance:
    input:
        direc="results/{date}/core-metrics-results",
    output:
        "results/{date}/visual/unweighted-unifrac-body-site-significance.qzv",
    params:
        metadata=config["metadata-parameters"]["beta-metadata-column"],  #"extract-group-no"#"swab-site", config["metadata-parameters"]["beta-metadata-column"]
        unifrac=(
            "results/{date}/core-metrics-results/unweighted_unifrac_distance_matrix.qza"
        ),
    log:
        "logs/{date}/visualisation/beta-significance.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta-group-significance "
        "--i-distance-matrix {params.unifrac} "
        "--m-metadata-file config/pep/sample.tsv "
        "--m-metadata-column {params.metadata} "
        "--o-visualization {output} "
        "--p-pairwise"


rule copy_emperor:
    input:
        direc="results/{date}/core-metrics-results",
        scatter="results/{date}/out/beta-correlation-scatter.qzv",
    output:
        bray="results/{date}/visual/bray-curtis-emperor.qzv",
        jaccard="results/{date}/visual/jaccard-emperor.qzv",
        unweighted_unifrac="results/{date}/visual/unweighted-unifrac-emperor.qzv",
        weighted_unifrac="results/{date}/visual/weighted-unifrac-emperor.qzv",
        scatter="results/{date}/visual/beta-correlation-scatter.qzv",
    params:
        bray="results/{date}/core-metrics-results/bray_curtis_emperor.qzv",
        jaccard="results/{date}/core-metrics-results/jaccard_emperor.qzv",
        unweighted_unifrac="results/{date}/core-metrics-results/unweighted_unifrac_emperor.qzv",
        weighted_unifrac="results/{date}/core-metrics-results/weighted_unifrac_emperor.qzv",
    log:
        "logs/{date}/visualisation/copy-emperor.log",
    conda:
        "../envs/python.yaml",
    shell:
        "cp {params.bray} {output.bray};"
        "cp {params.jaccard} {output.jaccard};"
        "cp {params.unweighted_unifrac} {output.unweighted_unifrac};"
        "cp {params.weighted_unifrac} {output.weighted_unifrac};"
        "cp {input.scatter} {output.scatter};"


rule emperor:
    input:
        direc="results/{date}/core-metrics-results",
    output:
        unifrac="results/{date}/visual/unweighted-unifrac-emperor-days-since-experiment-start.qzv",
        bray_curtis=(
            "results/{date}/visual/bray-curtis-emperor-days-since-experiment-start.qzv"
        ),
    params:
        unifrac_pcoa=(
            "results/{date}/core-metrics-results/unweighted_unifrac_pcoa_results.qza"
        ),
        bray_curtis_pcoa=(
            "results/{date}/core-metrics-results/bray_curtis_pcoa_results.qza"
        ),
    log:
        "logs/{date}/visualisation/emperor.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime emperor plot "
        "--i-pcoa {params.unifrac_pcoa} "
        "--m-metadata-file config/pep/sample.tsv "
        "--p-custom-axes year "
        "--o-visualization {output.unifrac} \n"

        "qiime emperor plot "
        "--i-pcoa {params.bray_curtis_pcoa} "
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
        max_depth=config["rarefaction"]["max-depth"],
        sampling_depth=config["rarefaction"]["sampling_depth"],
        metric=config["rarefaction"]["metric"],
        clustering_method=config["rarefaction"]["clustering_method"],
    log:
        "logs/{date}/visualisation/rarefaction.log",
    conda:
        "../envs/qiime-only-env.yaml"
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
        "--p-metric {params.metric} "
        "--p-clustering-method {params.clustering_method} "
        "--m-metadata-file config/pep/sample.tsv "
        "--p-sampling-depth {params.sampling_depth} "
        "--o-visualization {output.beta}"


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
        "../envs/qiime-only-env.yaml"
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


rule add_biom_files_featcount:
    input:
        biom="results/{date}/out/biom_table/",
        taxonomy="results/{date}/out/taxonomy_biom.tsv",
    output:
        biom="results/{date}/out/table.w-taxa-featcount.biom",
        txt="results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
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
            caption="../report/binary-heatmap.rst",
            category="1. Heatmap",
            subcategory="Presence/absence heatmap",
        ),
    log:
        "logs/{date}/visualisation/binary-heatmap.log",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/binaryheatmap.py"


rule absolute_taxa:
    input:
        "results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
    output:
        report(
            "results/{date}/visual/absolute-taxabar-plot.png",
            caption="../report/absolute-taxabar-plot.rst",
            category="2. Taxonomy",
            subcategory="Taxa Barplot",
        ),
    log:
        "logs/{date}/visualisation/absolute_taxabarplot.log",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/absolute_taxabarplot.py"


rule alpha:
    input:
        "results/{date}/out/table-taxa-filtered.qza",
    output:
        "results/{date}/out/alpha-diversity.qza",
    params:
        metric=config["diversity"]["alpha"]["diversity-metric"],
    log:
        "logs/{date}/visualisation/alpha-diversity.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha "
        "--i-table {input} "
        "--p-metric {params.metric} "
        "--o-alpha-diversity {output} "
        "--verbose 2> {log}"
    

rule beta:
    input:
        "results/{date}/out/table-taxa-filtered.qza",
    output:
        "results/{date}/out/beta-diversity-distance.qza",
    params:
        metric=config["diversity"]["beta"]["diversity-metric"],
        pseudocount=config["diversity"]["beta"]["diversity-pseudocount"],
        n_jobs=config["diversity"]["beta"]["diversity-n-jobs"],
    log:
        "logs/{date}/visualisation/beta-diversity.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta "
        "--i-table {input} "
        "--p-metric {params.metric} "
        "--p-pseudocount {params.pseudocount} "
        "--p-n-jobs {params.n_jobs} "
        "--o-distance-matrix {output} "
        "--verbose 2> {log}"


rule alpha_phylogeny:
    input:
        phylogeny="results/{date}/visual/rooted-tree.qza",
        table="results/{date}/out/table-taxa-filtered.qza",
    output:
        "results/{date}/out/alpha-phylogeny.qza",
    params:
        metric=config["diversity"]["alpha"]["phylogeny-metric"],
    log:
        "logs/{date}/visualisation/alpha-phylogeny.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha-phylogenetic "
        "--i-phylogeny {input.phylogeny} "
        "--i-table {input.table} "
        "--p-metric {params.metric} "
        "--o-alpha-diversity {output} "
        "--verbose 2> {log}"


rule beta_phylogeny:
    input:
        table="results/{date}/out/table-taxa-filtered.qza",
        phylogeny="results/{date}/visual/rooted-tree.qza",
    output:
        "results/{date}/out/beta-phylogeny.qza",
    params:
        metrics=config["diversity"]["beta"]["phylogeny-metric"],
        threads=config["threads"],
        variance_adjusted=config["diversity"]["beta"]["phylogeny-variance-adjusted"],
    log:
        "logs/{date}/visualisation/beta-phylogeny.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta-phylogenetic "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-metric {params.metrics} "
        "--p-threads {params.threads} "
        "--p-variance-adjusted {params.variance_adjusted} "
        "--o-distance-matrix {output}"


rule beta_correlation:
    input:
        "results/{date}/core-metrics-results/",
    output:
        distance_matrix="results/{date}/out/beta-correlation.qza",
        mantel_scatter_vis="results/{date}/out/beta-correlation-scatter.qzv",
    params:
        matrix="results/{date}/core-metrics-results/jaccard_distance_matrix.qza",
        metadata_file="config/pep/sample.tsv",
        metadata_column=config["diversity"]["beta"]["correlation-column"],
        method=config["diversity"]["beta"]["correlation-method"],
        permutations=config["diversity"]["beta"]["correlation-permutations"],
    log:
        "logs/{date}/visualisation/beta-correlation.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta-correlation "
        "--i-distance-matrix {params.matrix} "
        "--m-metadata-file {params.metadata_file} "
        "--m-metadata-column {params.metadata_column} "
        "--p-method {params.method} "
        "--p-permutations {params.permutations} "
        "--o-metadata-distance-matrix {output.distance_matrix} "
        "--o-mantel-scatter-visualization {output.mantel_scatter_vis} "
        "--verbose 2> {log}"
