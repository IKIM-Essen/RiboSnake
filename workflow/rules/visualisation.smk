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
        "--o-visualization {output} "
        "--verbose 2> {log}"


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
        "--o-visualization {output} "
        "--verbose 2> {log}"


rule table_compare_human:
    input:
        table_wh="results/{date}/out/derepl-table.qza",
        table_woh="results/{date}/out/derep-table-nonhum.qza",
    output:
        visual_wh="results/{date}/visual/table-whuman.qzv",
        visual_woh="results/{date}/visual/table-wohuman.qzv",
    log:
        "logs/{date}/visualisation/table-compare-human.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-table summarize "
        "--i-table {input.table_wh} "
        "--o-visualization {output.visual_wh} "
        "--verbose 2> {log} \n"
        "qiime feature-table summarize "
        "--i-table {input.table_woh} "
        "--o-visualization {output.visual_woh} "
        "--verbose 2> {log}"


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
        "--o-visualization {output} "
        "--verbose 2> {log}"


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
        "--verbose 2> {log}"


if config["data-type"] == "human":

    rule visual_humancount:
        input:
            "results/{date}/out/human.qza",
        output:
            "results/{date}/visual/human-count.qzv",
        log:
            "logs/{date}/visualisation/human-count.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime feature-table tabulate-seqs "
            "--i-data {input} "
            "--o-visualization {output} "
            "--verbose 2> {log}"

    rule unzip_human_count:
        input:
            "results/{date}/visual/human-count.qzv",
        output:
            human_count=report(
                directory("results/{date}/visual/report/human-count"),
                caption="../report/human-count.rst",
                category="4. Qualitycontrol",
                htmlindex="index.html",
            ),
        params:
            between="results/{date}/visual/report/human-count-unzipped",
        log:
            "logs/{date}/visualisation/human-count-unzip.log",
        conda:
            "../envs/qiime-only-env.yaml"
        script:
            "../scripts/extract_humancount.py"


if config["data-type"] == "environmental":

    rule unzip_human_dummy:
        output:
            directory("results/{date}/visual/report/human-count"),
        log:
            "logs/{date}/visualisation/human-count-dummy.log",
        conda:
            "../envs/snakemake.yaml"
        shell:
            "mkdir {output}"


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
        "--o-visualization {output} "
        "--verbose 2> {log}"


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
        "--o-visualization {output} "
        "--verbose 2> {log}"


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
        "--o-visualization {output} "
        "--verbose 2> {log}"


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
        "--output-path {output} "
        "--verbose 2> {log}"


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
        "--o-visualization {output.alpha} "
        "--verbose 2> {log} \n"
        "qiime diversity beta-rarefaction "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-metric {params.metric} "
        "--p-clustering-method {params.clustering_method} "
        "--m-metadata-file config/pep/sample.tsv "
        "--p-sampling-depth {params.sampling_depth} "
        "--o-visualization {output.beta} "
        "--verbose 2> {log}"


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
        "--o-visualization {output.heatmap_gneiss} "
        "--verbose 2> {log}"


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


rule ancom:
    input:
        "results/{date}/out/taxa_collapsed.qza",
    output:
        pseudocount_table="results/{date}/out/pseudocount_table-{metadata_column}.qza",
        ancom_output="results/{date}/visual/ancom-{metadata_column}.qzv",
    params:
        metadata_column="{metadata_column}",
        metadata_file="config/pep/sample.tsv",
    log:
        "logs/{date}/visualisation/ancom-{metadata_column}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime composition add-pseudocount "
        "--i-table {input} "
        "--o-composition-table {output.pseudocount_table} \n"
        "qiime composition ancom "
        "--i-table {output.pseudocount_table} "
        "--m-metadata-file {params.metadata_file} "
        "--m-metadata-column {params.metadata_column} "
        "--o-visualization {output.ancom_output} "
        "--verbose 2> {log}"


rule hum_filter_difference:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            "results/{date}/visual/sample_frequencys_difference.csv",
            caption="../report/hum_filter_difference.rst",
            category="4. Qualitycontrol",
        ),
    params:
        visual_wh="results/{date}/visual/unzipped/table-whuman/data/sample-frequency-detail.csv",
        visual_woh="results/{date}/visual/unzipped/table-wohuman/data/sample-frequency-detail.csv",
    log:
        "logs/{date}/visualisation/frequency_difference.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sample_freq_difference.py"
