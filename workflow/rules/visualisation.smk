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


rule visualise_trimmed:
    input:
        "results/{date}/out/trimmed-seqs.qza",
    output:
        "results/{date}/visual/trimmed-seqs.qzv",
    log:
        "logs/{date}/visualisation/visualise-trimmed.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output} "
        "--verbose 2> {log}"


rule visualise_joined:
    input:
        "results/{date}/out/joined-seqs.qza",
    output:
        "results/{date}/visual/joined-seqs.qzv",
    log:
        "logs/{date}/visualisation/visualise-joined.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime demux summarize "
        "--i-data {input} "
        "--o-visualization {output} "
        "--verbose 2> {log}"


rule unzip_samples:
    input:
        "results/{date}/visual/paired-seqs.qzv",
    output:
        temp(directory("results/{date}/visual/paired-seqs")),
    log:
        "logs/{date}/outputs/unzip-samples.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


rule unzip_trimmed:
    input:
        "results/{date}/visual/trimmed-seqs.qzv",
    output:
        temp(directory("results/{date}/visual/trimmed-seqs")),
    log:
        "logs/{date}/outputs/unzip-trimmed.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


rule unzip_joined:
    input:
        "results/{date}/visual/joined-seqs.qzv",
    output:
        temp(directory("results/{date}/visual/joined-seqs")),
    log:
        "logs/{date}/outputs/unzip-joined.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


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


rule unzip_frequency_length:
    input:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
    output:
        temp(directory("results/{date}/visual/lengthfilter_unzip")),
    log:
        "logs/{date}/outputs/unzip-length.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


rule visualise_afterab:
    input:
        "results/{date}/out/table-cluster-filtered.qza",
    output:
        "results/{date}/visual/table-cluster-filtered.qzv",
    log:
        "logs/{date}/visualisation/visualise-table.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-table summarize "
        "--i-table {input} "
        "--o-visualization {output} "
        "--verbose 2> {log} "


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
        "../envs/python.yaml"
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
        "../envs/python.yaml"
    shell:
        """
        (biom add-metadata -i {input.biom}/feature-table.biom -o {output.biom} --observation-metadata-fp {input.taxonomy}) > {log} 2>&1

        (biom convert -i {output.biom} -o {output.txt} --to-tsv --header-key taxonomy) > {log} 2>&1
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
        "../envs/python.yaml"
    shell:
        """
        (biom add-metadata -i {input.biom}/feature-table.biom -o {output.biom} --observation-metadata-fp {input.taxonomy}) > {log} 2>&1

        (biom convert -i {output.biom} -o {output.txt} --to-tsv --header-key taxonomy) > {log} 2>&1
        """


rule binary_heatmap:
    input:
        "results/{date}/out/table.from_biom_w_taxonomy.txt",
    output:
        report(
            "results/{date}/visual/heatmap_binary.html",
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
            "results/{date}/visual/absolute-taxabar-plot.html",
            caption="../report/absolute-taxabar-plot.rst",
            category="2. Taxonomy",
            subcategory="Taxa Barplot",
        ),
    params:
        samplename=config["metadata-parameters"]["absolute-taxa-name"],
        metadata="config/pep/sample.tsv",
    log:
        "logs/{date}/visualisation/absolute_taxabarplot.log",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/absolute_taxabarplot.py"


rule copy_diversity:
    input:
        expand(
            "results/{{date}}/out/beta-diversity-{metric}-normal.qza",
            metric=get_metric("beta"),
        ),
        expand(
            "results/{{date}}/out/beta-diversity-{metric}-phylogenetic.qza",
            metric=get_phylogenetic_metric("beta"),
        ),
    output:
        directory("results/{date}/out/distance_matrices/"),
    log:
        "logs/{date}/visualisation/distance_matrices_copy.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prepare_diversity.py"


rule create_heatmap:
    input:
        "results/{date}/out/distance_matrices/",
    output:
        report(
            "results/{date}/visual/beta-diversity-{metric}.html",
            caption="../report/distance-matrices.rst",
            category="3. Analysis",
            subcategory="Beta",
        ),
    log:
        "logs/{date}/visualisation/beta-diversity-{metric}.log",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/plot_distance.py"


rule rank_abundance:
    input:
        "results/{date}/out/taxa_collapsed_relative.qza",
    output:
        folder=report(
            directory("results/{date}/visual/report/rank-abundance/plots/"),
            caption="../report/rank-abundance.rst",
            htmlindex="rank-abundance.html",
            category="3. Analysis",
            subcategory="Rank abundance",
        ),
        file="results/{date}/visual/report/rank-abundance/plots/rank-abundance.html",
    params:
        "results/{date}/visual/report/rank-abundance/",
    log:
        "logs/{date}/visualisation/rank-abundance.log",
    conda:
        "../envs/plot.yaml"
    script:
        "../scripts/rank-abundance.py"


if config["DADA2"] == False:

    rule all_filter:
        input:
            samples="results/{date}/visual/paired-seqs",
            trimmed="results/{date}/visual/trimmed-seqs",
            joined="results/{date}/visual/joined-seqs/",
            first="results/{date}/visual/report/demux-joined-filter-stats/",
            length="results/{date}/visual/lengthfilter_unzip/",
            before_abundance="results/{date}/visual/table-cluster-lengthfilter/data/",
            final="results/{date}/visual/report/table-cluster-filtered/",
        output:
            report(
                "results/{date}/visual/allfilter.html",
                caption="../report/all-filter.rst",
                category="4. Qualitycontrol",
            ),
        log:
            "logs/{date}/visualisation/all-filter.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/complete_filter.py"


if config["DADA2"] == True:

    rule all_filter:
        input:
            dada2="results/{date}/visual/unzipped/",
            length="results/{date}/visual/lengthfilter_unzip/",
            before_abundance="results/{date}/visual/table-cluster-lengthfilter/data/",
            final="results/{date}/visual/report/table-cluster-filtered/",
        output:
            report(
                "results/{date}/visual/allfilter.html",
                caption="../report/all-filter.rst",
                category="4. Qualitycontrol",
            ),
        log:
            "logs/{date}/visualisation/all-filter.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/complete_filter_DADA2.py"


rule empress_tree:
    input:
        tree="results/{date}/visual/rooted-tree.qza",
        table="results/{date}/out/table-cluster-filtered.qza",
    output:
        "results/{date}/visual/empress-community.qzv",
    log:
        "logs/{date}/visualisation/empress-treeviewer.log",
    params:
        metadata="config/pep/sample.tsv",
        taxonomy="results/{date}/out/taxonomy.qza",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime empress community-plot "
        "--i-tree {input.tree} "
        "--i-feature-table {input.table} "
        "--m-sample-metadata-file {params.metadata} "
        "--m-feature-metadata-file {params.taxonomy} "
        "--p-filter-extra-samples "
        "--p-filter-missing-features "
        "--o-visualization {output} "


rule include_metadata:
    input:
        "config/pep/sample.tsv",
    output:
        report(
            "results/{date}/visual/report/sample.tsv",
            caption="../report/metadata.rst",
            category="4. Qualitycontrol",
        ),
    log:
        "logs/{date}/visualisation/metadata.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cp {input} {output}"
