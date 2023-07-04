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
        "../envs/python.yaml"
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


rule biom_file:
    input:
        table="results/{date}/out/table-cluster-filtered.qza",
        taxonomy="results/{date}/out/taxonomy.qza",
    output:
        table_binary="results/{date}/out/table_binary.qza",
        table_biom=directory("results/{date}/out/biom_table/"),
        taxa_biom=directory("results/{date}/out/taxonomy_biom/"),
        binary_biom=directory("results/{date}/out/binary_biom/"),
    log:
        "logs/{date}/outputs/biom-file.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools export "
        "--input-path {input.table} "
        "--output-path {output.table_biom} \n"
        "qiime tools export "
        "--input-path {input.taxonomy} "
        "--output-path {output.taxa_biom} \n"
        "qiime feature-table presence-absence "
        "--i-table {input.table} "
        "--o-presence-absence-table {output.table_binary} \n"
        "qiime tools export "
        "--input-path {output.table_binary} "
        "--output-path {output.binary_biom}"


rule unzip_reports:
    input:
        "results/{date}/visual/heatmap.qzv",
        "results/{date}/visual/taxa-bar-plots.qzv",
        "results/{date}/visual/rooted-tree.qza",
        "results/{date}/visual/taxonomy.qzv",
        "results/{date}/visual/table-whuman.qzv",
        "results/{date}/visual/table-wohuman.qzv",
        "results/{date}/visual/paired-seqs.qzv",
        "results/{date}/visual/fastq_stats.qzv",
        "results/{date}/visual/demux-joined-filter-stats.qzv",
    output:
        directory("results/{date}/visual/unzipped"),
    log:
        "logs/{date}/outputs/unzip-reports.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


rule report_files:
    input:
        "results/{date}/visual/unzipped/",
    output:
        heatmap=report(
            "results/{date}/visual/report/heatmap.svg",
            caption="../report/heatmap.rst",
            category="1. Heatmap",
            subcategory="Relative abundances",
        ),
        taxonomy_tsv=report(
            "results/{date}/visual/report/taxonomy.tsv",
            caption="../report/taxonomy-tsv.rst",
            category="2. Taxonomy",
            subcategory="Taxonomy Table",
        ),
        taxa_barplot=report(
            directory("results/{date}/visual/report/taxa_barplot_data"),
            caption="../report/taxa-barplot.rst",
            category="2. Taxonomy",
            subcategory="Taxa Barplot",
            htmlindex="index.html",
        ),
        paired_seqs=report(
            directory("results/{date}/visual/report/paired_seqs"),
            caption="../report/paired-seqs.rst",
            category="4. Qualitycontrol",
            htmlindex="index.html",
        ),
        fastq_stats=report(
            directory("results/{date}/visual/report/fastq_stats"),
            caption="../report/fastq-stats.rst",
            category="4. Qualitycontrol",
            htmlindex="index.html",
        ),
        demux_filter_stats=report(
            directory("results/{date}/visual/report/demux-joined-filter-stats"),
            caption="../report/demux-filter-stats.rst",
            category="4. Qualitycontrol",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-files.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_reports.py"


rule parameter_summary:
    output:
        report(
            "results/{date}/out/parameter-summary.csv",
            caption="../report/parameter-summary.rst",
            category="4. Qualitycontrol",
        ),
    log:
        "logs/{date}/outputs/parameter_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/parameter_summary.py"


rule snakemake_report:
    input:
        "results/{date}/visual/heatmap_binary.png",
        "results/{date}/visual/report/heatmap.svg",
        "results/{date}/visual/unzipped",
        "results/{date}/visual/report/multiqc.html",
        "results/{date}/visual/absolute-taxabar-plot.png",
        "results/{date}/visual/report/human-count",
    output:
        "results/{date}/out/report.zip",
    params:
        for_testing=get_if_testing("--snakefile ../workflow/Snakefile"),
    log:
        "logs/{date}/outputs/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} --report-stylesheet resources/custom-stylesheet.css "
        "{params.for_testing} "
        "> {log} 2>&1"


rule compress_kraken:
    input:
        expand(
            "results/{{date}}/out/kraken/{sample}.kreport2",
            sample=get_reads_for_kraken(),
        ),
    output:
        "results/{date}/out/kraken.tar.gz",
    params:
        directory="results/{date}/out/kraken/",
    log:
        "logs/{date}/outputs/kraken-compress.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "tar -czvf {output} {params.directory} "


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


rule export_parameters:
    input:
        "config/config.yaml",
    output:
        report(
            "results/{date}/out/config_parameters.html",
            caption="../report/parameter-summary.rst",
            category="4. Qualitycontrol",
        ),
    log:
        "logs/{date}/outputs/config_html.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/yaml_to_table.py"


rule zip_report:
    input:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
        "results/{date}/visual/rooted-tree.qza",
        "results/{date}/out/taxonomy.qza",
        "results/{date}/out/aligned-rep-seqs.qza",
        "results/{date}/out/biom_table/",
        "results/{date}/out/taxonomy_biom/",
        "results/{date}/out/binary_biom/",
        "results/{date}/out/multiqc.html",
        "results/{date}/visual/heatmap_binary.png",
        "results/{date}/visual/report/heatmap.svg",
        "results/{date}/visual/report/taxonomy.tsv",
        "results/{date}/out/report.zip",
        "results/{date}/visual/fastq_stats.qzv",
        "results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
        "results/{date}/visual/absolute-taxabar-plot.png",
        "results/{date}/out/kraken.tar.gz",
        "results/{date}/out/parameter-summary.csv",
        "results/{date}/visual/sample_frequencys_difference.csv",
        "results/{date}/out/config_parameters.html",
    output:
        "results/{date}/16S-report.tar.gz",
    log:
        "logs/{date}/outputs/zip-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        """
        mkdir results/{wildcards.date}/16S-report
        cp -r {input} results/{wildcards.date}/16S-report/
        tar -czvf results/{wildcards.date}/16S-report.tar.gz results/{wildcards.date}/16S-report/
        """
