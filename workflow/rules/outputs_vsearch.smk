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
        expand(
            "results/{{date}}/visual/alpha-correlation-{metric_alpha}-normal.qzv",
            metric_alpha=get_metric("alpha"),
        ),
        expand(
            "results/{{date}}/visual/alpha-significance-{metric_alpha}-normal.qzv",
            metric_alpha=get_metric("alpha"),
        ),
        expand(
            "results/{{date}}/visual/alpha-correlation-{metric_alpha}-phylogenetic.qzv",
            metric_alpha=get_phylogenetic_metric("alpha"),
        ),
        expand(
            "results/{{date}}/visual/alpha-significance-{metric_alpha}-phylogenetic.qzv",
            metric_alpha=get_phylogenetic_metric("alpha"),
        ),
        expand(
            "results/{{date}}/visual/beta-significance-{metric}-normal-{metadata_column}.qzv",
            metric=get_metric("beta"),
            metadata_column=get_metadata_categorical_columns(),
        ),
        expand(
            "results/{{date}}/visual/beta-significance-{metric}-phylogenetic-{metadata_column}.qzv",
            metric=get_phylogenetic_metric("beta"),
            metadata_column=get_metadata_categorical_columns(),
        ),
        expand(
            "results/{{date}}/visual/beta-correlation-scatter-{metric}-normal-{metadata_column}.qzv",
            metric=get_metric("beta"),
            metadata_column=get_metadata_columns(),
        ),
        expand(
            "results/{{date}}/visual/beta-correlation-scatter-{metric}-phylogenetic-{metadata_column}.qzv",
            metric=get_phylogenetic_metric("beta"),
            metadata_column=get_metadata_columns(),
        ),
        expand(
            "results/{{date}}/visual/emperor-{metric}-{diversity}.qzv",
            diversity="normal",
            metric=get_metric("beta"),
        ),
        expand(
            "results/{{date}}/visual/emperor-{metric}-{diversity}.qzv",
            metric=get_phylogenetic_metric("beta"),
            diversity="phylogenetic",
        ),
        expand(
            "results/{{date}}/visual/ancom-{metadata_column}.qzv",
            metadata_column=get_ancom_columns(),
        ),
        "results/{date}/visual/alpha-rarefaction.qzv",
        "results/{date}/visual/beta-rarefaction.qzv",
        "results/{date}/visual/heatmap.qzv",
        "results/{date}/visual/taxa-bar-plots.qzv",
        "results/{date}/visual/rooted-tree.qza",
        "results/{date}/visual/taxonomy.qzv",
        "results/{date}/visual/table-whuman.qzv",
        "results/{date}/visual/table-wohuman.qzv",
        "results/{date}/visual/paired-seqs.qzv",
        "results/{date}/visual/fastq_stats.qzv",
        "results/{date}/visual/demux-joined-filter-stats.qzv",
        "results/{date}/visual/heatmap_gneiss.qzv",
        "results/{date}/visual/empress-community.qzv",
    output:
        temp(directory("results/{date}/visual/unzipped")),
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
        beta_svg=report(
            "results/{date}/visual/report/beta-rarefaction.svg",
            caption="../report/beta-heatmap.rst",
            category="3. Analysis",
            subcategory="Beta",
        ),
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
        beta_html=report(
            directory("results/{date}/visual/report/beta_rarefaction"),
            caption="../report/beta-rarefaction.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
        alpha_html=report(
            directory("results/{date}/visual/report/alpha_rarefaction"),
            caption="../report/alpha-rarefaction.rst",
            category="3. Analysis",
            subcategory="Alpha",
            htmlindex="index.html",
        ),
        gneiss=report(
            "results/{date}/visual/heatmap_gneiss.svg",
            caption="../report/gneiss.rst",
            category="3. Analysis",
            subcategory="Gneiss",
        ),
        paired_seqs=report(
            directory("results/{date}/visual/report/paired-seqs"),
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


rule report_beta_correlation:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory(
                "results/{date}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}"
            ),
            caption="../report/beta-correlation-scatter.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-beta-correlation-{metric}-{diversity}-{metadata_column}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_beta_corr.py"


rule report_beta_significance:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory(
                "results/{date}/visual/report/beta-significance-{metric}-{diversity}-{metadata_column}"
            ),
            caption="../report/beta-significance.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-significance-{metric}-{diversity}-{metadata_column}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_significance.py"


rule report_alpha_correlation:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory(
                "results/{date}/visual/report/alpha-correlation-{metric_alpha}-{diversity}"
            ),
            caption="../report/alpha-correlation.rst",
            category="3. Analysis",
            subcategory="Alpha",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-alpha-correlation-{metric_alpha}-{diversity}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_significance.py"


rule report_alpha_significance:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory(
                "results/{date}/visual/report/alpha-significance-{metric_alpha}-{diversity}"
            ),
            caption="../report/alpha-significance.rst",
            category="3. Analysis",
            subcategory="Alpha",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-alpha-significance-{metric_alpha}-{diversity}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_significance.py"


rule report_emperor:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory("results/{date}/visual/report/emperor-{metric}-{diversity}"),
            caption="../report/jaccard-emperor.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-emperor-{metric}-{diversity}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_significance.py"


rule report_empress:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory("results/{date}/visual/report/empress-community"),
            caption="../report/empress.rst",
            category="2. Taxonomy",
            subcategory="Phylogenetic Tree",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-empress.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_significance.py"


rule report_ancom:
    input:
        "results/{date}/visual/unzipped/",
    output:
        report(
            directory("results/{date}/visual/report/ancom-{metadata_column}"),
            caption="../report/ancom.rst",
            category="3. Analysis",
            subcategory="Ancom",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/ancom-{metadata_column}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_beta_corr.py"


if config["longitudinal"] == False:

    rule snakemake_report:
        input:
            "results/{date}/visual/heatmap_binary.html",
            "results/{date}/visual/report/beta-rarefaction.svg",
            "results/{date}/visual/report/heatmap.svg",
            "results/{date}/visual/unzipped",
            "results/{date}/visual/report/multiqc.html",
            "results/{date}/visual/report/empress-community",
            "results/{date}/visual/absolute-taxabar-plot.html",
            "results/{date}/out/qurro_plot",
            "results/{date}/visual/report/rank-abundance/plots/",
            "results/{date}/visual/allfilter.html",
            "results/{date}/visual/report/sample.tsv",
            expand(
                "results/{{date}}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}",
                metric=get_metric("beta"),
                metadata_column=get_metadata_columns(),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}",
                metric=get_phylogenetic_metric("beta"),
                metadata_column=get_metadata_columns(),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/beta-significance-{metric}-{diversity}-{metadata_column}",
                metric=get_phylogenetic_metric("beta"),
                metadata_column=get_metadata_categorical_columns(),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/beta-significance-{metric}-{diversity}-{metadata_column}",
                metric=get_metric("beta"),
                metadata_column=get_metadata_categorical_columns(),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-correlation-{metric_alpha}-{diversity}",
                metric_alpha=get_metric("alpha"),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-correlation-{metric_alpha}-{diversity}",
                metric_alpha=get_phylogenetic_metric("alpha"),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-significance-{metric_alpha}-{diversity}",
                metric_alpha=get_metric("alpha"),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-significance-{metric_alpha}-{diversity}",
                metric_alpha=get_phylogenetic_metric("alpha"),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/emperor-{metric}-{diversity}",
                metric=get_phylogenetic_metric("beta"),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/emperor-{metric}-{diversity}",
                metric=get_metric("beta"),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/ancom-{metadata_column}",
                metadata_column=get_ancom_columns(),
            ),
            expand(
                "results/{{date}}/visual/beta-diversity-{metric}.html",
                metric=get_complete_beta_metric(),
            ),
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


if config["longitudinal"]:

    rule snakemake_report:
        input:
            "results/{date}/visual/heatmap_binary.html",
            "results/{date}/visual/report/beta-rarefaction.svg",
            "results/{date}/visual/report/heatmap.svg",
            "results/{date}/visual/unzipped",
            "results/{date}/visual/longitudinal_unzipped",
            "results/{date}/visual/report/multiqc.html",
            "results/{date}/visual/absolute-taxabar-plot.html",
            "results/{date}/out/qurro_plot",
            "results/{date}/visual/report/feature",
            "results/{date}/visual/report/accuracy",
            "results/{date}/visual/report/volatility",
            "results/{date}/visual/report/empress-community",
            "results/{date}/visual/report/lme",
            "results/{date}/visual/report/rank-abundance/plots/",
            "results/{date}/visual/allfilter.html",
            "results/{date}/visual/report/sample.tsv",
            expand(
                "results/{{date}}/visual/beta-diversity-{metric}.html",
                metric=get_metric("beta"),
            ),
            expand(
                "results/{{date}}/visual/beta-diversity-{metric}.html",
                metric=get_phylogenetic_metric("beta"),
            ),
            expand(
                "results/{{date}}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}",
                metric=get_metric("beta"),
                metadata_column=get_metadata_columns(),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}",
                metric=get_phylogenetic_metric("beta"),
                metadata_column=get_metadata_columns(),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/beta-significance-{metric}-{diversity}-{metadata_column}",
                metric=get_phylogenetic_metric("beta"),
                metadata_column=get_metadata_categorical_columns(),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/beta-significance-{metric}-{diversity}-{metadata_column}",
                metric=get_metric("beta"),
                metadata_column=get_metadata_categorical_columns(),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-correlation-{metric_alpha}-{diversity}",
                metric_alpha=get_metric("alpha"),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-correlation-{metric_alpha}-{diversity}",
                metric_alpha=get_phylogenetic_metric("alpha"),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-significance-{metric_alpha}-{diversity}",
                metric_alpha=get_metric("alpha"),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/alpha-significance-{metric_alpha}-{diversity}",
                metric_alpha=get_phylogenetic_metric("alpha"),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/emperor-{metric}-{diversity}",
                metric=get_phylogenetic_metric("beta"),
                diversity="phylogenetic",
            ),
            expand(
                "results/{{date}}/visual/report/emperor-{metric}-{diversity}",
                metric=get_metric("beta"),
                diversity="normal",
            ),
            expand(
                "results/{{date}}/visual/report/ancom-{metadata_column}",
                metadata_column=get_ancom_columns(),
            ),
            expand(
                "results/{{date}}/visual/beta-diversity-{metric}.html",
                metric=get_complete_beta_metric(),
            ),
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
        "results/{date}/visual/report/multiqc.html",
        "results/{date}/visual/heatmap_binary.html",
        "results/{date}/visual/report/beta-rarefaction.svg",
        "results/{date}/visual/report/heatmap.svg",
        "results/{date}/visual/report/taxonomy.tsv",
        "results/{date}/out/report.zip",
        #"results/{date}/visual/fastq_stats.qzv",
        "results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
        "results/{date}/visual/absolute-taxabar-plot.html",
        "results/{date}/out/kraken.tar.gz",
        "results/{date}/out/qurro_plot/",
        "results/{date}/visual/report/rank-abundance/plots/",
        "results/{date}/visual/allfilter.html",
        expand(
            "results/{{date}}/visual/beta-diversity-{metric}.html",
            metric=get_metric("beta"),
        ),
        expand(
            "results/{{date}}/visual/beta-diversity-{metric}.html",
            metric=get_phylogenetic_metric("beta"),
        ),
        expand(
            "results/{{date}}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}",
            metric=get_metric("beta"),
            metadata_column=get_metadata_columns(),
            diversity="normal",
        ),
        expand(
            "results/{{date}}/visual/report/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}",
            metric=get_phylogenetic_metric("beta"),
            metadata_column=get_metadata_columns(),
            diversity="phylogenetic",
        ),
        "results/{date}/out/songbird/",
        "results/{date}/out/differentials_taxonomy.tsv",
        "results/{date}/out/config_parameters.html",
    output:
        "results/{date}/16S-report.tar.gz",
    params:
        outpath=config["output"],
    log:
        "logs/{date}/outputs/zip-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        """
        mkdir results/{wildcards.date}/16S-report
        cp -r {input} results/{wildcards.date}/16S-report/
        tar -czvf results/{wildcards.date}/16S-report.tar.gz results/{wildcards.date}/16S-report/
        cp results/{wildcards.date}/16S-report.tar.gz {params.outpath}
        rm -r results/{wildcards.date}/16S-report
        """
