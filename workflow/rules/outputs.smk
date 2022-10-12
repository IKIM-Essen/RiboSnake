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
        "results/{date}/visual/alpha-rarefaction.qzv",
        "results/{date}/visual/beta-rarefaction.qzv",
        "results/{date}/visual/heatmap.qzv",
        "results/{date}/visual/taxa-bar-plots.qzv",
        "results/{date}/visual/rooted-tree.qza",
        "results/{date}/visual/taxonomy.qzv",
        "results/{date}/visual/faith-pd-group-significance.qzv",
        "results/{date}/visual/evenness-group-significance.qzv",
        expand(
            "results/{{date}}/visual/unweighted-unifrac-significance-{metadata_column}.qzv",
            metadata_column=get_metadata_categorical_columns(),
        ),
        "results/{date}/visual/bray-curtis-emperor.qzv",
        "results/{date}/visual/jaccard-emperor.qzv",
        "results/{date}/visual/unweighted-unifrac-emperor.qzv",
        "results/{date}/visual/weighted-unifrac-emperor-plot.qzv",
        expand(
            "results/{{date}}/visual/beta-correlation-scatter-{metadata_column}.qzv",
            metadata_column=get_metadata_columns(),
        ),
        "results/{date}/visual/heatmap_gneiss.qzv",
        expand(
            "results/{{date}}/visual/ancom-{metadata_column}.qzv",
            metadata_column=get_ancom_columns(),
        ),
        "results/{date}/visual/alpha_correlation.qzv",
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
        alpha_significance=report(
            directory("results/{date}/visual/report/evenness-group-significance"),
            caption="../report/alpha-significance.rst",
            category="3. Analysis",
            subcategory="Alpha",
            htmlindex="index.html",
        ),
        faith_pd=report(
            directory("results/{date}/visual/report/faith-pd-group-significance"),
            caption="../report/faith-pd-significance.rst",
            category="3. Analysis",
            subcategory="Alpha",
            htmlindex="index.html",
        ),
        bray_curtis_emperor=report(
            directory("results/{date}/visual/report/bray-curtis-emperor"),
            caption="../report/bray-curtis-emperor.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
        jaccard_emperor=report(
            directory("results/{date}/visual/report/jaccard-emperor"),
            caption="../report/jaccard-emperor.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
        unweighted_unifrac_emperor=report(
            directory("results/{date}/visual/report/unweighted-unifrac-emperor"),
            caption="../report/unweighted-unifrac-emperor.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
        weighted_unifrac_emperor=report(
            directory("results/{date}/visual/report/weighted-unifrac-emperor"),
            caption="../report/weighted-unifrac-emperor.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
        gneiss=report(
            "results/{date}/visual/heatmap_gneiss.svg",
            caption="../report/gneiss.rst",
            category="3. Analysis",
            subcategory="Gneiss",
        ),
        alpha_correlation=report(
            directory("results/{date}/visual/report/alpha_correlation"),
            caption="../report/alpha_correlation.rst",
            category="3. Analysis",
            subcategory="Alpha",
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
                "results/{date}/visual/report/beta-correlation-scatter-{metadata_column}"
            ),
            caption="../report/beta-correlation-scatter.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-beta-correaltion-{metadata_column}.log",
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
                "results/{date}/visual/report/unweighted-unifrac-significance-{metadata_column}"
            ),
            caption="../report/beta-significance.rst",
            category="3. Analysis",
            subcategory="Beta",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-significance-{metadata_column}.log",
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
        "results/{date}/visual/report/beta-rarefaction.svg",
        "results/{date}/visual/report/heatmap.svg",
        "results/{date}/visual/unzipped",
        "results/{date}/visual/report/multiqc.html",
        "results/{date}/visual/absolute-taxabar-plot.png",
        expand(
            "results/{{date}}/visual/report/beta-correlation-scatter-{metadata_column}",
            metadata_column=get_metadata_columns(),
        ),
        expand(
            "results/{{date}}/visual/report/unweighted-unifrac-significance-{metadata_column}",
            metadata_column=get_metadata_categorical_columns(),
        ),
        expand(
            "results/{{date}}/visual/report/ancom-{metadata_column}",
            metadata_column=get_ancom_columns(),
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
        "results/{date}/visual/report/beta-rarefaction.svg",
        "results/{date}/visual/report/heatmap.svg",
        "results/{date}/visual/report/taxonomy.tsv",
        "results/{date}/out/report.zip",
        "results/{date}/visual/fastq_stats.qzv",
        "results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
        "results/{date}/visual/absolute-taxabar-plot.png",
        "results/{date}/out/kraken.tar.gz",
        #"results/{date}/out/alpha-diversity.qza",
        "results/{date}/out/beta-diversity-distance.qza",
        "results/{date}/out/parameter-summary.csv",
        expand(
            "results/{{date}}/out/beta-correlation-{metadata_column}.qza",
            metadata_column=get_metadata_columns(),
        ),
        #"results/{date}/out/beta-phylogeny.qza",
        #"results/{date}/out/alpha-phylogeny.qza",
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
