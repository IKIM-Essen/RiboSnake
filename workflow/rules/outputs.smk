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
            category="3. Rarefaction",
            subcategory="Beta rarefaction",
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
            category="3. Rarefaction",
            subcategory="Beta",
            htmlindex="index.html",
        ),
        alpha_html=report(
            directory("results/{date}/visual/report/alpha_rarefaction"),
            caption="../report/alpha-rarefaction.rst",
            category="3. Rarefaction",
            subcategory="Alpha",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-files.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_reports.py"


rule snakemake_report:
    input:
        "results/{date}/visual/heatmap_binary.png",
        "results/{date}/visual/report/beta-rarefaction.svg",
        "results/{date}/visual/report/heatmap.svg",
        "results/{date}/visual/unzipped",
        "results/{date}/visual/report/multiqc.html",
    output:
        "results/{date}/out/report.zip",
    params:
        for_testing=get_if_testing("--snakefile ../workflow/Snakefile"),
    log:
        "logs/{date}/outputs/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} --report-stylesheet resources/custom-stylesheet.css {input}"
        " --report {output} {params.for_testing}"


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
        #"results/{date}/visual/demux-joined-filter-stats.qzv",
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
