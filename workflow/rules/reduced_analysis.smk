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


rule visualise_beforeChimera:
    input:
        "results/{date}/out/table-nonchimeric-wo-borderline.qza",
    output:
        "results/{date}/out/table-nonchimeric-wo-borderline.qzv",
    log:
        "logs/{date}/visualisation/visualise-chimera.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-table summarize "
        "--i-table {input} "
        "--o-visualization {output} "
        "--verbose 2> {log}"


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


if config["data-type"] == "human" and config["bowtie"] == False:

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


if config["data-type"] == "environmental" or config["bowtie"] == True:

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


# if config["DADA2"] == False:


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


if config["DADA2"] == True:

    rule visualize_dada2_stats:
        input:
            "results/{date}/out/dada2-stats.qza",
        output:
            "results/{date}/visual/dada2-stats-visual.qzv",
        log:
            "logs/{date}/outputs/dada2-stats-visual.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime metadata tabulate "
            "--m-input-file {input} "
            "--o-visualization {output} "
            "--verbose 2> {log}"

    rule unzip_reports:
        input:
            "results/{date}/visual/heatmap.qzv",
            "results/{date}/visual/taxa-bar-plots.qzv",
            "results/{date}/visual/rooted-tree.qza",
            "results/{date}/visual/taxonomy.qzv",
            "results/{date}/visual/paired-seqs.qzv",
            "results/{date}/visual/fastq_stats.qzv",
            "results/{date}/visual/dada2-stats-visual.qzv",
            "results/{date}/visual/empress-community.qzv",
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
            dada2=report(
                directory("results/{date}/visual/report/dada2-stats-visual"),
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

    rule snakemake_report:
        input:
            "results/{date}/visual/heatmap_binary.html",
            "results/{date}/visual/report/heatmap.svg",
            "results/{date}/visual/unzipped",
            "results/{date}/visual/report/multiqc.html",
            "results/{date}/visual/absolute-taxabar-plot.html",
            "results/{date}/visual/allfilter.html",
            "results/{date}/visual/report/sample.tsv",
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

    rule zip_report:
        input:
            "results/{date}/visual/table-cluster-lengthfilter.qzv",
            "results/{date}/visual/rooted-tree.qza",
            "results/{date}/out/taxonomy.qza",
            "results/{date}/out/aligned-rep-seqs.qza",
            "results/{date}/out/biom_table/",
            "results/{date}/out/binary_biom/",
            "results/{date}/visual/report/multiqc.html",
            "results/{date}/visual/heatmap_binary.html",
            "results/{date}/visual/report/heatmap.svg",
            "results/{date}/visual/report/taxonomy.tsv",
            "results/{date}/out/report.zip",
            "results/{date}/visual/fastq_stats.qzv",
            "results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
            "results/{date}/visual/absolute-taxabar-plot.html",
            "results/{date}/out/kraken.tar.gz",
            "results/{date}/out/config_parameters.html",
            "results/{date}/visual/report/rank-abundance/plots",
            "results/{date}/visual/allfilter.html",
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


if config["DADA2"] == False:

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
            "results/{date}/visual/empress-community.qzv",
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

    rule snakemake_report:
        input:
            "results/{date}/visual/heatmap_binary.html",
            "results/{date}/visual/report/heatmap.svg",
            "results/{date}/visual/unzipped",
            "results/{date}/visual/report/multiqc.html",
            "results/{date}/visual/absolute-taxabar-plot.html",
            "results/{date}/visual/allfilter.html",
            "results/{date}/visual/report/empress-community",
            "results/{date}/visual/report/sample.tsv",
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

    rule all_filter:
        input:
            first="results/{date}/visual/report/demux-joined-filter-stats/",
            human="results/{date}/visual/sample_frequencys_difference.csv",
            wo_chimera="results/{date}/visual/chimera_unzipped/",
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
            "results/{date}/visual/report/heatmap.svg",
            "results/{date}/visual/report/taxonomy.tsv",
            "results/{date}/out/report.zip",
            "results/{date}/visual/fastq_stats.qzv",
            "results/{date}/out/table.from_biom_w_taxonomy-featcount.txt",
            "results/{date}/visual/absolute-taxabar-plot.html",
            "results/{date}/out/kraken.tar.gz",
            "results/{date}/out/config_parameters.html",
            "results/{date}/visual/report/rank-abundance/plots",
            "results/{date}/visual/allfilter.html",
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


if config["bowtie"] == False:

    rule hum_filter_difference:
        input:
            "results/{date}/visual/unzipped/",
        output:
            "results/{date}/visual/sample_frequencys_difference.csv",
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
