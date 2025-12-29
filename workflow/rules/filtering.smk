if (
    config["datatype"] == "SampleData[PairedEndSequencesWithQuality]"
    and config["DADA2"] == False
):

    rule fastq_score:
        input:
            "results/{date}/out/joined-seqs.qza",
        output:
            filtering="results/{date}/out/demux-joined-filtered.qza",
            stats="results/{date}/out/demux-joined-filter-stats.qza",
        params:
            date=get_date(),
            min_quality=config["filtering"]["phred-score"],
            min_length_frac=config["filtering"]["min-length-frac"],
            max_ambig=config["filtering"]["max-ambiguity"],
        log:
            "logs/{date}/filtering/fastq-score.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime quality-filter q-score "
            "--i-demux {input} "
            "--p-min-quality {params.min_quality} "
            "--p-min-length-fraction {params.min_length_frac} "
            "--p-max-ambiguous {params.max_ambig} "
            "--o-filtered-sequences {output.filtering} "
            "--o-filter-stats {output.stats} "
            "--verbose 2> {log}"


if (
    config["datatype"] == "SampleData[SequencesWithQuality]"
    and config["DADA2"] == False
):

    rule fastq_score:
        input:
            "results/{date}/out/trimmed-seqs.qza",
        output:
            filtering="results/{date}/out/demux-joined-filtered.qza",
            stats="results/{date}/out/demux-joined-filter-stats.qza",
        params:
            date=get_date(),
            min_quality=config["filtering"]["phred-score"],
            min_length_frac=config["filtering"]["min-length-frac"],
            max_ambig=config["filtering"]["max-ambiguity"],
        log:
            "logs/{date}/filtering/fastq-score.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime quality-filter q-score "
            "--i-demux {input} "
            "--p-min-quality {params.min_quality} "
            "--p-min-length-fraction {params.min_length_frac} "
            "--p-max-ambiguous {params.max_ambig} "
            "--o-filtered-sequences {output.filtering} "
            "--o-filter-stats {output.stats} "
            "--verbose 2> {log}"


rule abundance_frequency:
    input:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
    output:
        abundance="results/{date}/out/abundance.txt",
        feature_table=directory("results/{date}/visual/table-cluster-lengthfilter/data"),
    params:
        relative_abundance=config["filtering"]["relative-abundance-filter"],
    log:
        "logs/{date}/filtering/abundance-frequency.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/relative_abundance.py"


rule filter_frequency:
    input:
        table="results/{date}/out/table-cluster-lengthfilter.qza",  #"results/{date}/out/table-cluster.qza", 
        seqs="results/{date}/out/seq-cluster-lengthfilter.qza",  #"results/{date}/out/seq-cluster.qza", 
        abundance="results/{date}/out/abundance.txt",
    output:
        table="results/{date}/out/table-cluster-filtered.qza",  # "results/{date}/out/table-cluster-freq.qza"
        seqs="results/{date}/out/seq-cluster-filtered.qza",  # "results/{date}/out/seq-cluster-freq.qza"
    log:
        "logs/{date}/filtering/filter-frequency.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "value=$(<{input.abundance}) \n"
        "echo $value \n"
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--p-min-frequency $value "
        "--o-filtered-table {output.table} "
        "--verbose 2> {log} \n"
        "qiime feature-table filter-seqs "
        "--i-data {input.seqs} "
        "--i-table {output.table} "
        "--p-no-exclude-ids "
        "--o-filtered-data {output.seqs} "
        "--verbose 2> {log} "


rule unzip_frequency:
    input:
        "results/{date}/visual/table-cluster-filtered.qzv",
    output:
        temp(directory("results/{date}/visual/frequency_unzipped")),
    log:
        "logs/{date}/outputs/unzip-frequency.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


if config["DADA2"] == False:

    rule unzip_frequency_chimera:
        input:
            "results/{date}/out/table-nonchimeric-wo-borderline.qzv",
        output:
            temp(directory("results/{date}/visual/chimera_unzipped")),
        log:
            "logs/{date}/outputs/unzip-chimera.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/rename_qzv.py"


if config["reduced-analysis"] == True:

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


rule frequency_after_abundancefilter:
    input:
        "results/{date}/visual/frequency_unzipped",
    output:
        directory("results/{date}/visual/report/table-cluster-filtered"),
    log:
        "logs/{date}/filtering/after_abundance-frequency.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_significance.py"
        

rule taxa_collapse:
    input:
        table="results/{date}/out/table-cluster-filtered.qza",
        taxonomy="results/{date}/out/taxonomy.qza",
    output:
        "results/{date}/out/taxa_collapsed.qza",
    log:
        "logs/{date}/filtering/taxa-collapse.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime taxa collapse "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--p-level 6 "
        "--o-collapsed-table {output} "
        "--verbose 2> {log} "


rule filter_taxonomy:
    input:
        table="results/{date}/out/table-cluster-filtered.qza",
        seq="results/{date}/out/seq-cluster-filtered.qza",
        taxonomy="results/{date}/out/taxonomy.qza",
    output:
        table="results/{date}/out/table-taxa-filtered.qza",
        seq="results/{date}/out/seq-taxa-filtered.qza",
    log:
        "logs/{date}/filtering/filter-taxonomy.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime taxa filter-table "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-table {output.table} "
        "--verbose 2> {log} \n"
        "qiime taxa filter-seqs "
        "--i-sequences {input.seq} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-sequences {output.seq} "
        "--verbose 2> {log} "
