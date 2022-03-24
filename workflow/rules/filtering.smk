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
        "--o-filter-stats {output.stats}"


rule chimera_filtering:
    input:
        table="results/{date}/out/table-cluster.qza",
        seqs="results/{date}/out/seq-cluster.qza",
    output:
        direc=directory("results/{date}/out/uchime-dn-out"),
        table="results/{date}/out/table-nonchimeric-wo-borderline.qza",
        seqs="results/{date}/out/rep-seqs-nonchimeric-wo-borderline.qza",
    params:
        minh=config["filtering"]["chimera-minh"],
    log:
        "logs/{date}/filtering/chimera-filtering.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime vsearch uchime-denovo "
        "--i-table {input.table} "
        "--i-sequences {input.seqs} "
        "--p-minh {params.minh} "
        "--output-dir {output.direc} \n"
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file {output.direc}/chimeras.qza "
        "--p-exclude-ids "
        "--o-filtered-table {output.table} \n"
        "qiime feature-table filter-seqs "
        "--i-data {input.seqs} "
        "--m-metadata-file {output.direc}/chimeras.qza "
        "--p-exclude-ids "
        "--o-filtered-data {output.seqs}"


rule filter_seq_length:
    input:
        seq="results/{date}/out/rep-seqs-nonchimeric-wo-borderline.qza",  #results/{date}/out/seq-cluster.qza",
        table="results/{date}/out/table-nonchimeric-wo-borderline.qza",  #"results/{date}/out/table-cluster.qza"
    output:
        seq="results/{date}/out/seq-cluster-lengthfilter.qza",
        table="results/{date}/out/table-cluster-lengthfilter.qza",
    params:
        min_length=config["filtering"]["min-seq-length"],
    log:
        "logs/{date}/filtering/filter-seq-length.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-table filter-seqs "
        "--i-data {input.seq} "
        "--m-metadata-file {input.seq} "
        "--p-where 'length(sequence) > {params.min_length}' "
        "--o-filtered-data {output.seq} \n"
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file {output.seq} "
        "--o-filtered-table {output.table} "


rule abundance_frequency:
    input:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
    output:
        "results/{date}/out/abundance.txt",
    params:
        relative_abundance=config["filtering"]["relative-abundance-filter"],
    log:
        "logs/{date}/filtering/abundance-frequency.log",
    conda:
        "../envs/abundancefiltering.yaml"
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
        "--o-filtered-table {output.table} \n"
        "qiime feature-table filter-seqs "
        "--i-data {input.seqs} "
        "--i-table {output.table} "
        "--p-no-exclude-ids "
        "--o-filtered-data {output.seqs}"


# rule filter_abundance:
#    input:
#        table = "results/{date}/out/table-cluster-freq.qza",
#        seqs = "results/{date}/out/seq-cluster-freq.qza"
#    output:
#        table = "results/{date}/out/table-cluster-filtered.qza",
#        seqs = "results/{date}/out/seq-cluster-filtered.qza"
#    params:
#        abundance = float(0.015),
#        prevalence = float(0.001)
#    shell:
#        "qiime feature-table filter-features-conditionally "
#            "--i-table {input.table} "
#            "--p-abundance {params.abundance} "
#            "--p-prevalence {params.prevalence} "
#            "--o-filtered-table {output.table} \n"
#        "qiime feature-table filter-seqs "
#            "--i-data {input.seqs} "
#            "--i-table {output.table} "
#            "--p-no-exclude-ids "
#            "--o-filtered-data {output.seqs}"


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
        "--o-collapsed-table {output}"


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
        "--o-filtered-table {output.table} \n"
        "qiime taxa filter-seqs "
        "--i-sequences {input.seq} "
        "--i-taxonomy {input.taxonomy} "
        "--p-exclude mitochondria,chloroplast "
        "--o-filtered-sequences {output.seq} "
