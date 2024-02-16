rule chimera_only:
    input:
        table="results/{date}/out/table-cluster.qza",
        seqs="results/{date}/out/seq-cluster.qza",
        direc="results/{date}/out/uchime-dn-out",
    output:
        table="results/{date}/out/table-chimeras.qza",
        seqs="results/{date}/out/rep-seqs-chimeras.qza",
    log:
        "logs/{date}/chimera-checking/chimera-only.log",
    conda:
        "../envs/qiime-only-env.yaml"  #"../envs/qiime-chimerafilter.yaml"
    shell:
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file {input.direc}/chimeras.qza "
        "--p-no-exclude-ids "
        "--o-filtered-table {output.table} \n"
        "qiime feature-table filter-seqs "
        "--i-data {input.seqs} "
        "--m-metadata-file {input.direc}/chimeras.qza "
        "--p-no-exclude-ids "
        "--o-filtered-data {output.seqs}"


rule chimera_taxonomy:
    input:
        query="results/{date}/out/rep-seqs-chimeras.qza",
        reference_reads="resources/ref-seqs.qza",
        reference_taxonomy="resources/ref-taxa.qza",
    output:
        "results/{date}/out/chimera_taxonomy.qza",
    params:
        perc_identity=config["classification"]["perc-identity"],
        maxaccepts=config["classification"]["maxaccepts"],
        maxrejects=config["classification"]["maxrejects"],
    log:
        "logs/{date}/chimera-checking/chimera-taxonomy.log",
    conda:
        "../envs/qiime-only-env.yaml"  #"../envs/qiime-classifiers.yaml"
    shell:
        "qiime feature-classifier classify-consensus-vsearch "
        "--i-query {input.query} "
        "--i-reference-reads {input.reference_reads} "
        "--i-reference-taxonomy {input.reference_taxonomy} "
        "--p-maxaccepts {params.maxaccepts} "
        "--p-maxrejects {params.maxrejects} "
        "--p-perc-identity {params.perc_identity} "
        "--p-threads 4 "
        "--o-classification {output} "
        "--verbose"


rule taxa_barplot_chimera:
    input:
        table="results/{date}/out/table-chimeras.qza",  #"results/{date}/out/table-cluster-filtered.qza",
        taxonomy="results/{date}/out/chimera_taxonomy.qza",
    output:
        "results/{date}/visual/taxa-bar-plots-chimeras.qzv",
    log:
        "logs/{date}/chimera-checking/taxa-barplot-chimera.log",
    conda:
        "../envs/qiime-only-env.yaml"  #"../envs/qiime-taxonomy.yaml"
    shell:
        "qiime taxa barplot "
        "--i-table {input.table} "
        "--i-taxonomy {input.taxonomy} "
        "--m-metadata-file config/pep/sample.tsv "
        "--o-visualization {output}"
