if config["DADA2"] == False:

    rule dereplication:
        input:
            "results/{date}/out/demux-joined-filtered.qza",
        output:
            table="results/{date}/out/derepl-table.qza",
            seqs="results/{date}/out/derepl-seq.qza",
        log:
            "logs/{date}/classification/dereplication.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime vsearch dereplicate-sequences "
            "--i-sequences {input} "
            "--o-dereplicated-table {output.table} "
            "--o-dereplicated-sequences {output.seqs} "
            "--verbose 2> {log}"

    rule de_novo_clustering:
        input:
            table="results/{date}/out/derep-table-nonhum.qza",
            seqs="results/{date}/out/derep-seq-nonhum.qza",
        output:
            table="results/{date}/out/table-cluster.qza",
            seqs="results/{date}/out/seq-cluster.qza",
        params:
            perc_identity=config["clustering"]["perc-identity"],
            threads=config["threads"],
        log:
            "logs/{date}/classification/de-novo-clustering.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime vsearch cluster-features-de-novo "
            "--i-table {input.table} "
            "--i-sequences {input.seqs} "
            "--p-perc-identity {params.perc_identity} "
            "--p-threads {params.threads} "
            "--o-clustered-table {output.table} "
            "--o-clustered-sequences {output.seqs} "
            "--verbose 2> {log}"


rule classification:
    input:
        query="results/{date}/out/seq-cluster-filtered.qza",
        reference_reads="resources/ref-seqs.qza",
        reference_taxonomy="resources/ref-taxa.qza",
    output:
        tax="results/{date}/out/taxonomy.qza",
        search="results/{date}/out/blast-search-results.qza",
    params:
        perc_identity=config["classification"]["perc-identity"],
        maxaccepts=config["classification"]["maxaccepts"],
        maxrejects=config["classification"]["maxrejects"],
        threads=config["threads"],
        query_cov=config["classification"]["query-cov"],
        min_consensus=config["classification"]["min-consensus"],
    log:
        "logs/{date}/classification/classification.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-classifier classify-consensus-vsearch "
        "--i-query {input.query} "
        "--i-reference-reads {input.reference_reads} "
        "--i-reference-taxonomy {input.reference_taxonomy} "
        "--p-maxaccepts {params.maxaccepts} "
        "--p-maxrejects {params.maxrejects} "
        "--p-perc-identity {params.perc_identity} "
        "--p-query-cov {params.query_cov} "
        "--p-min-consensus {params.min_consensus} "
        "--p-threads {params.threads} "
        "--o-classification {output.tax} "
        "--o-search-results {output.search} "
        "--verbose 2> {log} "


rule phylogenetic_tree:
    input:
        "results/{date}/out/seq-taxa-filtered.qza",  #"results/{date}/out/seq-cluster-filtered.qza"
    output:
        alignment="results/{date}/out/aligned-rep-seqs.qza",
        masked_alignment="results/{date}/out/masked-aligned-rep-seqs.qza",
        tree="results/{date}/visual/unrooted-tree.qza",
        rooted_tree="results/{date}/visual/rooted-tree.qza",
    log:
        "logs/{date}/classification/phylogenetic-tree.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime phylogeny align-to-tree-mafft-fasttree "
        "--i-sequences {input} "
        "--p-n-threads auto "
        "--o-alignment {output.alignment} "
        "--o-masked-alignment {output.masked_alignment} "
        "--o-tree {output.tree} "
        "--o-rooted-tree {output.rooted_tree} "
        "--verbose 2> {log}"


rule repeat_rarefaction:
    input:
        table="results/{date}/out/taxa_collapsed.qza",
    output:
        table="results/{date}/out/average-rarefied-table.qza",
    params:
        sampling_depth=config["rarefaction"]["sampling_depth"],
        repeats=config["rarefaction"]["repeats"],
    log:
        "logs/{date}/classification/repeat_rarefaction.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime repeat-rarefy repeat-rarefy "
        "--i-table {input.table} "
        "--p-sampling-depth {params.sampling_depth} "
        "--p-repeat-times {params.repeats} "
        "--o-rarefied-table {output.table} "
        "--verbose 2> {log}"


rule relative_collapsed_taxa:
    input:
        "results/{date}/out/taxa_collapsed.qza",
    output:
        "results/{date}/out/taxa_collapsed_relative.qza",
    log:
        "logs/{date}/outputs/taxa-collapse-relative.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-table relative-frequency "
        "--i-table {input} "
        "--o-relative-frequency-table {output} "
        "--verbose 2> {log} "
