rule get_database:
    output:
        genomic=temp("resources/ref-genome.fna.gz"),
        #kraken=temp("resources/filtering-database.tgz"),
    params:
        genomic=str(config["database"]["ref-genome"]),
        #kraken=str(config["database"]["kraken-db"]),
    log:
        "logs/prep_database.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cd resources; "
        "wget -O ref-genome.fna.gz {params.genomic}; "
        #"wget -O filtering-database.tgz {params.kraken}; "

rule unzip_ref_gen:
    input:
        "resources/ref-genome.fna.gz",
    output:
        temp("resources/ref-genome.fna"),
    log:
        "logs/unzip_ref_gen.log",
    conda:
        "../envs/python.yaml"
    shell:
        "gzip -dc {input} > {output} 2> {log}; "

rule lower_to_upper:
    input:
        "resources/ref-genome.fna",
    output:
        temp("resources/ref-genome_upper.fna"),
    log:
        "logs/lower_to_upper_fasta.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/lower_to_upper_fasta.py"

rule import_ref_genome:
    input:
        "resources/ref-genome_upper.fna",
    output:
        temp("resources/ref-genome_upper.qza"),
    log:
        "logs/import_ref_gen.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools import "
        "--input-path {input} "
        "--output-path {output} "
        "--type 'FeatureData[Sequence]' "
        "2> {log} "

#rule get_SILVA:
#    output:
#        seq="resources/ref-seqs.qza",
#        tax="resources/ref-taxa.qza",
#    params:
#        seq=str(config["database"]["download-path-seq"]),
#        tax=str(config["database"]["download-path-tax"]),
#    log:
#        "logs/prep_SILVA.log",
#    conda:
#        "../envs/python.yaml"
#    shell:
#        "cd resources; "
#        "wget -O ref-seqs.qza {params.seq}; "
#        "wget -O ref-taxa.qza {params.tax}; "

rule read_samples:
    input:
        tsv="config/pep/sample.tsv",
        info="config/pep/sample_info.txt",
    output:
        temp("results/{date}/out/demux-paired-end.qza"),
    params:
        direc=get_data_dir(),
        datatype=config["datatype"],
    log:
        "logs/{date}/preprocessing/read-samples.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools import "
        "--type {params.datatype} "
        "--input-path {params.direc} "
        "--input-format CasavaOneEightSingleLanePerSampleDirFmt "
        "--output-path {output} "
        "2> {log} "

rule trim_paired:
    input:
        "results/{date}/out/demux-paired-end.qza",
    output:
        temp("results/{date}/out/trimmed-seqs.qza"),
    params:
        datatype=config["datatype"],
        adapter1=config["adapter1"],
        adapter2=config["adapter2"],
        primer1=config["primertrimming"]["forward"],
        primer2=config["primertrimming"]["reverse"],
        error_rate=config["primertrimming"]["error_rate"],
        rep_times=config["primertrimming"]["rep_times"],
        overlap=config["primertrimming"]["overlap"],
        min_length=config["primertrimming"]["min_length"],
        threads=config["threads"],
    log:
        "logs/{date}/preprocessing/trim-paired.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        """
        if [[ '${params.datatype}' == '$SampleData[PairedEndSequencesWithQuality]' ]] 
        then 
            qiime cutadapt trim-paired \
            --i-demultiplexed-sequences {input} \
            --p-cores {params.threads} \
            --p-adapter-f {params.adapter1} \
            --p-front-f {params.primer1} \
            --p-front-r {params.primer2} \
            --p-adapter-r {params.adapter2} \
            --p-error-rate {params.error_rate} \
            --p-times {params.rep_times} \
            --p-overlap {params.overlap} \
            --p-minimum-length {params.min_length} \
            --o-trimmed-sequences {output} \
            --verbose 2> {log}
        else 
            qiime cutadapt trim-single \
            --i-demultiplexed-sequences {input} \
            --p-cores {params.threads} \
            --p-adapter {params.primer1} \
            --p-front {params.primer2} \
            --p-error-rate {params.error_rate} \
            --p-times {params.rep_times} \
            --p-overlap {params.overlap} \
            --p-minimum-length {params.min_length} \
            --o-trimmed-sequences {output} \
            --verbose 2> {log}
        fi
        """

rule join_ends:
    input:
        "results/{date}/out/trimmed-seqs.qza",
    output:
        temp("results/{date}/out/joined-seqs.qza"),
    params:
        minovlen=config["sequence_joining"]["seq_join_length"],
        minlen=config["sequence_joining"]["minlen"],
        maxdiffs=config["sequence_joining"]["maxdiffs"],
        #qmin=config["sequence_joining"]["qmin"],
        #qminout=config["sequence_joining"]["qminout"],
        #qmax=config["sequence_joining"]["qmax"],
        #qmaxout=config["sequence_joining"]["qmaxout"],
        threads=config["sequence_joining"]["threads"],
    log:
        "logs/{date}/preprocessing/join-ends.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime vsearch merge-pairs "
        "--i-demultiplexed-seqs {input} "
        "--p-allowmergestagger "
        "--p-minovlen {params.minovlen} "
        "--p-minlen {params.minlen} "
        "--p-maxdiffs {params.maxdiffs} "
        "--p-threads {params.threads} "
        "--o-merged-sequences {output} "
        "--verbose 2> {log}"

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
        "--o-filtered-data {output.seqs} "
        "--verbose 2> {log}"

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
        "--verbose 2> {log}"

rule visualize_beforeChimera:
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

rule filter_human:
    input:
        seq="results/{date}/out/derepl-seq.qza",
        table="results/{date}/out/derepl-table.qza",
        ref_seq="resources/ref-genome_upper.qza",
    output:
        seq="results/{date}/out/derep-seq-nonhum.qza",
        table="results/{date}/out/derep-table-nonhum.qza",
        human_hit="results/{date}/out/human.qza",
    params:
        threads=config["threads"],
        perc_identity=config["filtering"]["perc-identity"],
        perc_query_aligned=config["filtering"]["perc-query-aligned"],
    threads: config["threads"]
    log:
        "logs/{date}/filtering/filter-human.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime quality-control exclude-seqs "
        "--i-query-sequences {input.seq} "
        "--i-reference-sequences {input.ref_seq} "
        "--p-threads {params.threads} "
        "--p-perc-identity {params.perc_identity} "
        "--p-perc-query-aligned {params.perc_query_aligned} "
        "--o-sequence-hits {output.human_hit} "
        "--o-sequence-misses {output.seq} "
        "--verbose 2> {log} \n"
        "qiime feature-table filter-features "
        "--i-table {input.table} "
        "--m-metadata-file {output.seq} "
        "--o-filtered-table {output.table} "
        "--verbose 2> {log} "

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

rule visualize_trimmed:
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

rule export_taxa_collapsed_relative:
    input:
        "results/{date}/out/taxa_collapsed_relative.qza",
    output:
        directory("results/{date}/visual/report/taxa_collapsed_relative/"),
    log:
        "logs/{date}/visualisation/export_taxa_collapsed_relative.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "2> {log}"

rule convert_taxa_collapsed_relative_tsv:
    input:
        "results/{date}/visual/report/taxa_collapsed_relative/",
    output:
        report(
            "results/{date}/visual/report/taxa_collapsed_relative.tsv",
            caption="../report/relative-taxa.rst",
            category="1. Abundances",
        ),
    params:
        export_dir="results/{date}/visual/report/taxa_collapsed_relative/",
    log:
        "logs/{date}/visualisation/taxa_collapsed_relative.log",
    conda:
        "../envs/python.yaml"
    shell:
        "biom convert "
        "-i {params.export_dir}/feature-table.biom "
        "-o {output} "
        "--to-tsv --header-key taxonomy "
        "2>> {log}"

rule export_taxa_collapsed_absolute:
    input:
        "results/{date}/out/taxa_collapsed.qza",
    output:
        directory("results/{date}/visual/report/taxa_collapsed_absolute/"),
    log:
        "logs/{date}/visualisation/export_taxa_collapsed_absolute.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools export "
        "--input-path {input} "
        "--output-path {output} "
        "2> {log}"

rule convert_taxa_collapsed_absolute_tsv:
    input:
        "results/{date}/visual/report/taxa_collapsed_absolute/",
    output:
        report(
            "results/{date}/visual/report/taxa_collapsed_absolute.tsv",
            caption="../report/absolute-taxa.rst",
            category="1. Abundances",
        ),
    log:
        "logs/{date}/visualisation/convert_taxa_collapsed_absolute.log",
    conda:
        "../envs/python.yaml"
    shell:
        "biom convert "
        "-i {input}/feature-table.biom "
        "-o {output} "
        "--to-tsv --header-key taxonomy "
        "2> {log}"

rule separate_samples:
    input:
        abs="results/{date}/visual/report/taxa_collapsed_absolute.tsv",
        rel="results/{date}/visual/report/taxa_collapsed_relative.tsv",
    output:
        directory("results/{date}/visual/report/sep_sample"),
    log:
        "logs/{date}/visualisation/separate_samples.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/separate_samples.py"


rule compress_samples:
    input:
        "results/{date}/visual/report/sep_sample"
    output:
        report(
            "results/{date}/visual/report/sep_sample.tar.gz",
            caption="../report/all-filter.rst",
            category="1. Abundances",
        )
    log:
        "logs/{date}/outputs/separate-samples.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        tar -czf {output} {input}
        """


rule visualize_samples:
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

rule visualize_joined:
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

rule all_filter:
    input:
        samples="results/{date}/visual/paired-seqs",
        trimmed="results/{date}/visual/trimmed-seqs",
        joined="results/{date}/visual/joined-seqs",
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
            category="1. Abundances",
        ),
    log:
        "logs/{date}/visualisation/all-filter.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/complete_filter.py"

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
        "results/{date}/visual/paired-seqs.qzv",
        "results/{date}/visual/fastq_stats.qzv",
        "results/{date}/visual/table-whuman.qzv",
        "results/{date}/visual/table-wohuman.qzv",
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
        paired_seqs=report(
            directory("results/{date}/visual/report/paired-seqs"),
            caption="../report/paired-seqs.rst",
            category="2. Qualitycontrol",
            htmlindex="index.html",
        ),
        demux_filter_stats=report(
            directory("results/{date}/visual/report/demux-joined-filter-stats"),
            caption="../report/demux-filter-stats.rst",
            category="2. Qualitycontrol",
            htmlindex="index.html",
        ),
        fastq_stats=report(
            directory("results/{date}/visual/report/fastq_stats"),
            caption="../report/fastq-stats.rst",
            category="2. Qualitycontrol",
            htmlindex="index.html",
        ),
    log:
        "logs/{date}/outputs/report-files.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_reports.py"

rule include_metadata:
    input:
        "config/pep/sample.tsv",
    output:
        report(
            "results/{date}/visual/report/sample.tsv",
            caption="../report/metadata.rst",
            category="2. Qualitycontrol",
        ),
    log:
        "logs/{date}/visualisation/metadata.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cp {input} {output}"

rule zip_report:
    input:
        "results/{date}/visual/table-cluster-lengthfilter.qzv",
        "results/{date}/visual/fastq_stats.qzv",
        "results/{date}/visual/allfilter.html",
        "results/{date}/visual/report/sep_sample.tar.gz",
        report="results/{date}/out/report.zip",
    output:
        "results/{date}/{date}.tar.gz",
    params:
        outpath=config["output"],
    log:
        "logs/{date}/outputs/zip-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        """
        set -euo pipefail
        exec > {log} 2>&1

        echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] Starting zip_report for {wildcards.date}"
        mkdir -p results/{wildcards.date}/16S-report/
        mkdir -p results/{wildcards.date}/16S-report/additional/
        cp -r {input} results/{wildcards.date}/16S-report/additional/
        rm -f results/{wildcards.date}/16S-report/additional/report.zip || true
        cp {input.report} results/{wildcards.date}/16S-report/
        tar -czvf {output} results/{wildcards.date}/16S-report/
        cp {output} {params.outpath}
        rm -r results/{wildcards.date}/16S-report
        echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] Finished zip_report for {wildcards.date}"
        """

rule snakemake_report:
    input:
        "results/{date}/visual/unzipped",
        "results/{date}/visual/allfilter.html",
        "results/{date}/visual/report/sample.tsv",
        "results/{date}/visual/report/sep_sample",
        "results/{date}/visual/report/sep_sample.tar.gz",
    output:
        "results/{date}/out/report.zip",
    params:
        for_testing=get_if_testing("--snakefile ../workflow/Snakefile"),
    log:
        "logs/{date}/outputs/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --profile '' --nolock --report {output} "
        "{params.for_testing} "
        "> {log} 2>&1"

rule concatenate_logs:
    input:
        "results/{date}/{date}.tar.gz",
    conda:
        "../envs/python.yaml"
    log:
        "logs/{date}/outputs/logs.log",
    output:
        "logs/{date}_logs.tar.gz",
    shell:
        """
        tar -czvf {output} logs/{wildcards.date}/
        rm -r logs/{wildcards.date}
        """