rule get_database:
    output:
        genomic=temp("resources/ref-genome.fna.gz"),
        kraken=temp("resources/filtering-database.tgz"),
    params:
        genomic=str(config["database"]["ref-genome"]),
        kraken=str(config["database"]["kraken-db"]),
    log:
        "logs/prep_database.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cd resources; "
        "wget -O ref-genome.fna.gz {params.genomic}; "
        "wget -O filtering-database.tgz {params.kraken}; "


if (
    config["SILVA"] == True
    config["Greengenes"] == False
    config["NCBI"] == False
):

    rule get_SILVA:
        output:
            seq="resources/ref-seqs.qza",
            tax="resources/ref-taxa.qza",
        params:
            seq=str(config["database"]["download-path-seq"]),
            tax=str(config["database"]["download-path-tax"]),
        log:
            "logs/prep_SILVA.log",
        conda:
            "../envs/python.yaml"
        shell:
            "cd resources; "
            "wget -O ref-seqs.qza {params.seq}; "
            "wget -O ref-taxa.qza {params.tax}; "


if (
    config["Greengenes"] == True
    config["SILVA"] == False
    config["NCBI"] == False
):

    rule get_greengenes:
        output:
            seq="resources/ref-seqs.qza",
            tax="resources/ref-taxa.qza",
        params:
            seq=str(config["database"]["gg2-seq"]),
            tax=str(config["database"]["gg2-tax"]),
        log:
            "logs/prep_Greengenes.log",
        conda:
            "../envs/python.yaml"
        shell:
            "cd resources; "
            "wget -O ref-seqs.qza {params.seq}; "
            "wget -O ref-taxa.qza {params.tax}; "


if (
    config["NCBI"] == True
    config["NCBI"] == False
    config["Greengenes"] == False
):

    rule get_NCBI_ref:
        output:
            seq="resources/ref-seqs.qza",
            tax="resources/ref-taxa.qza",
        params:
            query=config["database"]["NCBI-query"],
        log:
            "logs/prep_NCBI.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime rescript get-ncbi-data "
            "--p-query {params.query} "
            "--o-sequences {output.seq} "
            "--o-taxonomy {output.tax} "
            "--verbose 2> {log} "


# rule get_SILVA:
#    output:
#        seq_rna=temp("resources/silva-138.1-ssu-nr99-rna-seqs.qza"),
#        tax=temp("resources/silva-138-99-tax.qza"),
#    params:
#        version="138.1",
#        target="SSURef_NR99",
#    log:
#        "logs/prerp_SILVA.log",
#    conda:
#        "../envs/qiime-only-env.yaml"
#    shell:
#        "qiime rescript get-silva-data "
#        "--p-version {params.version} "
#        "--p-target {params.target} "
#        "--p-include-species-labels "
#        "--o-silva-sequences {output.seq_rna} "
#        "--o-silva-taxonomy {output.tax} "
#        "2> {log}"


# rule rna_to_dna_SILVA:
#    input:
#        "resources/silva-138.1-ssu-nr99-rna-seqs.qza",
#   output:
#       temp("resources/silva-138-99-seqs.qza"),
#   log:
#       "logs/prerp_SILVA_toDNA.log",
#   conda:
#       "../envs/qiime-only-env.yaml"
#   shell:
#       "qiime rescript reverse-transcribe "
#       "--i-rna-sequences {input} "
#       "--o-dna-sequences {output} "
#       "2> {log}"


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


rule unzip_kraken:
    input:
        "resources/filtering-database.tgz",
    output:
        temp(directory("resources/filtering-database")),
    log:
        "logs/unzip_kraken_db.log",
    conda:
        "../envs/python.yaml"
    shell:
        "tar -zxvf {input} -C resources;"


if config["bowtie"] == False:

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


if (
    config["bowtie"] == True
    and config["datatype"] == "SampleData[PairedEndSequencesWithQuality]"
):

    rule read_samples:
        input:
            expand(
                "results/{{date}}/bowtie/data/{names}_L001_R1_001.fastq.gz",
                names=get_reads_for_kraken(),
            ),
            expand(
                "results/{{date}}/bowtie/data/{names}_L001_R2_001.fastq.gz",
                names=get_reads_for_kraken(),
            ),
            tsv="config/pep/sample.tsv",
            info="config/pep/sample_info.txt",
        output:
            temp("results/{date}/out/demux-paired-end.qza"),
        params:
            dirc="results/{date}/bowtie/data/",
            datatype=config["datatype"],
        log:
            "logs/{date}/preprocessing/read-samples.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime tools import "
            "--type {params.datatype} "
            "--input-path {params.dirc} "
            "--input-format CasavaOneEightSingleLanePerSampleDirFmt "
            "--output-path {output} "
            "2> {log} "


if (
    config["bowtie"] == True
    and config["datatype"] == "SampleData[SequencesWithQuality]"
):

    rule read_samples:
        input:
            expand(
                "results/{{date}}/bowtie/data/{names}_L001_R1_001.fastq.gz",
                names=get_reads_for_kraken(),
            ),
            tsv="config/pep/sample.tsv",
            info="config/pep/sample_info.txt",
        output:
            temp("results/{date}/out/demux-paired-end.qza"),
        params:
            dirc="results/{date}/bowtie/data/",
            datatype=config["datatype"],
        log:
            "logs/{date}/preprocessing/read-samples.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime tools import "
            "--type {params.datatype} "
            "--input-path {params.dirc} "
            "--input-format CasavaOneEightSingleLanePerSampleDirFmt "
            "--output-path {output} "
            "2> {log} "


if config["DADA2"] == False:

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
                --p-cores 10 \
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


if (
    config["datatype"] == "SampleData[PairedEndSequencesWithQuality]"
    and config["DADA2"] == False
):

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
