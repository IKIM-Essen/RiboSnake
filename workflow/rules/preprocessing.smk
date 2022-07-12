rule get_database:
    output:
        seq="resources/silva-138-99-seqs.qza",
        tax="resources/silva-138-99-tax.qza",
        genomic="resources/GRCh38_latest_genomic.fna.gz",
        kraken="resources/minikraken2_v2_8GB_201904.tgz",
    params:
        seq=str(config["database"]["download-path-seq"]),
        tax=str(config["database"]["download-path-tax"]),
        genomic=str(config["database"]["ref-genome"]),
        kraken=str(config["database"]["kraken-db"]),
    log:
        "logs/prep_database.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cd resources; "
        "wget {params.seq}; "
        "wget {params.tax}; "
        "wget {params.genomic}; "
        "wget {params.kraken}; "


rule unzip_ref_gen:
    input:
        "resources/GRCh38_latest_genomic.fna.gz",
    output:
        fasta="resources/GRCh38_latest_genomic.fna",
    log:
        "logs/unzip_ref_gen.log",
    conda:
        "../envs/python.yaml"
    shell:
        "gzip -dc {input} > {output.fasta}; "


rule lower_to_upper:
    input:
        "resources/GRCh38_latest_genomic.fna",
    output:
        "resources/GRCh38_latest_genomic_upper.fna",
    log:
        "logs/lower_to_upper_fasta.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/lower_to_upper_fasta.py"


rule import_ref_genome:
    input:
        "resources/GRCh38_latest_genomic_upper.fna",
    output:
        "resources/GRCh38_latest_genomic_upper.qza",
    log:
        "logs/import_ref_gen.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools import "
        "--input-path {input} "
        "--output-path {output} "
        "--type 'FeatureData[Sequence]' "


rule unzip_kraken:
    input:
        "resources/minikraken2_v2_8GB_201904.tgz",
    output:
        directory("resources/minikraken2_v2_8GB_201904_UPDATE"),
    log:
        "logs/unzip_kraken_db.log",
    conda:
        "../envs/python.yaml"
    shell:
        "tar -zxvf {input} -C resources;"


rule read_samples:
    input:
        tsv="config/pep/sample.tsv",
        info="config/pep/sample_info.txt",
    output:
        "results/{date}/out/demux-paired-end.qza",
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
        "--output-path {output}"


rule trim_paired:
    input:
        "results/{date}/out/demux-paired-end.qza",
    output:
        "results/{date}/out/trimmed-seqs.qza",
    params:
        datatype=config["datatype"],
        primer1=config["primer1"],
        primer2=config["primer2"],
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
            --p-adapter-f {params.primer1} \
            --p-front-f {params.primer2} \
            --p-error-rate {params.error_rate} \
            --p-times {params.rep_times} \
            --p-overlap {params.overlap} \
            --p-minimum-length {params.min_length} \
            --o-trimmed-sequences {output} 
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
            --o-trimmed-sequences {output}
        fi
        """


if config["datatype"] == "SampleData[PairedEndSequencesWithQuality]":

    rule join_ends:
        input:
            "results/{date}/out/trimmed-seqs.qza",
        output:
            "results/{date}/out/joined-seqs.qza",
        params:
            minovlen=config["sequence_joining"]["seq_join_length"],
            minlen=config["sequence_joining"]["minlen"],
            maxdiffs=config["sequence_joining"]["maxdiffs"],
            qmin=config["sequence_joining"]["qmin"],
            qminout=config["sequence_joining"]["qminout"],
            qmax=config["sequence_joining"]["qmax"],
            qmaxout=config["sequence_joining"]["qmaxout"],
            threads=config["sequence_joining"]["threads"],
        log:
            "logs/{date}/preprocessing/join-ends.log",
        conda:
            "../envs/qiime-only-env.yaml"
        shell:
            "qiime vsearch join-pairs "
            "--i-demultiplexed-seqs {input} "
            "--p-minovlen {params.minovlen} "
            "--p-minlen {params.minlen} "
            "--p-maxdiffs {params.maxdiffs} "
            "--p-qmin {params.qmin} "
            "--p-qminout {params.qminout} "
            "--p-qmax {params.qmax} "
            "--p-qmaxout {params.qmaxout} "
            "--p-threads {params.threads} "
            "--o-joined-sequences {output}"
