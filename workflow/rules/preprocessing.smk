rule read_samples:
    input:
        tsv="config/pep/sample.tsv",
        info="config/pep/sample_info.txt",
    output:
        "results/{date}/out/demux-paired-end.qza",
    params:
        direc=get_data_dir(),
    log:
        "logs/{date}/preprocessing/read-samples.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools import "
        "--type 'SampleData[PairedEndSequencesWithQuality]' "
        "--input-path {params.direc} "
        "--input-format CasavaOneEightSingleLanePerSampleDirFmt "
        "--output-path {output}"


rule trim_paired:
    input:
        "results/{date}/out/demux-paired-end.qza",
    output:
        "results/{date}/out/trimmed-seqs.qza",
    params:
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
        "qiime cutadapt trim-paired "
        "--i-demultiplexed-sequences {input} "
        "--p-adapter-f {params.primer1} "
        "--p-front-f {params.primer2} "
        "--p-error-rate {params.error_rate} "
        "--p-times {params.rep_times} "
        "--p-overlap {params.overlap} "
        "--p-minimum-length {params.min_length} "
        "--o-trimmed-sequences {output}"


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
