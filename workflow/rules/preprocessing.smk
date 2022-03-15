rule read_samples:
    input:
        tsv = "config/pep/sample.tsv",
        info = "config/pep/sample_info.txt",
    output:
        "results/{date}/out/demux-paired-end.qza"
    params:
        direc = get_data_dir()
    shell:
        "qiime tools import "
            "--type 'SampleData[PairedEndSequencesWithQuality]' "
            "--input-path {params.direc} "
            "--input-format CasavaOneEightSingleLanePerSampleDirFmt "
            "--output-path {output}"

rule trim_paired:
    input:
        "results/{date}/out/demux-paired-end.qza"
    output:
        "results/{date}/out/trimmed-seqs.qza"
    params:
        primer1 = config["primer1"],
        primer2 = config["primer2"],
        min_length = 6
    shell:
        "qiime cutadapt trim-paired "
            "--i-demultiplexed-sequences {input} "
            "--p-adapter-f {params.primer1} "
            "--p-front-f {params.primer2} "
            "--p-minimum-length {params.min_length} "
            "--o-trimmed-sequences {output}"


rule join_ends:
    input:
        "results/{date}/out/trimmed-seqs.qza"
    output:
        "results/{date}/out/joined-seqs.qza"
    params:
        minlen = 30
    shell:
        "qiime vsearch join-pairs "
            "--i-demultiplexed-seqs {input} "
            "--p-minovlen {params.minlen} "
            "--o-joined-sequences {output}"