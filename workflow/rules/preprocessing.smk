rule read_samples:
    input:
        tsv="config/pep/sample.tsv",
        info="config/pep/sample_info.txt",
    output:
        temp("results/{date}/out/demux-paired-end.qza"),
    params:
        direc=lambda wc, input: os.path.dirname(pd.read_csv(input.info)["path1"].iloc[0]),
        datatype=config["datatype"],
    log:
        "logs/{date}/preprocessing/read-samples.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        """
        direc=$(cut -d',' -f3 config/pep/sample_info.txt | sed -n '2p' | xargs dirname)
        qiime tools import \
            --type {params.datatype} \
            --input-path "$direc" \
            --input-format CasavaOneEightSingleLanePerSampleDirFmt \
            --output-path {output} \
            2> {log}
        """

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


if (
    config["datatype"] == "SampleData[PairedEndSequencesWithQuality]"
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
