if (
    config["database"]["Silva"] == True,
    config["database"]["Greengenes"] == False,
    config["database"]["NCBI"] == False,
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
