if (
    config["database"]["Greengenes"] == True,
    config["database"]["Silva"] == False,
    config["database"]["NCBI"] == False,
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
