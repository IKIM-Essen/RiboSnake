if (
    config["database"]["NCBI"] == True,
    config["database"]["Silva"] == False,
    config["database"]["Greengenes"] == False,
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
