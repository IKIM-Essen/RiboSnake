rule data_prep:
    input:
        config["metadata"],
    output:
        metadata="config/pep/sample.tsv",
        sample_info="config/pep/sample_info.txt",
    priority: 50
    log:
        "logs/data_prep.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_sample_metadata.py"


rule get_database:
    output:
        seq="resources/silva-138-99-seqs.qza",
        tax="resources/silva-138-99-tax.qza",
    params:
        seq=str(config["database"]["download-path-seq"]),
        tax=str(config["database"]["download-path-tax"]),
    log:
        "logs/prep_database.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cd resources; "
        "wget {params.seq}; "
        "wget {params.tax}; "
