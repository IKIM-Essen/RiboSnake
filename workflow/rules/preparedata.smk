rule data_prep:
    input:
        config["metadata"],
    output:
        sample_tsv="config/pep/sample.tsv",
        sample_info="config/pep/sample_info.txt",
    priority: 50
    log:
        "logs/data_prep.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_sample_metadata.py"
