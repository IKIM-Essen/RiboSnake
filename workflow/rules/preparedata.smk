rule data_prep:
    input:
        config["metadata"],
    output:
        metadata="config/pep/sample.tsv",
        sample_info="config/pep/sample_info.txt",
    params:
        datatype=config["datatype"],
    priority: 50
    log:
        "logs/data_prep.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_sample_metadata.py"
