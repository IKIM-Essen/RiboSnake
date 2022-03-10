rule data_prep:
    input:
        incoming = config["input"],
        data = config["data"],
        metadata = "/local/work/16S/snakemake_qiime/16S/config/pep/metadata.txt",
    output:
        "config/pep/sample.tsv"
    params:
        date = get_date()
    script:
        "../scripts/create_sample_metadata.py"