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


rule get_database:
    output:
        seq="resources/silva-138-99-seqs.qza",
        tax="resources/silva-138-99-tax.qza",
        ref_gen_packed="resources/GRCh38_latest_rna.fna.gz",
    params:
        seq=str(config["database"]["download-path-seq"]),
        tax=str(config["database"]["download-path-tax"]),
        ref_seq=str(config["database"]["ref-genome"])
    log:
        "logs/prep_database.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cd resources; "
        "wget {params.seq}; "
        "wget {params.tax}; "
        "wget {params.ref_seq}; "


rule unzip_ref_gen:
    input:
        "resources/GRCh38_latest_rna.fna.gz"
    output:
        "resources/GRCh38_latest_rna.fna"
    log:
        "logs/unzip_ref_gen.log",
    conda:
        "../envs/python.yaml"
    shell:
        "gzip -d {input}"


rule import_ref_genome:
    input:
        "resources/GRCh38_latest_genomic_upper.fna"
    output:
        "resources/GRCh38_latest_genomic_upper.qza"
    log:
        "logs/import_ref_gen.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime tools import "
        "--input-path {input} "
        "--output-path {output} "
        "--type 'FeatureData[Sequence]' "