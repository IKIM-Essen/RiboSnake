rule fastqc:
    input:
        "data/{date}/{names}",
    output:
        html="results/{date}/out/fastqc/{names}.html",
        zip="results/{date}/out/fastqc/{names}_fastqc.zip",
    log:
        "logs/{date}/fastqc/{names}.log"
    threads: 8
    wrapper:
        "v1.2.0/bio/fastqc"


rule multiqc:
    input:
        expand(
            "results/{{date}}/out/fastqc/{names}_fastqc.zip",
            names=get_filenames(),
        ),
    output:
        "results/{date}/out/multiqc.html",
    log:
        "logs/{date}/multiqc.log"
    wrapper:
        "0.69.0/bio/multiqc"
