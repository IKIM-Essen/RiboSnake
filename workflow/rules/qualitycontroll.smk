rule fastqc:
    input:
        "data/{date}/{names}",
    output:
        html="results/{date}/out/fastqc/{names}.html",
        zip="results/{date}/out/fastqc/{names}_fastqc.zip",
    log:
        "logs/{date}/fastqc/{names}.log",
    threads: 8
    wrapper:
        "v1.3.1/bio/fastqc"


rule multiqc:
    input:
        expand(
            "results/{{date}}/out/fastqc/{names}_fastqc.zip",
            names=get_filenames(),
        ),
    output:
        "results/{date}/out/multiqc.html",
    log:
        "logs/{date}/multiqc.log",
    wrapper:
        "v1.3.1/bio/multiqc"


rule multiqc_report:
    input:
        expand(
            "results/{{date}}/out/fastqc/{names}_fastqc.zip",
            names=get_filenames(),
        ),
    output:
        report(
            "results/{date}/visual/report/multiqc.html",
            htmlindex="multiqc.htm",
            category="4. Qualitycontrol",
        ),
    log:
        "logs/{date}/multiqc.log",
    wrapper:
        "v1.3.1/bio/multiqc"
