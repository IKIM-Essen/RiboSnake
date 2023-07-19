rule create_bowtie_db:
    input:
        "resources/GRCh38_latest_genomic_upper.fna"
    output:
        "resources/bowtie_host_DB"
    log:
        "logs/bowtie_db.log"
    conda:
        "../envs/python.yaml"
    shell:
        "bowtie2-build {input} {output}"

rule map_sequences:
    input:
        db="resources/bowtie_host_DB",
        read1="data/{date}/{sample}_L001_R1_001.fastq.gz",
        read2="data/{date}/{sample}_L001_R2_001.fastq.gz",
    output:
        "results/{date}/out/{sample}_mapped_and_unmapped.sam"
    log:
        "logs/{date}/bowtie/{sample}_mapping.log"
    conda:
        "../envs/python.yaml"
    shell:
        "bowtie2 -p 8 -x {input.db} "
        "-1 {input.read1} "
        "-2 {input.read2} "
        "-S {output}"

rule sam_to_bam:
    input:
        "results/{date}/out/{sample}_mapped_and_unmapped.sam"
    output:
        "results/{date}/out/{sample}_mapped_and_unmapped.bam"
    log:
        "logs/{date}/bowtie/{sample}_sam_to_bam.log"
    conda:
        "../envs/python.yaml"
    shell:
        "samtools view -bS {input} > {output}"

rule filter_unmapped:
    input:
        "results/{date}/out/{sample}_mapped_and_unmapped.bam"
    output:
        "results/{date}/out/{sample}_botReadsUnmapped.bam"
    log:
        "logs/{date}/bowtie/{sample}_filter_unmapped.log"
    conda:
        "../envs/python.yaml"
    shell:
        "samtools view -b -f 12 -F 256 "
        "{input} > {output}"

rule split_paired:
    input:
        "results/{date}/out/{sample}_botReadsUnmapped.bam"
    output:
        expand{
            read1="results/{{date}}/bowtie/{sample}_L001_R1_001.fastq.gz",
            read2="results/{{date}}/bowtie/{sample}_L001_R2_001.fastq.gz",
            sample=get_reads_for_kraken(),
        },
    params:
        "sorted=SAMPLE_bothReadsUnmapped_sorted.bam"
    log:
        "logs/{date}/bowtie/{sample}_split_paired.log"
    conda:
        "../envs/python.yaml"
    shell:
        "samtools sort -n -m 5G -@ 2 {input} -o {params.sorted}"
        "samtools fastq -@ 8 {params.sorted} "
        "-1 {output.read1} "
        "-2 {output.read2} "
        "-0 /dev/null -s /dev/null -n"
