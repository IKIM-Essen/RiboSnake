rule test:
    output:
        "results/test.txt",
    log:
        "logs/test.log",
    shell:
        "echo 'hello world' > {output} 2> {log}"
