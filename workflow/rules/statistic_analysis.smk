rule alpha:
    input:
        "results/{date}/out/average-rarefied-table.qza",
    output:
        "results/{date}/out/alpha-diversity-{metric}-normal.qza",
    params:
        metric="{metric}",
    log:
        "logs/{date}/visualisation/alpha-diversity-{metric}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha "
        "--i-table {input} "
        "--p-metric {params.metric} "
        "--o-alpha-diversity {output} "
        "--verbose 2> {log}"


rule alpha_significance:
    input:
        "results/{date}/out/alpha-diversity-{metric}-{diversity}.qza",
    output:
        "results/{date}/visual/alpha-significance-{metric}-{diversity}.qzv",
    log:
        "logs/{date}/visualisation/alpha-significance-{metric}-{diversity}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha-group-significance "
        "--i-alpha-diversity {input} "
        "--m-metadata-file config/pep/sample.tsv "
        "--o-visualization {output} "
        "--verbose 2> {log} "


rule alpha_correlation:
    input:
        "results/{date}/out/alpha-diversity-{metric}-{diversity}.qza",
    output:
        "results/{date}/visual/alpha-correlation-{metric}-{diversity}.qzv",
    params:
        metadata="config/pep/sample.tsv",
        method=config["diversity"]["alpha"]["correlation-method"],
    log:
        "logs/{date}/visualisation/alpha-correlation-{metric}-{diversity}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha-correlation "
        "--i-alpha-diversity {input} "
        "--m-metadata-file {params.metadata} "
        "--p-method {params.method} "
        "--p-intersect-ids "
        "--o-visualization {output} "
        "--verbose 2> {log}"


rule beta:
    input:
        "results/{date}/out/average-rarefied-table.qza",
    output:
        "results/{date}/out/beta-diversity-{metric}-normal.qza",
    params:
        metric="{metric}",
        pseudocount=config["diversity"]["beta"]["diversity-pseudocount"],
        n_jobs=config["diversity"]["beta"]["diversity-n-jobs"],
    log:
        "logs/{date}/visualisation/beta-diversity-{metric}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta "
        "--i-table {input} "
        "--p-metric {params.metric} "
        "--p-pseudocount {params.pseudocount} "
        "--p-n-jobs {params.n_jobs} "
        "--o-distance-matrix {output} "
        "--verbose 2> {log}"


rule beta_significance:
    input:
        "results/{date}/out/beta-diversity-{metric}-{diversity}.qza",
    output:
        out="results/{date}/visual/beta-significance-{metric}-{diversity}-{metadata_column}.qzv",
    params:
        metadata=config["metadata-parameters"]["beta-metadata-column"],
    log:
        "logs/{date}/visualisation/beta-significance-{metric}-{diversity}-{metadata_column}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta-group-significance "
        "--i-distance-matrix {input} "
        "--m-metadata-file config/pep/sample.tsv "
        "--m-metadata-column {wildcards.metadata_column} "
        "--o-visualization {output.out} "
        "--verbose 2> {log}"


rule beta_correlation:
    input:
        "results/{date}/out/beta-diversity-{metric}-{diversity}.qza",
    output:
        distance_matrix="results/{date}/out/beta-correlation-{metric}-{diversity}-{metadata_column}.qza",
        mantel_scatter_vis="results/{date}/visual/beta-correlation-scatter-{metric}-{diversity}-{metadata_column}.qzv",
    params:
        metadata_file="config/pep/sample.tsv",
        method=config["diversity"]["beta"]["correlation-method"],
        permutations=config["diversity"]["beta"]["correlation-permutations"],
        metadata="{metadata_column}",
    log:
        "logs/{date}/visualisation/beta-correlation-{metric}-{diversity}-{metadata_column}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta-correlation "
        "--i-distance-matrix {input} "
        "--m-metadata-file {params.metadata_file} "
        "--m-metadata-column {wildcards.metadata_column} "
        "--p-method {params.method} "
        "--p-permutations {params.permutations} "
        "--p-intersect-ids "
        "--p-label1 jaccard-distance-matrix "
        "--p-label2 {params.metadata} "
        "--o-metadata-distance-matrix {output.distance_matrix} "
        "--o-mantel-scatter-visualization {output.mantel_scatter_vis} "
        "--verbose 2> {log}"


rule alpha_phylogeny:
    input:
        phylogeny="results/{date}/visual/rooted-tree.qza",
        table="results/{date}/out/table-taxa-filtered.qza",
    output:
        "results/{date}/out/alpha-diversity-faith_pd-phylogenetic.qza",
    log:
        "logs/{date}/visualisation/alpha-phylogeny.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity alpha-phylogenetic "
        "--i-phylogeny {input.phylogeny} "
        "--i-table {input.table} "
        "--p-metric faith_pd "
        "--o-alpha-diversity {output} "
        "--verbose 2> {log} "


rule beta_phylogeny:
    input:
        table="results/{date}/out/table-taxa-filtered.qza",
        phylogeny="results/{date}/visual/rooted-tree.qza",
    output:
        "results/{date}/out/beta-diversity-{metric}-phylogenetic.qza",
    params:
        metrics="{metric}",
        threads=config["threads"],
        variance_adjusted=config["diversity"]["beta"]["phylogeny-variance-adjusted"],
    log:
        "logs/{date}/visualisation/beta-phylogeny-{metric}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity beta-phylogenetic "
        "--i-table {input.table} "
        "--i-phylogeny {input.phylogeny} "
        "--p-metric {params.metrics} "
        "--p-threads {params.threads} "
        "--p-variance-adjusted {params.variance_adjusted} "
        "--o-distance-matrix {output} "
        "--verbose 2> {log}"


rule PCoA:
    input:
        "results/{date}/out/beta-diversity-{metric}-{diversity}.qza",
    output:
        "results/{date}/visual/beta-pcoa-{metric}-{diversity}.qza",
    log:
        "logs/{date}/visualisation/beta-pcoa-{metric}-{diversity}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime diversity pcoa "
        "--i-distance-matrix {input} "
        "--o-pcoa {output} "
        "--verbose 2> {log} "


rule emperor:
    input:
        "results/{date}/visual/beta-pcoa-{metric}-{diversity}.qza",
    output:
        "results/{date}/visual/emperor-{metric}-{diversity}.qzv",
    params:
        metadata="config/pep/sample.tsv",
    log:
        "logs/{date}/visualisation/emperor-{metric}-{diversity}.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime emperor plot "
        "--i-pcoa {input} "
        "--m-metadata-file {params.metadata} "
        "--o-visualization {output} "
        "--verbose 2> {log} "
