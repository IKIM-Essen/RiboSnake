rule relative_collapsed_taxa:
    input:
        "results/{date}/out/taxa_collapsed.qza",
    output:
        "results/{date}/out/taxa_collapsed_relative.qza",
    log:
        "logs/{date}/outputs/taxa-collapse-relative.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime feature-table relative-frequency "
        "--i-table {input} "
        "--o-relative-frequency-table {output} "
        "--verbose 2> {log} "


rule volatility:
    input:
        table="results/{date}/out/taxa_collapsed_relative.qza",
        alpha=expand(
            "results/{{date}}/out/alpha-diversity-{metric}-normal.qza",
            metric=get_metric("alpha"),
        ),
        alpha_phylo=expand(
            "results/{{date}}/out/alpha-diversity-{metric}-phylogenetic.qza",
            metric=get_phylogenetic_metric("alpha"),
        ),
    output:
        "results/{date}/visual/volatility.qzv",
    params:
        metadata="config/pep/sample.tsv",
        state_column=config["longitudinal-params"]["state_column"],
        individual_id_column=config["longitudinal-params"]["individual_id_column"],
    log:
        "logs/{date}/visualisation/volatility.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime longitudinal volatility "
        "--i-table {input.table} "
        "--m-metadata-file {params.metadata} "
        "--m-metadata-file {input.alpha} "
        "--m-metadata-file {input.alpha_phylo} "
        "--p-state-column {params.state_column} "
        "--p-individual-id-column {params.individual_id_column} "
        "--o-visualization {output} "
        "--verbose 2> {log} "


rule feature_volatility:
    input:
        table="results/{date}/out/taxa_collapsed.qza",
        alpha=expand(
            "results/{{date}}/out/alpha-diversity-{metric}-normal.qza",
            metric=get_metric("alpha"),
        ),
        alpha_phylo=expand(
            "results/{{date}}/out/alpha-diversity-{metric}-phylogenetic.qza",
            metric=get_phylogenetic_metric("alpha"),
        ),
    output:
        table="results/{date}/visual/important_long.qza",
        feature_importance="results/{date}/visual/importance.qza",
        volatility_plot="results/{date}/visual/feature.qzv",
        accuracy_results="results/{date}/visual/accuracy.qzv",
        sample_estimator="results/{date}/visual/estimator.qza",
    params:
        metadata="config/pep/sample.tsv",
        state_column=config["longitudinal-params"]["state_column"],
        individual_id_column=config["longitudinal-params"]["individual_id_column"],
        p_n_jobs=config["threads"],
    log:
        "logs/{date}/visualisation/feature_volatility.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime longitudinal feature-volatility "
        "--i-table {input.table} "
        "--m-metadata-file {params.metadata} "
        "--m-metadata-file {input.alpha} "
        "--m-metadata-file {input.alpha_phylo} "
        "--p-state-column {params.state_column} "
        "--p-individual-id-column {params.individual_id_column} "
        "--p-n-jobs 5 "
        "--o-filtered-table {output.table} "
        "--o-feature-importance {output.feature_importance} "
        "--o-volatility-plot {output.volatility_plot} "
        "--o-accuracy-results {output.accuracy_results} "
        "--o-sample-estimator {output.sample_estimator} "
        "--verbose 2> {log} "


rule linear_mixed_effects:
    input:
        table="results/{date}/out/taxa_collapsed_relative.qza",
        alpha=expand(
            "results/{{date}}/out/alpha-diversity-{metric}-normal.qza",
            metric=get_metric("alpha"),
        ),
        alpha_phylo=expand(
            "results/{{date}}/out/alpha-diversity-{metric}-phylogenetic.qza",
            metric=get_phylogenetic_metric("alpha"),
        ),
    output:
        "results/{date}/visual/lme.qzv",
    params:
        metadata="config/pep/sample.tsv",
        state_column=config["longitudinal-params"]["state_column"],
        individual_id_column=config["longitudinal-params"]["individual_id_column"],
        #formula=config["longitudinal"]["formula"],
        metric=config["longitudinal-params"]["metric"],
    log:
        "logs/{date}/visualisation/linearmixedeff.log",
    conda:
        "../envs/qiime-only-env.yaml"
    shell:
        "qiime longitudinal linear-mixed-effects "
        "--i-table {input.table} "
        "--m-metadata-file {params.metadata} "
        "--m-metadata-file {input.alpha} "
        "--m-metadata-file {input.alpha_phylo} "
        "--p-state-column {params.state_column} "
        "--p-individual-id-column {params.individual_id_column} "
        "--p-metric {params.metric} "
        "--o-visualization {output} "
        "--verbose 2> {log} "


rule unzip_longitudinal:
    input:
        #"results/{date}/visual/lme.qzv",
        "results/{date}/visual/volatility.qzv",
        "results/{date}/visual/feature.qzv",
        "results/{date}/visual/accuracy.qzv",
    output:
        directory("results/{date}/visual/longitudinal_unzipped"),
    log:
        "logs/{date}/outputs/unzip-longitudinal.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rename_qzv.py"


rule report_longitudinal:
    input:
        "results/{date}/visual/longitudinal_unzipped/",
    output:
        feature=report(
            directory("results/{date}/visual/report/feature"),
            caption="../report/feature-volatility.rst",
            category="3. Analysis",
            subcategory="Longitudinal",
            htmlindex="index.html",
        ),
        accuracy=report(
            directory("results/{date}/visual/report/accuracy"),
            caption="../report/feature-accuracy.rst",
            category="3. Analysis",
            subcategory="Longitudinal",
            htmlindex="index.html",
        ),
        general=report(
            directory("results/{date}/visual/report/volatility"),
            caption="../report/volatility.rst",
            category="3. Analysis",
            subcategory="Longitudinal",
            htmlindex="index.html",
        ),
        #lme=report(
        #    "results/{date}/visual/report/lme",
        #    caption="../report/lme.rst",
        #    category="3. Analysis",
        #    subcategory="Longitudinal",
        #    htmlindex="index.html",
        #),
    log:
        "logs/{date}/outputs/report-longitudinal.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_longitudinal.py"
