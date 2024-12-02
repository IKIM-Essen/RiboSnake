# RiboSnake: 16S rRNA analysis workflow with QIIME2 and Snakemake

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.10-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/16S.svg?branch=master)](https://travis-ci.org/snakemake-workflows/16S)

Qiime2 workflow for 16S analysis created with snakemake.

## Authors

* Ann-Kathrin Dörr (@AKBrueggemann)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

If you want to use the workflow, please obtain a copy of it by either:
[Cloning](https://help.github.com/en/articles/cloning-a-repository) the repository to your local system, into the place where you want to perform the data analysis or
Downloading a zip-file of the repository to your local machine.

When you have the folder structure added on your local machine, please add a "data" folder manually.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `metadata.txt` to specify your sample setup. If you want to use one of our many preset config-files holding parameters specific for certain sample types or kits, please choose the most fitting for you and rename it into `config.yaml`.

Some important parameters you should check and set according to your own FASTQ-files in the `config.yaml` are primers for the forward and reverse reads, the `datatype`, that should be used by QIIME2 and the `min-seq-length`. Based on the sequencing, the length of the reads can vary.

The default parameters for filtering and truncation were validated with the help of a MOCK community and fitted to retrieve all bacteria from that community.

In addition to that, you need to fit the metadata-parameters to your data. Please change the names of the used metadata-columns according to your information.
Take special care of the "remove-columns" information. Here you can add the columns you don't want to have analyzed or the workflow can't anlyse. This can happen when
all of the values in one column are unique or all the same. You should also look out for the information under "metadata-parameters" and "songbird" as well as "ancom".
In every case you have to specify the column names based on your own data.

If your metadata is not containing numeric values, please use the "reduced-analysis" option in the config file to run the workflow, as the workflow is currently not able to run only on categorical metadata for the full analysis version. We are going to fix that in the future.

The workflow is able to perform clustering and denoising either with vsearch, leading to OTU creation, or with DADA2, creating ASVs. You can decide which modus to use by setting the variable "DADA2" to `True` (DADA2 usage) or `False` (vsearch).

Please make sure, that the names of your FASTQ files are correctly formatted. They should look like this:

    samplename_SNumber_Lane_R1/R2_001.fastq.gz

In the config file you can also set the input and output directory. You can either create a specific directory for your input data and then put that filepath in the config file, or you can put the path to an existing directory where the data is located.
The data will then be copied to the workflow's data directory. The compressed and final file holding the results will be copied to the directory you specified in "output". It will also stay in the local "results" folder together with important intermediate results.
The "data" folder is also not provided by the repository. It is the folder the fastq files are copied to before being used in the workflow. It is best if you create the folder inside the workflows folder structure. It must definitely be created on the machine, the workflow is running on.

### Step 3: Install Snakemake

Create a snakemake environment using [mamba](https://mamba.readthedocs.io/en/latest/) via:

    mamba create -c conda-forge -c bioconda -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Fill up the `metadata.txt` with the information of your samples:

    Please be careful to not include spaces between the commas. If there is a column, that you don't have any information about, please leave it empty and simply go on with the next column.

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Executing the workflow takes two steps:

    Data preparation: snakemake --cores $N --use-conda data_prep
    Workflow execution: snakemake --cores $N --use-conda

using `$N` cores.

When running on snakemake > 8.0 we recommend setting the --shared-fs-usage none as well as setting the environment variable TEMP to a local directory to prevent problems with the usage of the fs-storage system.
The environment variable can be set like this:

    conda activate your_environment_name
    export TEMP=/path/to/local/tmp

Then run the snakemake command like above with the addition of the storage flag:

    snakemake --cores $N --use-conda --shared-fs-usage none

### Step 5: Investigate results

After successful execution, the workflow provides you with a compressed folder, holding all interesting results ready to decompress or to download to your local machine.
The compressed file 16S-report.tar.gz holds several qiime2-artifacts that can be inspected via qiime-view. In the zipped folder report.zip is the snakemake html
report holding graphics as well as the DAG of the executed jobs and html files leading you directly to the qiime2-results, without the need of using qiime-view.

This report can, e.g., be forwarded to your collaborators.

### Step 6: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/16S.git` or `git remote add -f upstream https://github.com/snakemake-workflows/16S.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.

## Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

### Step 1: Forking the repository

[Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.

### Step 2: Cloning

[Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.

### Step 3: Add changes

1. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
2. Commit and push your changes to your fork.
3. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.
4. If you want to add your config file and the parameters as a new default parameter sets, please do this by opening a pull request adding the file to the "contributions" folder.

## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).
If you want to test the RiboSnake functions yourself, you can use the same data used for the CI/CD tests. The used fastq files can be downloaded [here](https://data.qiime2.org/2022.2/tutorials/importing/casava-18-paired-end-demultiplexed.zip). They have been published by Neilson et al., mSystems, 2017.

### Example

1. First clone teh repository to your local machine as described above.
2. Download a dataset of your liking, or the data used for testing the pipeline. The FASTQ files can be downloaded with:
    curl -sL \
          "https://data.qiime2.org/2022.2/tutorials/importing/casava-18-paired-end-demultiplexed.zip"
3. Unzip the data into a folder of your liking, it can be called "incoming" but it does not have to be.
If you name your folder differently, please change the "input" path in the config file.
4. If you don't want to use the whole dataset for testing, remove some of the FASTQ files from the folder:
    rm PAP*
    rm YUN*
    rm Rep*
    rm blank*
5. Use the information that can be found in [this](https://data.qiime2.org/2024.5/tutorials/atacama-soils/sample_metadata.tsv) file from the Qiime2 tutorial, to fill out your metadata.txt file for the samples starting with "BAQ".
6. The default-parameters to be used in the config file can be found in the provided file "PowerSoil-Illumina-soil.yaml" in the config folder.
7. With these parameters and the previous steps, you should be able to execute the workflow.

## Tools

A list of the tools used in this pipeline:

| Tool         | Link                                              |
|--------------|---------------------------------------------------|
| QIIME2       | www.doi.org/10.1038/s41587-019-0209-9             |
| Snakemake    | www.doi.org/10.12688/f1000research.29032.1        |
| FastQC       | www.bioinformatics.babraham.ac.uk/projects/fastqc |
| MultiQC      | www.doi.org/10.1093/bioinformatics/btw354         |
| pandas       | pandas.pydata.org                                 |
| kraken2      | www.doi.org/10.1186/s13059-019-1891-0             |
| vsearch      | www.github.com/torognes/vsearch                   |
| DADA2        | www.doi.org/10.1038/nmeth.3869                    |
| songbird     | www.doi.org/10.1038/s41467-019-10656-5            |
| bowtie2      | www.doi.org/10.1038/nmeth.1923                    |
| Ancom        | www.doi.org/10.3402/mehd.v26.27663                |
| cutadapt     | www.doi.org/10.14806/ej.17.1.200                  |
| BLAST        | www.doi.org/10.1016/S0022-2836(05)80360-2         |
| gneiss       | www.doi.org/10.1128/mSystems.00162-16             |
| qurro        | www.doi.org/10.1093/nargab/lqaa023                |
| Rescript     | www.doi.org/10.1371/journal.pcbi.1009581          |
| EMPeror      | www.doi.org/10.1186/2047-217X-2-16                |
| EMPress      | www.doi.org/10.1128/msystems.01216-20             |

## Citation

If you use RiboSnake in your work, please cite the paper:

Ann-Kathrin Dörr,Josefa Welling,Adrian Dörr,Jule Gosch,Hannah Möhlen,Ricarda Schmithausen,Jan Kehrmann,Folker Meyer,Ivana Kraiselburd,RiboSnake – a user-friendly, robust, reproducible, multipurpose and documentation-extensive pipeline for 16S rRNA gene microbiome analysis,Gigabyte,2024 www.doi.org/10.46471/gigabyte.132
