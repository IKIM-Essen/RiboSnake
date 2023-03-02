# General configuration settings

To configure the workflow to your needs, change the settings in config/config.yaml. The default parameters were determined by using test data holding a MOCK community as positive
control. The parameters were set to values that allowed to retrieve all bacterial genera in the MOCK community.

## Metadata sheets

You need to fill out the <em>metadata.txt</em> according to your needs. It holds all numeric and categorical metadata information of your data. Please be careful while filling out
this file, and make sure you don't miss a tab to separate the columns and look for spelling mistakes. Those can lead to problems in the analysis further down and are quite annoying
to look for.

The two files <em>sample_info.txt</em> and <em>sample.tsv</em>, are filled out automatically, if the variable `include-data-prep` is set to <em>true</em> in the config file and you
are running the command `snakemake --cores $N --use-conda data_prep`.
Please set this variable to <em>false</em> when the data_prep step is done.