#!/bin/bash
# Bash script for qiime2

INPUTDIR=/local/work/16S/qiime2/16SHands/Optional
OUTDIR=/local/work/16S/qiime2/single_fastqs/results/2022-02-17-silva/out
VISUALISE=/local/work/16S/qiime2/single_fastqs/results/2022-02-17-silva/visual
PRIMERSEQ1=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
PRIMERSEQ2=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
start=`date +%s`

# Importing the sequences and creating a qiime artifact
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $INPUTDIR/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $OUTDIR/demux-paired-end.qza

qiime demux summarize \
  --i-data $OUTDIR/demux-paired-end.qza  \
  --o-visualization $VISUALISE/demux-seqs.qzv

#Cutting the primers from the sequences
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $OUTDIR/demux-paired-end.qza \
  --p-adapter-f $PRIMERSEQ1 \
  --p-front-f $PRIMERSEQ2 \
  --p-minimum-length 6 \
  --o-trimmed-sequences $OUTDIR/trimmed_seqs.qza
# Join paired ends 

qiime vsearch join-pairs \
  --i-demultiplexed-seqs $OUTDIR/trimmed_seqs.qza \
  --p-minovlen 30 \
  --o-joined-sequences $OUTDIR/demux-joined.qza

#Filtering the sequences for quality
qiime quality-filter q-score \
  --i-demux $OUTDIR/demux-joined.qza \
  --p-min-quality 20 \
  --p-min-length-fraction 0.75 \
  --p-max-ambiguous 50 \
  --o-filtered-sequences $OUTDIR/demux-joined-filtered.qza \
  --o-filter-stats $OUTDIR/demux-joined-filter-stats.qza

qiime vsearch fastq-stats \
  --i-sequences $OUTDIR/demux-paired-end.qza \
  --o-visualization $VISUALISE/fastq_stats.qzv
#FastQC and MultiQC for basic sample information 
fastqc $INPUTDIR/*.fastq.gz -o $OUTDIR/fastQC
#cd ${WORKDIR}/fastqc_analysis
#multiqc $OUTDIR/fastQC/fastqc_analysis

#Dereplicating sequences
qiime vsearch dereplicate-sequences \
  --i-sequences $OUTDIR/demux-joined-filtered.qza \
  --o-dereplicated-table $OUTDIR/table.qza \
  --o-dereplicated-sequences $OUTDIR/rep-seqs.qza

#de-novo clustering sequences
qiime vsearch cluster-features-de-novo \
  --i-table $OUTDIR/table.qza \
  --i-sequences $OUTDIR/rep-seqs.qza \
  --p-perc-identity 0.99 \
  --p-threads 4 \
  --o-clustered-table $OUTDIR/table-dn-99.qza \
  --o-clustered-sequences $OUTDIR/rep-seqs-dn-99.qza

# creating a feature table to visualise
qiime feature-table tabulate-seqs \
  --i-data $OUTDIR/rep-seqs-dn-99.qza \
  --o-visualization $VISUALISE/rep-seqs-dn-99.qzv

#Finding chimeras
qiime vsearch uchime-denovo \
  --i-table $OUTDIR/table-dn-99.qza \
  --i-sequences $OUTDIR/rep-seqs-dn-99.qza \
  --output-dir $OUTDIR/uchime-dn-out

# Removing chimeras from feature table and rep-seqs
qiime feature-table filter-features \
  --i-table $OUTDIR/table-dn-99.qza \
  --m-metadata-file $OUTDIR/uchime-dn-out/nonchimeras.qza \
  --o-filtered-table $OUTDIR/uchime-dn-out/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
  --i-data $OUTDIR/rep-seqs-dn-99.qza \
  --m-metadata-file $OUTDIR/uchime-dn-out/nonchimeras.qza \
  --o-filtered-data $OUTDIR/uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza

#visualising the output of the chimera removal
qiime metadata tabulate \
  --m-input-file $OUTDIR/uchime-dn-out/stats.qza \
  --o-visualization $VISUALISE/uchime_dn_out_stats.qzv

qiime feature-table tabulate-seqs \
  --i-data $OUTDIR/uchime-dn-out/nonchimeras.qza \
  --o-visualization $VISUALISE/nonchimeras.qzv

qiime feature-table tabulate-seqs \
  --i-data $OUTDIR/uchime-dn-out/chimeras.qza \
  --o-visualization $VISUALISE/chimeras.qzv

# Abundance filtering of the data (if data from some sequences is very rare) -> do we want that?
# more filtering for the number/proportion of one found feature in one sample
qiime feature-table filter-features \
  --i-table $OUTDIR/uchime-dn-out/table-nonchimeric-wo-borderline.qza \
  --p-min-frequency 20 \
  --o-filtered-table $OUTDIR/uchime-dn-out/table-nonchimeric-filtered.qza

qiime feature-table filter-seqs \
  --i-data $OUTDIR/uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza \
  --i-table $OUTDIR/uchime-dn-out/table-nonchimeric-filtered.qza \
  --p-no-exclude-ids \
  --o-filtered-data $OUTDIR/uchime-dn-out/rep-seqs-nonchimeric-filtered.qza

qiime feature-table summarize \
  --i-table $OUTDIR/uchime-dn-out/table-nonchimeric-filtered.qza \
  --o-visualization $VISUALISE/table.qzv

qiime feature-table tabulate-seqs \
  --i-data $OUTDIR/uchime-dn-out/rep-seqs-nonchimeric-filtered.qza \
  --o-visualization $VISUALISE/rep-seqs.qzv

qiime feature-table filter-features-conditionally \
  --i-table $OUTDIR/uchime-dn-out/table-nonchimeric-filtered.qza \
  --p-abundance 0.015 \
  --p-prevalence 0.001 \
  --o-filtered-table $OUTDIR/table_abundance_filter.qza

qiime feature-table filter-seqs \
  --i-data $OUTDIR/uchime-dn-out/rep-seqs-nonchimeric-filtered.qza \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --p-no-exclude-ids \
  --o-filtered-data $OUTDIR/rep-seqs-abundance-filtered.qza

# Classifying the data via vsearch
# importing the fasta files of the refrenece database
"""
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /local/work/16S/qiime2/single_fastqs/Databases/SILVA_138_SSURef_tax_silva.fasta \
  --output-path $OUTDIR/99_otus.qza

# importing the reference taxonomy 
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /local/work/16S/qiime2/single_fastqs/Databases/taxmap_slv_ssu_ref_138.txt \
  --output-path $OUTDIR/ref-taxonomy.qza
"""
qiime feature-classifier classify-consensus-vsearch \
  --i-query $OUTDIR/rep-seqs-abundance-filtered.qza \
  --i-reference-reads /local/work/16S/qiime2/single_fastqs/Databases/silva-138-99-seqs.qza \
  --i-reference-taxonomy /local/work/16S/qiime2/single_fastqs/Databases/silva-138-99-tax.qza \
  --p-maxaccepts 1 \
  --p-maxrejects 1 \
  --p-perc-identity 0.97 \
  --p-threads 4 \
  --o-classification $OUTDIR/class_table.qza \
  --verbose

qiime taxa collapse \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --i-taxonomy $OUTDIR/class_table.qza \
  --p-level 6 \
  --o-collapsed-table $OUTDIR/collapsed_taxa.qza

qiime feature-table heatmap \
  --i-table $OUTDIR/collapsed_taxa.qza \
  --m-sample-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --m-sample-metadata-column swab-site \
  --p-cluster "features" \
  --o-visualization $VISUALISE/heatmap.qzv

qiime metadata tabulate \
  --m-input-file $OUTDIR/collapsed_taxa.qza \
  --o-visualization $VISUALISE/taxonomy.qzv

qiime taxa barplot \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --i-taxonomy $OUTDIR/class_table.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --o-visualization $VISUALISE/taxa-bar-plots.qzv

classification=`date +%s`
echo Half execution time was `expr $classification - $start` seconds.


#Generate tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $OUTDIR/rep-seqs-abundance-filtered.qza \
  --o-alignment $OUTDIR/aligned-rep-seqs.qza \
  --o-masked-alignment $OUTDIR/masked-aligned-rep-seqs.qza \
  --o-tree $VISUALISE/unrooted-tree.qza \
  --o-rooted-tree $VISUALISE/rooted-tree.qza

#Exporting in newick-format
qiime tools export \
  --input-path $VISUALISE/rooted-tree.qza \
  --output-path $VISUALISE

#Alpha and beta diversity
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $VISUALISE/rooted-tree.qza \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --p-sampling-depth 100 \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --output-dir $OUTDIR/core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity $OUTDIR/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --o-visualization $VISUALISE/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $OUTDIR/core-metrics-results/evenness_vector.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --o-visualization $VISUALISE/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $OUTDIR/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --m-metadata-column swab-site \
  --o-visualization $VISUALISE/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

qiime emperor plot \
  --i-pcoa $OUTDIR/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --p-custom-axes year \
  --o-visualization $VISUALISE/unweighted-unifrac-emperor-days-since-experiment-start.qzv

qiime emperor plot \
  --i-pcoa $OUTDIR/core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --p-custom-axes year \
  --o-visualization $VISUALISE/bray-curtis-emperor-days-since-experiment-start.qzv

#Alpha and beta rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --i-phylogeny $VISUALISE/rooted-tree.qza \
  --p-max-depth 500 \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --o-visualization $VISUALISE/alpha-rarefaction.qzv

qiime diversity beta-rarefaction \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --i-phylogeny $VISUALISE/rooted-tree.qza \
  --p-metric euclidean \
  --p-clustering-method nj \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --p-sampling-depth 100 \
  --o-visualization $VISUALISE/beta-rarefaction.qzv

# Abundance filtering with gneiss
qiime gneiss correlation-clustering \
  --i-table $OUTDIR/collapsed_taxa.qza \
  --o-clustering $OUTDIR/hirarchy2.qza

qiime gneiss dendrogram-heatmap \
  --i-table $OUTDIR/collapsed_taxa.qza \
  --i-tree $OUTDIR/hirarchy2.qza \
  --m-metadata-file /local/work/16S/qiime2/16SHands/sample_metadata.tsv \
  --m-metadata-column subject \
  --p-color-map seismic \
  --o-visualization $VISUALISE/heatmap_gneiss2.qzv

# exporting feature table as biom file and taxonomy as .tsv
qiime tools export \
  --input-path $OUTDIR/table_abundance_filter.qza \
  --output-path $OUTDIR

qiime tools export \
  --input-path $OUTDIR/class_table.qza \
  --output-path $OUTDIR

qiime feature-table presence-absence \
  --i-table $OUTDIR/table_abundance_filter.qza \
  --o-presence-absence-table $OUTDIR/table_binary.qza

qiime tools export \
  --input-path $OUTDIR/table_binary.qza \
  --output-path $OUTDIR

end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
