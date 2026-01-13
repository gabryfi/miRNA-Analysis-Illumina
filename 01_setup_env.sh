#!/bin/bash

# SETUP AMBIENTE DI LAVORO
###################################

# AMBIENTE CONDA APPOSITO
conda create -n mirna_qiaseq -y
conda activate mirna_qiaseq

# INSTALLAZIONE DEI SOFTWARE NECESSARI
conda install -y -c bioconda fastqc cutadapt bowtie samtools umi_tools
conda install -y -c conda-forge pigz

R -e 'install.packages("BiocManager")'
R -e 'BiocManager::install(c("edgeR","limma","pheatmap","EnhancedVolcano","clusterProfiler","org.Hs.eg.db"))'

# DOWNLOAD DEI FILE DI RIFERIMENTO miRNA
wget https://mirbase.org/download/hairpin.fa.gz
wget https://mirbase.org/download/mature.fa.gz
gunzip *.gz
grep -A1 "hsa-" mature.fa > hsa_mature.fa
grep -A1 "hsa-" hairpin.fa > hsa_hairpin.fa

# INDICIZZAZIONE DEI FILE DI RIFERIMENTO CON BOWTIE
bowtie-build hsa_mature.fa hsa_mature
