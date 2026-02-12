###################################
# CREAZIONE AMBIENTE CONDA
###################################

conda create -n mirna_qiaseq -y
conda activate mirna_qiaseq

###################################
# INSTALLAZIONE SOFTWARE
###################################

conda install -y -c bioconda \
  fastqc \
  cutadapt \
  bowtie \
  samtools \
  umi_tools \
  subread

conda install -y -c conda-forge pigz

###################################
# R + PACCHETTI BIOC
###################################

R -e 'install.packages("BiocManager")'
R -e 'BiocManager::install(c("edgeR","limma","pheatmap","EnhancedVolcano","UpSetR","clusterProfiler","org.Hs.eg.db"))'

###################################
# DOWNLOAD miRBase
###################################

# reference miRNA (mirBASE)

wget https://mirbase.org/download/mature.fa
wget https://mirbase.org/download/hairpins.fa

# https://www.mirbase.org/download/ #

cat mature.fa hairpins.fa > hsa_mature.fa          # concateno i 2 file (miRNA maturi e hairpins)

# solo Homo sapiens
grep -A1 "^>hsa-" hsa_mature.fa > hsa_miRNA.fa

###################################
# INDICIZZAZIONE BOWTIE
###################################

bowtie-build hsa_miRNA.fa hsa_miRNA

###################################
###################################
# VARIABILI
###################################

THREADS=16

RAW=raw_data
TRIM=trimmed
ALIGN=alignment
COUNT=counts
LOGS=logs

REF_IDX=/home/user/reference/hsa_miRNA

###################################
##### ADAPTER QIAseq miRNA #####
##### VERIFICA SUL MANUALE #####
##### E MODIFICA QUI !!!!! #####
###################################

ADAPTER="INSERISCI_QUI_ADAPTER_QIASEQ_CORRETTO" ################################################################################

mkdir -p fastqc $TRIM $ALIGN $COUNT $LOGS

fastqc -o fastqc -t $THREADS $RAW/*.fastq.gz
