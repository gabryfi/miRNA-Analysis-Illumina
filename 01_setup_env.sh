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

wget https://mirbase.org/download/mature.fa.gz

# https://www.mirbase.org/download/ #

wget https://mirgenedb.org/fasta/ALL?mat=1

# https://mirgenedb.org/download

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz

####################################

# solo Homo sapiens
grep -A1 "^>hsa-" hsa_mature.fa > hsa_miRNA.fa

###################################
# INDICIZZAZIONE BOWTIE
###################################

bowtie-build hsa_miRNA.fa hsa_miRNA

bowtie-build hsa_mature.fas hsa_miRNA_MirGeneDB

bowtie-build GRCh38_latest_genomic.fna hsa_miRNA_NCBI_GRCh38

###################################
###################################
# VARIABILI
###################################

THREADS=30

RAW=raw_data
TRIM=trimmed
ALIGN=alignment
COUNT=counts
LOGS=logs

REF_IDX=/home/user/reference/hsa_miRNA

REF_IDX=/home/user/reference/hsa_miRNA_MirGeneDB

REF_IDX=/home/user/reference/hsa_miRNA_NCBI_GRCh38

###################################
##### ADAPTER QIAseq miRNA #####
##### VERIFICA SUL MANUALE #####
##### E MODIFICA QUI !!!!! #####
###################################

ADAPTER="INSERISCI_QUI_ADAPTER_QIASEQ_CORRETTO" ################################################################################

mkdir -p fastqc $TRIM $ALIGN $COUNT $LOGS

fastqc -o fastqc -t $THREADS $RAW/*.fastq.gz
