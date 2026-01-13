#!/bin/bash


THREADS=16
RAW=raw_data # cartella con i file grezzi fastq.gz di input
TRIM=trimmed # cartella output per i file fastq.gz dopo trimming
UMI=umi_processed  # cartella output per i file fastq.gz dopo estrazione UMI
ALIGN=alignment # cartella output per i file sam dopo allineamento
COUNT=counts # cartella output per i file di conteggio
REF=reference/hsa_mature  # path riferimento bowtie

#####################################################
# PIPELINE DI ANALISI
#####################################################

set -euo pipefail

# CREO CARTELLA OUTPUT
mkdir -p fastqc $TRIM $UMI $ALIGN $COUNT logs  # creo una cartella e varie sottocartelle all'interno # "-p" non blocca il comando se la cartella esiste già

# VERIFICA LUNGHEZZA READS, QUALITA', ADATTATORI
fastqc -o fastqc $RAW/*.fastq.gz

# TRIMMING ADATTATORI 
ADAPTER="TGGAATTCTCGGGTGCCAAGG"  # SEQUENZA TIPICA DEGLI ADATTATORI 3' QIAseq

for fq in $RAW/*.fastq.gz; do
    base=$(basename $fq .fastq.gz)
    cutadapt -a $ADAPTER \
        -o $TRIM/${base}.trimmed.fastq.gz \
        $fq > logs/${base}_cutadapt.log
done

# ESTRAZIONE UMI
for fq in $TRIM/*.trimmed.fastq.gz; do
    base=$(basename $fq .trimmed.fastq.gz)
    umi_tools extract \
        --bc-pattern=NNNNNNNNNNNN \
        --stdin=$fq \
        --stdout=$UMI/${base}.umi.fastq.gz \
        --log=logs/${base}_umi_extract.log
done

# ALLINEAMENTO CON BOWTIE
for fq in "$UMI"/*.umi.fastq.gz; do
    base=$(basename "$fq" .umi.fastq.gz)

    bowtie -v 1 -k 1 --best --strata \
        reference/hsa_mature \
        "$fq" \
        "$ALIGN/${base}.sam"
done

# CONTEGGIO CON UMI-TOOLS

for sam in $ALIGN/*.sam; do
    base=$(basename $sam .sam)
    umi_tools count \
        --per-gene \
        --gene-tag=XT \
        --assigned-status-tag=XS \
        -I $sam \
        -S $COUNT/${base}_counts.tsv
done

mkdir -p results  # serve per i risultati di R

################################################
# CARTELLA METADATA.TSV x R

# VERIFICO I FILE "COUNTS" PRESENTI
ls counts/  

# CREO FILE METADATA.TSV
nano metadata.tsv

# ESEMPIO CONTENUTO METADATA.TSV
# Tra "sample" e "condition" devi premere TAB, non spazi
# NO commenti dentro il file 

sample	condition  
S01_sample.counts.tsv	control
S02_sample.counts.tsv	control
S03_sample.counts.tsv	atopic_dermatitis
S04_sample.counts.tsv	atopic_dermatitis
S05_sample.counts.tsv	pemphigoid
S06_sample.counts.tsv	psoriasis

# ogni riga = 1 campione biologico  
# Il valore in "sample" deve corrispondere ESATTAMENTE al nome del file dei conteggi

# CTRL + O # → salva
# INVIO # → conferma nome
# CTRL + X # → esci
