###################################
# TRIMMING CUTADAPT
###################################

for fq in $RAW/*.fastq.gz; do
  # Estrai il nome base del file senza percorso e estensione
  base=$(basename "$fq" .fastq.gz)

  cutadapt \
    -a $ADAPTER \               # Rimuove l'adattatore 3' con tolleranza mismatches del 10%
    -e 0.1 \                     # Error rate massimo = 10% (come negli articoli QIAseq)
    -m 17 -M 30 \                # Mantiene solo reads con lunghezza tra 17 e 30 nt
    --trim-n \                   # Rimuove eventuali basi N alle estremità delle reads
    -j $THREADS \                # Numero di thread per parallelizzazione
    -o $TRIM/${base}.trimmed.fastq.gz \  # Output trimmed
    "$fq" > $LOGS/${base}_cutadapt.log   # Salva log dettagliato
done

###################################
# ALIGNMENT miRNA
###################################

##########
# ALIGNMENT miRNA - pipeline QIAseq standard
##########

for fq in $TRIM/*.trimmed.fastq.gz; do
  # Estrai il nome base del file senza percorso e estensione
  base=$(basename "$fq" .trimmed.fastq.gz)

  bowtie \
    -v 1 \                         # Permette fino a 1 mismatch nella mappatura
    -k 1 \                         # Riporta fino a n posizioni per reads multimappanti (1 reads= + loci) !!! attenzione DUPLICATI (1 read conta più volte)
    --best \                       # Restituisce solo l'allineamento migliore
    --strata \                     # Seleziona l'allineamento migliore in termini di score
    -p $THREADS \                  # Numero di thread per parallelizzare l'analisi
    $REF_IDX \                     # Indice bowtie del reference miRNA (es. miRBase)
    "$fq" \                        # Input FASTQ trimmed dal passo precedente
    "$ALIGN/${base}.sam"           # Output SAM (poi convertito in BAM con samtools)
done


###################################
# ALIGNMENT miRNA - pipeline stile articolo UCSC/Ensembl
###################################

for fq in $TRIM/*.trimmed.fastq.gz; do
  # Estrai il nome base del file senza percorso e estensione
  base=$(basename "$fq" .trimmed.fastq.gz)

  bowtie \
    -v 1 \                           # Permette solo 1 mismatch per read, più stringente
    -m 1 \                           # Riporta fino a n posizioni se read multimappante
    --best \                         # Restituisce l’allineamento con miglior punteggio
    --strata \                       # Assicura che solo l’allineamento migliore sia restituito
    -p $THREADS \                    # Parallelizza il mapping
    $REF_IDX \                       # Indice Bowtie creato dal reference genome + miRNA/ncRNA combinati
    "$fq" \                          # Input FASTQ trim-mato
    --sam > "$ALIGN/${base}_article.sam"  # Output SAM separato per questa pipeline
done

