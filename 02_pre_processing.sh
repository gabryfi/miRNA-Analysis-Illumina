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
    -v 2 \                        # Permette fino a 2 mismatch nella mappatura
    -k 10 \                       # Riporta fino a 10 posizioni per reads multimappanti
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
    -v 1 \                        # Permette solo 1 mismatch per read, più stringente
    -m 20 \                        # Riporta fino a 20 posizioni se read multimappante
    --best \                        # Restituisce l’allineamento con miglior punteggio
    --strata \                       # Assicura che solo l’allineamento migliore sia restituito
    -p $THREADS \                    # Parallelizza il mapping
    $REF_IDX \                       # Indice Bowtie creato dal reference genome + miRNA/ncRNA combinati
    "$fq" \                          # Input FASTQ trim-mato
    --sam > "$ALIGN/${base}_article.sam"  # Output SAM separato per questa pipeline
done




##########
# OPTIONAL: separazione reads mapped / unmapped
#########

for sam in $ALIGN/*_article.sam; do
  base=$(basename "$sam" _article.sam)

  # Estrae le reads mappate
  bowtie --sam --all "$sam" > "$ALIGN/${base}_mapped.fastq"

  # Estrae le reads non mappate
  bowtie --sam --un "$sam" > "$ALIGN/${base}_unmapped.fastq"
done

# http://hgdownload.cse.ucsc.edu/downloads.html

###################################
# CONVERSIONE BAM
###################################

for sam in $ALIGN/*.sam; do
  # Estrai il nome base del file senza percorso e estensione
  base=$(basename "$sam" .sam)

  # Converte SAM in BAM binario e ordina per coordinate
  samtools view -bS "$sam" | \
  samtools sort -@ $THREADS -o "$ALIGN/${base}.sorted.bam"

  # Crea indice .bai per accesso veloce al BAM
  samtools index "$ALIGN/${base}.sorted.bam"

done

###################################
# UMI DEDUPLICATION + GENE COUNT
###################################

for bam in $ALIGN/*.sorted.bam; do
  # Estrai il nome base del file BAM senza estensione
  base=$(basename "$bam" .sorted.bam)

  umi_tools count \
    --per-gene \                     # Conta le reads raggruppandole per gene/miRNA
    --gene-tag=XT \                   # Indica il tag BAM dove è annotato il gene/miRNA (XT è il tag usato in bowtie o GTF)
    --assigned-status-tag=XS \        # Usa il tag XS per determinare se la read è stata assegnata al gene correttamente
    -I "$bam" \                       # Input BAM ordinato dal passo precedente
    -S "$COUNT/${base}_counts.tsv" \  # Output: tabella dei conteggi per ogni gene/miRNA
    --log="$LOGS/${base}_umi_count.log"  # Log dettagliato del processo
done

