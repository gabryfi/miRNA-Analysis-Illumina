###################################
# TRIMMING ADAPTER
###################################

for fq in $RAW/*.fastq.gz; do
  base=$(basename "$fq" .fastq.gz)

  cutadapt \
    -a $ADAPTER \
    -m 18 -M 30 \
    --trim-n \
    -j $THREADS \
    -o $TRIM/${base}.trimmed.fastq.gz \
    "$fq" > $LOGS/${base}_cutadapt.log
done

###################################
# ALIGNMENT miRNA
###################################

for fq in $TRIM/*.trimmed.fastq.gz; do
  base=$(basename "$fq" .trimmed.fastq.gz)

  bowtie \
    -v 2 \
    -k 10 \
    --best \
    --strata \
    -p $THREADS \
    $REF_IDX \
    "$fq" \
    "$ALIGN/${base}.sam"
done

###################################
# CONVERSIONE BAM
###################################

for sam in $ALIGN/*.sam; do
  base=$(basename "$sam" .sam)

  samtools view -bS "$sam" |
  samtools sort -@ $THREADS -o "$ALIGN/${base}.sorted.bam"

  samtools index "$ALIGN/${base}.sorted.bam"
  rm "$sam"
done

###################################
# UMI DEDUP + COUNT
###################################

for bam in $ALIGN/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)

  umi_tools count \
    --per-gene \
    --gene-tag=XT \
    --assigned-status-tag=XS \
    -I "$bam" \
    -S "$COUNT/${base}_counts.tsv" \
    --log="$LOGS/${base}_umi_count.log"
done
