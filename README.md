# miRNA-Analysis-Illumina

This script is designed for miRNA analysis in human skin. We have four conditions: healthy, atopic dermatitis, psoriasis, and bullous pephigus. Sequencing is performed on the Illumina NextSeq2000 platform, flowcell: P1.
The script is built using for loops in the Bash language section, while the R section attempts to follow traditional DEG analysis. This analysis focuses solely on miRNAs. The annotation reference file was downloaded from mirBase (https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.mirbase.org/&ved=2ahUKEwjC8fim3dOSAxULhv0HHQjpEJcQFnoECA0QAQ&usg=AOvVaw09u2bgiBrAGVkaj_JYqV03). Both the mature miRNA file (mature.fa) and the hairpin file (haiprpins.fa) were downloaded. These two files are combined and filtered for those belonging to Homo sapiens.
The variables used later by the for loops are set:

THREADS=16
RAW=raw_data
TRIM=trimmed
ALIGN=alignment
COUNT=counts
LOGS=logs
REF_IDX=/home/user/reference/hsa_miRNA
ADAPTER="INSERT_CORRECT_ADAPTER_QIASEQ_HERE"

mkdir -p fastqc $TRIM $ALIGN $COUNT $LOGS

fastqc -o fastqc -t $THREADS $RAW/*.fastq.gz

These must be checked before running the script, in particular: REF_IDX and ADAPTER.

Parametri:

Passo/Strumento	Parametro	Valore	Note
cutadapt	Adattatore 3'	$ADAPTER	Sequenza adattatore da rimuovere
cutadapt	Lunghezza minima read	18	Conserva solo reads ≥18 nt
cutadapt	Lunghezza massima read	30	Conserva solo reads ≤30 nt
cutadapt	Trim basi N	sì	Rimuove eventuali N alle estremità
cutadapt	Thread	$THREADS	Parallelizzazione
cutadapt	Output	$TRIM/${base}.trimmed.fastq.gz	File fastq trim-mato
cutadapt	Log	$LOGS/${base}_cutadapt.log	Log dettagliato
bowtie (pipeline standard)	Mismatches	2	Permette fino a 2 mismatch
bowtie (pipeline standard)	Max multimapping	10	Riporta fino a 10 posizioni
bowtie (pipeline standard)	--best	sì	Restituisce migliore allineamento
bowtie (pipeline standard)	--strata	sì	Considera solo il livello migliore
bowtie (pipeline standard)	Threads	$THREADS	Parallelizzazione
bowtie (pipeline standard)	Indice	$REF_IDX	Indice reference genome / miRNA
bowtie (pipeline standard)	Output	$ALIGN/${base}.sam	File SAM
samtools view / sort / index	Input SAM	$ALIGN/*.sam	Input da bowtie
samtools	Conversione in BAM	-bS	Converte SAM in BAM binario
samtools	Sort	-@ $THREADS	Ordina BAM per coordinate
samtools	Output BAM	$ALIGN/${base}.sorted.bam	BAM ordinato
samtools	Index	sì	Crea file .bai
umi_tools count	Input BAM	$ALIGN/*.sorted.bam	BAM ordinato deduplicato
umi_tools	--per-gene	sì	Conta reads raggruppando per gene/miRNA
umi_tools	--gene-tag	XT	Tag BAM contenente annotazione gene
umi_tools	--assigned-status-tag	XS	Filtra solo reads assegnate
umi_tools	Output counts	$COUNT/${base}_counts.tsv	Tabella counts per miRNA
umi_tools	Log	$LOGS/${base}_umi_count.log	Log dettagliato
R - Lettura counts	Path	counts/	Cartella con counts.tsv
R - Lettura counts	Pattern	_counts.tsv$	Seleziona solo files counts
R - Merge counts	NA replacement	0	Sostituisce NA con 0
R - Metadata	File	metadata.tsv	Contiene sample e condition
R - Metadata	Gruppi	SANO, PB, PSO, DA	Livelli factor
DGEList / edgeR	Counts	counts	Matrice counts miRNA
DGEList / edgeR	Group	group	Factor con condizioni
DGEList / edgeR	Filtro	filterByExpr	Rimuove espressione bassa
DGEList / edgeR	Normalizzazione	calcNormFactors (TMM)	Standard in miRNA-seq
R - PCA / MDS	logCPM	cpm(dge, log=TRUE, prior.count=1)	Trasformazione logaritmica
PCA / MDS	Colori	group	Visualizzazione per condizione
R - voom / limma	Design matrix	~0 + group	Linear model senza intercept
voom	QualityWeights	sì	voomWithQualityWeights
Limma	Contrasts	PB_vs_SANO, PSO_vs_SANO, DA_vs_SANO	Comparazioni specifiche
Limma	eBayes	sì	Moderated t-statistics
R - Volcano plot	x-axis	logFC	Fold change log2
Volcano plot	y-axis	adj.P.Val	FDR-corrected p-value
Volcano plot	pCutoff	0.05	Threshold FDR
Volcano plot	FCcutoff	1	Threshold log2 fold change
R - Heatmap	Top miRNA	top 30	MiRNA significativi
Heatmap	Scaling	row	Z-score per riga
Heatmap	Annotation	Condition	Colonna meta per colore
R - Analisi globale pruritus	Gruppi	SANO vs PRURITUS	Combinazione PB+PSO+DA
Global DGEList	Filter	filterByExpr	Rimuove espressione bassa
Global DGEList	Normalization	calcNormFactors	TMM
Global voom/limma	Design	~0 + group2	Factor binario SANO vs PRURITUS
Global contrasts	PRURITUS_vs_SANO	PRURITUS - SANO	DE per firma globale
R - UpSet plot	Input	significant miRNA per contrasto	PB, PSO, DA
UpSet plot	Order	freq	Frequenza di intersezioni
Output grafici	Path	results/plots/	PNG plots PCA, MDS, Volcano, Heatmap, UpSet
