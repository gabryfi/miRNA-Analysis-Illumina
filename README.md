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
