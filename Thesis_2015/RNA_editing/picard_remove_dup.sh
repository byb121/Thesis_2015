#! /bin/bash
#$ -cwd
#$ -j y
#$ -V
# #$ -pe smp 4

module load apps/picard/1.85/noarch
module load apps/samtools/0.1.18/gcc-4.4.6

SAMPLE_DIR=$1
SAMPLE_ID=$2
IN_BAM="${SAMPLE_DIR}/tophat2/accepted_hits.bam"

WRKGDIR="/users/a5907529/scratch/${SAMPLE_ID}.GATKpipe.picard"
PICARD_OUTDIR="${SAMPLE_DIR}/picard"
OUT_BAM="${SAMPLE_DIR}/picard/${SAMPLE_ID}_no_dups.bam"
REF_FILE="/users/a5907529/lustre/Yaobo/GenomeData/GATK_bundle_2.5/ucsc.hg19.4GATK.fasta"

Picard_nodups="java -jar $PICARDDIR/java/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT"
PICARD_TEMP="$WRKGDIR/Picard_Temp"
PICARD_LOG="$PICARD_OUTDIR/${SAMPLE_ID}_picard.log"

mkdir $PICARD_OUTDIR
$Picard_nodups INPUT=${IN_BAM} OUTPUT=${OUT_BAM} METRICS_FILE=$PICARD_LOG REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=$PICARD_TEMP
samtools index ${OUT_BAM}

