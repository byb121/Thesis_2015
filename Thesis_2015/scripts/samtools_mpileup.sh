#! /bin/bash
#$ -cwd
#$ -j y
#$ -V
# #$ -pe smp 4

module load apps/samtools/0.1.18/gcc-4.4.6

BAM_LIST=$1
OUTPUT_BAM=$2
CHR=$3

REF="/users/a5907529/lustre/Yaobo/GenomeData/hg19.fa"

echo `date`
echo "mpileup is starting."

samtools mpileup -B -d 8000 -b ${BAM_LIST} -r $CHR -D -f $REF > ${OUTPUT_BAM}
gzip ${OUTPUT_BAM}

echo `date` 
echo "pileup for $CHR is done."

