#! /bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe smp 4

#bowtie2-2.1.0
#tophat-2.0.10.Linux_x86_64

module load apps/samtools/0.1.18/gcc-4.4.6

READ_1_FILE=$1
READ_2_FILE=$2
SAMPLE_DIR=$3

OUTPUT_DIR="${SAMPLE_DIR}/tophat2"
OUTPUT_BAM="${SAMPLE_DIR}/tophat2/accepted_hits.bam"

echo `date`
echo "Tophat is starting."
mkdir ${OUTPUT_DIR}
tophat2 \
-p 4 \
--library-type fr-unstranded \
-o ${OUTPUT_DIR} \
-r 150 \
-g 1 \
--solexa1.3-quals \
--coverage-search \
--microexon-search \
--bowtie-n \
--b2-sensitive \
--transcriptome-index /users/a5907529/lustre/Yaobo/GenomeData/tophat_hg19_gencode14_index/hg19_gencode.V14 \
-x 20 \
/users/a5907529/lustre/Yaobo/GenomeData/bowtie_index/hg19 ${READ_1_FILE} ${READ_2_FILE}
echo `date` 
echo "Tophat is done."

echo `date` 
echo "Indexing bam file...."
samtools index ${OUTPUT_BAM}
echo `date` 
echo "Indexing is done."

