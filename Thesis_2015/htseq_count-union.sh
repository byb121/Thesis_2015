#!/bin/bash
#$ -cwd
#$ -j y
#$ -V
 
module load apps/samtools/0.1.18/gcc-4.4.6
module load libs/numpy/1.6.2/gcc-4.4.6+python-2.7.3+atlas-3.10.0
module load libs/pysam/0.7.7/gcc-4.4.6+python-2.7.3

IN_BAM=$1
GTF=$2
OUT_COUNT=$3

#for tophat alignment mapQ > 3 means unique alignment 
#samtools sort -no $IN_BAM "$IN_BAM.temp" \
#| \
python -m HTSeq.scripts.count -f bam -r name -s no -a 4 -t exon -i gene_id -m union \
$IN_BAM \
$GTF \
> $OUT_COUNT

echo "Done."

