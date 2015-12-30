#!/bin/bash
#$ -cwd
#$ -j y
#$ -V

module load apps/samtools/0.1.18/gcc-4.4.6
module load libs/numpy/1.6.2/gcc-4.4.6+python-2.7.3+atlas-3.10.0

IN_BAM=$1

NSORTED_PREFIX="${IN_BAM}.nsorted"
NSORTED="${NSORTED_PREFIX}.bam"

samtools sort -n ${IN_BAM} ${NSORTED_PREFIX}
#samtools index ${NSORTED}
