#!/bin/bash
#$ -V -cwd -j y -t 1-8

i=$(expr $SGE_TASK_ID)

source /etc/profile.d/modules.sh 
module load apps/samtools/0.1.18/gcc-4.4.6

samtools view -h '/users/a5907529/archive/RNA_seq/FlowCell1/s_'$i'.gsnap20120620.trim.sorted.bam' > '/users/a5907529/lustre/Yaobo/DiffSplice/FlowCell1/s_'$i'.gsnap20120620.trim.sorted.sam'
 
