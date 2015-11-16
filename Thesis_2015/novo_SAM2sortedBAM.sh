#!/bin/bash

#$ -V -cwd -j y -t 1-8
i=$(expr $SGE_TASK_ID)

~/biosofts/samtools-0.1.18/samtools view -S -b 's_'$i'.trim.novoalign.transcripts.sam' | ~/biosofts/samtools-0.1.18/samtools sort - 's_'$i'.trim.novoalign.transcripts.sorted'
~/biosofts/samtools-0.1.18/samtools index 's_'$i'.trim.novoalign.transcripts.sorted.bam'
#~/biosofts/samtools-0.1.18/samtools pileup -s -f ~/GenomeData/hg19/hg19.fa 's_'$i'.final_alignment.v1.sorted.bam' > 's_'$i'.final_alignment.v1.sorted.pileup'

