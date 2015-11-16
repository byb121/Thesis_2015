#!/bin/bash

#$ -V -cwd -j y -M yaobo.xu@ncl.ac.uk -m bes
#$ -t 1-8

i=$(expr $SGE_TASK_ID)
novoalign -d ~/GenomeData/hg19/hg19_rna_seq -f 's_'$i'_1_sequence.txt.trim' 's_'$i'_2_sequence.txt.trim' -F ILMFQ -K 'mismatch_statistics_transcrpts_s'$i'.txt' -o SAM -r Random -i PE 190,50 > 's_'$i'.trim.novoalign.transcripts.sam'


