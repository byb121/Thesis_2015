#!/bin/bash

#$ -V -cwd -j y -M yaobo.xu@ncl.ac.uk -m bes
#$ -t 1-8

i=$(expr $SGE_TASK_ID)

gsnap -d hg19 -m 5 -a off --trim-mismatch-score 0 --trim-indel-score 0 -v hg19.snp --novelsplicing 1 --use-splicing hg19.splicesites --pairexpect 190 --pairdev 50 --quality-protocol illumina -n 1 -Q --nofails -A sam raw_data/'s_'$i'_1_sequence.txt.trim' raw_data/'s_'$i'_2_sequence.txt.trim' > 's_'$i'.gsnap20120507.trim.sam'

echo "alignment that outputs sam format is done!"

gsnap -d hg19 -m 5 -a off --trim-mismatch-score 0 --trim-indel-score 0 -v hg19.snp --novelsplicing 1 --use-splicing hg19.splicesites --pairexpect 190 --pairdev 50 --quality-protocol illumina -n 1 -Q --nofails raw_data/'s_'$i'_1_sequence.txt.trim' raw_data/'s_'$i'_2_sequence.txt.trim' > 's_'$i'.gsnap20120507.gsnap_format.trim.txt'

echo "alignment that outputs gsnap format is done!"

echo "all done!"
