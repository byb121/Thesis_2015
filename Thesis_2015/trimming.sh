#! /bin/bash
#$ -cwd
#$ -j y
#$ -V

module load apps/python
module load apps/perl

READ_1_FILE=$1
READ_2_FILE=$2
OUTPUT_DIR=$3

#trim_galore to remove adaptors and low quality tails and perform fastqc on the output
/users/a5907529/archive/biosofts/trim_galore/trim_galore -q 20 --paired -o ${OUTPUT_DIR} --phred64 --length 20 --stringency 5 --fastqc --dont_gzip ${READ_1_FILE} ${READ_2_FILE}

