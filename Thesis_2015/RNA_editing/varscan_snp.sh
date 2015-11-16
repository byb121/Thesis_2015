#! /bin/bash
#$ -cwd
#$ -j y
# #$ -pe smp 2

# Load Modules

module load apps/samtools/0.1.18/gcc-4.4.6
module load apps/python/2.7.3/gcc-4.4.6
module load apps/perl/5.16.1/gcc-4.4.6

IN_Pileup_gzipped=$1
OUTPUT_VCF=$2

VarScan="java -jar /users/a5907529/archive/biosofts/VarScan.v2.3.6.jar mpileup2snp"

gunzip -c ${IN_Pileup_gzipped} | $VarScan \
--min-coverage 3 \
--min-reads2 2 \
--min-avg-qual 20 \
--min-var-freq 0.01 \
--p-value 0.99 \
--strand-filter 0 \
--output-vcf 1 \
> ${OUTPUT_VCF}

gzip ${OUTPUT_VCF}

echo "Done!"

