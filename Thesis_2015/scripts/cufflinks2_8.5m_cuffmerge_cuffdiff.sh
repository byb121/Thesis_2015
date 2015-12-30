#!/bin/bash
#$ -V -cwd -j y
#$ -pe mpi 1


module purge

module load mpi/openmpi/1.6.2/gcc-4.4.6
module load apps/samtools/0.1.18/gcc-4.4.6
module load apps/cufflinks/2.0.2/gcc-4.4.6+boost-1.49.0+samtools-0.1.18+eigen-3.0.5

cuffmerge -o tophat2_cuff202_ref_8.5m_cuffmerge_all -p 8 \
-g ~/lustre/Yaobo/GenomeData/ensembl66_transcripts_hg19.gtf \
-s ~/lustre/Yaobo/GenomeData/hg19.fa \
tophat_cufflinks202_no_ref_8.5m_assembly_sample_list.txt

cuffdiff -o tophat2_cufflinks_8.5m_cuffdiff -L OA,NOF -p 8 -N -b ~/lustre/Yaobo/GenomeData/hg19.fa \
-u \
-M ~/lustre/Yaobo/GenomeData/ensembl66_rRNA_transcripts_hg19.gtf \
tophat2_cuff202_ref_8.5m_cuffmerge_all/merged.gtf \
tophat_out_s1/accepted_hits.bam,\
tophat_out_s3/accepted_hits.bam,\
tophat_out_s4/accepted_hits.bam,\
tophat_out_s6/accepted_hits.bam,\
tophat_out_s8/accepted_hits.bam,\
tophat_out_s9/accepted_hits.bam,\
tophat_out_s11/accepted_hits.bam,\
tophat_out_s12/accepted_hits.bam,\
tophat_out_s15/accepted_hits.bam,\
tophat_out_s16/accepted_hits.bam \
tophat_out_s2/accepted_hits.bam,\
tophat_out_s5/accepted_hits.bam,\
tophat_out_s7/accepted_hits.bam,\
tophat_out_s10/accepted_hits.bam,\
tophat_out_s13/accepted_hits.bam,\
tophat_out_s14/accepted_hits.bam
