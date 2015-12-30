library(BitSeq)

id <- Sys.getenv("SGE_TASK_ID")
taskid <- as.numeric(id)

trSeqFile <- "ensemblGenes.fasta" 
samFile <- paste ("/users/a5907529/RNA_seq/BitSeq/f1_s", taskid, "_bowtie12.7.sam", sep = "")
alignmentProbabilityFile <- paste ("/users/a5907529/RNA_seq/BitSeq/f1_s", taskid, "_bowtie12.7.R.prob", sep = "")
trInfoFile = paste ("/users/a5907529/RNA_seq/BitSeq/f1_s", taskid, "_bowtie12.7.R.ensemblGenes.tr", sep = "")

## <for test
#taskid <- 1
#samFile <- paste ("/users/a5907529/RNA_seq/BitSeq/f1_s", taskid, "_bowtie12.7.sam", sep = "")
##for test>

parseAlignment( samFile,
outFile = alignmentProbabilityFile,
trSeqFile = trSeqFile,
trInfoFile = trInfoFile,
verbose = TRUE)

expressionFile <- paste ("f1_s", taskid, "_bowtie12.7.R", sep = "")

estimateExpression(alignmentProbabilityFile,
outFile = expressionFile, 
outputType = "RPKM",
trInfoFile = trInfoFile, 
MCMC_burnIn = 1000,
MCMC_samplesN = 1000, 
MCMC_samplesSave = 1000,
MCMC_scaleReduction = 1.2,
MCMC_chainsN = 4)

