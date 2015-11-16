

library(BitSeq)



NOFfiles <- c("f1_s2_bowtie12.7.R.rpkm", "f1_s5_bowtie12.7.R.rpkm", "f1_s7_bowtie12.7.R.rpkm", "f2_s2_bowtie12.7.R.rpkm", "f2_s5_bowtie12.7.R.rpkm", "f2_s6_bowtie12.7.R.rpkm")

 

OAfiles <- c("f1_s1_bowtie12.7.R_2nd.rpkm", "f1_s3_bowtie12.7.R.rpkm", "f1_s4_bowtie12.7.R.rpkm", "f1_s6_bowtie12.7.R.rpkm", "f1_s8_bowtie12.7.R.rpkm", "f2_s1_bowtie12.7.R.rpkm", "f2_s3_bowtie12.7.R.rpkm", "f2_s4_bowtie12.7.R.rpkm", "f2_s7_bowtie12.7.R.rpkm" , "f2_s8_bowtie12.7.R_2nd.rpkm")



getMeanVariance(c(NOFfiles,OAfiles),outFile = "data.means",log = TRUE)



estimateHyperPar( outFile = "data.par", cond1 = NOFfiles, cond2 = OAfiles, meanFile = "data.means", verbose = TRUE )



estimateDE(NOFfiles, OAfiles, outFile = "data", parFile = "data.par", samples = TRUE )

table <- read.table("data.pplr", sep = " ", header = F)

# seems any .tr file will do
annotations <- read.table("f1_s1_bowtie12.7.R_2nd.ensemblGenes.tr", header =F, sep = " ")  

y <- as.character(annotations$V2)

rownames(table) <- lapply(y, function(x) strsplit(x, "_")[[1]][3])


table$Transcripts <- rownames(table)


colnames(table) <- c("PPLR", "Confidence_low", "Confidence_high", "log2_FoldChange", "Mean_C1", "Mean_C2", "Transcripts_ID")


sig_table <- table[(table$PPLR >= 0.95 | table$PPLR <= 0.05) & (table$log2_FoldChange >= 1 | table$log2_FoldChange <= -1),]
length(sig_table$PPLR)

sig_sig_table <- table[(table$PPLR >= 0.95 | table$PPLR <= 0.05) & (table$log2_FoldChange >= 1 | table$log2_FoldChange <= -1) & (table$Mean_C1 >=1 | table$Mean_C2 >= 1),]
length(sig_sig_table$PPLR)


#write.table(table, "BitSeq_output.txt", row.name = T, quote = F, col.name = F, sep = "\t")
