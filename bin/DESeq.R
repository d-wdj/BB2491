directory <- "../data/processed_raw_counts.tsv"
data <- read.table(directory, sep='\t', header=TRUE, row.names=1)
N <- nrow(data)

healthy <- data[,grepl("H",names(data))]
steatosis <- data[,grepl("S",names(data))]
nash <- data[,grepl("N",names(data))]
cancer <- data[,grepl("C",names(data))]

source("http://www.bioconductor.org/biocLite.R")
library("DESeq")

total_data <- cbind(healthy, steatosis, nash, cancer)
total_data <- round(total_data,0)
col_names <- colnames(total_data)
condition <- c(rep("H", length(healthy)), rep("S", length(steatosis)), rep("N", length(nash)), rep("C", length(cancer)))
(coldata <- data.frame(row.names=colnames(total_data), condition))

#dds <- DESeqDataSetFromMatrix(countData=total_data, colData=coldata, design=~condition)
#dds <- DESeq(dds)

cds <- newCountDataSet(total_data, condition)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res_HS <- nbinomTest(cds, "H", "S")
res_SN <- nbinomTest(cds, "S", "N")
res_NC <- nbinomTest(cds, "N", "C")

write.table(res_HS, file="../results/HS_DESeq.tsv", sep='\t')
write.table(res_SN, file="../results/SN_DESeq.tsv", sep='\t')
write.table(res_NC, file="../results/SC_DESeq.tsv", sep='\t')

###ARCHIVED 2018.01.15
