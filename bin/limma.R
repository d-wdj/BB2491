directory <- "../data/processed_raw_counts.tsv"
data <- read.table(directory, sep='\t', header=TRUE, row.names=1)
N <- nrow(data)

library("limma")
library("edgeR")

healthy <- data[,grepl("H",names(data))]
steatosis <- data[,grepl("S",names(data))]
nash <- data[,grepl("N",names(data))]
cancer <- data[,grepl("C",names(data))]

HS <- cbind(healthy, steatosis)
dge <- DGEList(counts=HS)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM)
fit <- eBayes(fit, trend=TRUE)
HS_fit <- data.frame(fit)
write.table(HS_fit, file="../results/HS_limma.tsv", sep='\t')

SN <- cbind(steatosis, nash)
dge <- DGEList(counts=SN)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM)
fit <- eBayes(fit, trend=TRUE)
SN_fit <- data.frame(fit)
write.table(SN_fit, file="../results/SN_limma.tsv", sep='\t')

NC <- cbind(nash, cancer)
dge <- DGEList(counts=NC)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM)
fit <- eBayes(fit, trend=TRUE)
NC_fit <- data.frame(fit)
write.table(NC_fit, file="../results/NC_limma.tsv", sep='\t')