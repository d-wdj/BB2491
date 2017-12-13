library("piano")
kegg <- loadGSC("../reference/c2.cp.kegg.v6.1.symbols.gmt")
go <- loadGSC("../reference/c5.bp.v6.1.symbols.gmt")
onco <- loadGSC("../reference/c6.all.v6.1.symbols.gmt")


#### DESeq
print ("DESeq")
### Healthy-Steatosis
print ("HEALTHY-STEATOSIS")
HS_dir <- "../results/DESeq/HS_DESeq_genesymb.tsv"
HS <- read.table(HS_dir, sep='\t', header=TRUE) 
HS <- HS[complete.cases(HS), ]
HS <- HS[!duplicated(HS$gene_symbol), ]
HS <- data.frame(HS)

fc <- HS$log2FoldChange
HS_stat <- HS$pval
names(HS_stat) <- paste(HS$gene_symbol)
names(fc) <- paste(HS$gene_symbol)

gsaRes_HS_kegg <- runGSA(HS_stat, fc, gsc=kegg, geneSetStat="reporter", 
                 signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_HS_kegg, save=TRUE, file="../results/DESeq/HS_kegg_DESeq_piano.xls")
gsaRes_HS_go <- runGSA(HS_stat, fc, gsc=go, geneSetStat="reporter", 
                 signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_HS_go, save=TRUE, file="../results/DESeq/HS_go_DESeq_piano.xls")
gsaRes_HS_onco <- runGSA(HS_stat, fc, gsc=onco, geneSetStat="reporter", 
                 signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_HS_onco, save=TRUE, file="../results/DESeq/HS_onco_DESeq_piano.xls")

### Steatosis-NASH
print ("HEALTHY-STEATOSIS")
SN_dir <- "../results/DESeq/SN_DESeq_genesymb.tsv"
SN <- read.table(SN_dir, sep='\t', header=TRUE) 
SN <- SN[complete.cases(SN), ]
SN <- SN[!duplicated(SN$gene_symbol), ]
SN <- data.frame(SN)

fc <- SN$log2FoldChange
SN_stat <- SN$pval
names(SN_stat) <- paste(SN$gene_symbol)
names(fc) <- paste(SN$gene_symbol)

gsaRes_SN_kegg <- runGSA(SN_stat, fc, gsc=kegg, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_SN_kegg, save=TRUE, file="../results/DESeq/SN_kegg_DESeq_piano.xls")
gsaRes_SN_go <- runGSA(SN_stat, fc, gsc=go, geneSetStat="reporter", 
                       signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_SN_go, save=TRUE, file="../results/DESeq/SN_go_DESeq_piano.xls")
gsaRes_SN_onco <- runGSA(SN_stat, fc, gsc=onco, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_SN_onco, save=TRUE, file="../results/DESeq/SN_onco_DESeq_piano.xls")

### NASH-HCC
print ("HEALTHY-STEATOSIS")
NC_dir <- "../results/DESeq/NC_DESeq_genesymb.tsv"
NC <- read.table(NC_dir, sep='\t', header=TRUE) 
NC <- NC[complete.cases(NC), ]
NC <- NC[!duplicated(NC$gene_symbol), ]
NC <- data.frame(NC)

fc <- NC$log2FoldChange
NC_stat <- NC$pval
names(NC_stat) <- paste(NC$gene_symbol)
names(fc) <- paste(NC$gene_symbol)

gsaRes_NC_kegg <- runGSA(NC_stat, fc, gsc=kegg, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_NC_kegg, save=TRUE, file="../results/DESeq/NC_kegg_DESeq_piano.xls")
gsaRes_NC_go <- runGSA(NC_stat, fc, gsc=go, geneSetStat="reporter", 
                       signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_NC_go, save=TRUE, file="../results/DESeq/NC_go_DESeq_piano.xls")
gsaRes_NC_onco <- runGSA(NC_stat, fc, gsc=onco, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_NC_onco, save=TRUE, file="../results/DESeq/NC_onco_DESeq_piano.xls")


##### Limma
print ("Limma")
HS_dir <- "../results/Limma/HS_Limma_genesymb.tsv"
HS <- read.table(HS_dir, sep='\t', header=TRUE) 
HS <- HS[complete.cases(HS), ]
HS <- HS[!duplicated(HS$gene_symbol), ]
HS <- data.frame(HS)

fc <- HS$logFC
HS_stat <- HS$P.Value
names(HS_stat) <- paste(HS$gene_symbol)
names(fc) <- paste(HS$gene_symbol)

gsaRes_HS_kegg <- runGSA(HS_stat, fc, gsc=kegg, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_HS_kegg, save=TRUE, file="../results/Limma/HS_kegg_Limma_piano.xls")
gsaRes_HS_go <- runGSA(HS_stat, fc, gsc=go, geneSetStat="reporter", 
                       signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_HS_go, save=TRUE, file="../results/Limma/HS_go_Limma_piano.xls")
gsaRes_HS_onco <- runGSA(HS_stat, fc, gsc=onco, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_HS_onco, save=TRUE, file="../results/Limma/HS_onco_Limma_piano.xls")

### Steatosis-NASH
print ("HEALTHY-STEATOSIS")
SN_dir <- "../results/Limma/SN_Limma_genesymb.tsv"
SN <- read.table(SN_dir, sep='\t', header=TRUE) 
SN <- SN[complete.cases(SN), ]
SN <- SN[!duplicated(SN$gene_symbol), ]
SN <- data.frame(SN)

fc <- SN$logFC
SN_stat <- SN$P.Value
names(SN_stat) <- paste(SN$gene_symbol)
names(fc) <- paste(SN$gene_symbol)

gsaRes_SN_kegg <- runGSA(SN_stat, fc, gsc=kegg, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_SN_kegg, save=TRUE, file="../results/Limma/SN_kegg_Limma_piano.xls")
gsaRes_SN_go <- runGSA(SN_stat, fc, gsc=go, geneSetStat="reporter", 
                       signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_SN_go, save=TRUE, file="../results/Limma/SN_go_Limma_piano.xls")
gsaRes_SN_onco <- runGSA(SN_stat, fc, gsc=onco, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_SN_onco, save=TRUE, file="../results/Limma/SN_onco_Limma_piano.xls")

### NASH-HCC
print ("HEALTHY-STEATOSIS")
NC_dir <- "../results/Limma/NC_Limma_genesymb.tsv"
NC <- read.table(NC_dir, sep='\t', header=TRUE) 
NC <- NC[complete.cases(NC), ]
NC <- NC[!duplicated(NC$gene_symbol), ]
NC <- data.frame(NC)

fc <- NC$logFC
NC_stat <- NC$P.Value
names(NC_stat) <- paste(NC$gene_symbol)
names(fc) <- paste(NC$gene_symbol)

gsaRes_NC_kegg <- runGSA(NC_stat, fc, gsc=kegg, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_NC_kegg, save=TRUE, file="../results/Limma/NC_kegg_Limma_piano.xls")
gsaRes_NC_go <- runGSA(NC_stat, fc, gsc=go, geneSetStat="reporter", 
                       signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_NC_go, save=TRUE, file="../results/Limma/NC_go_Limma_piano.xls")
gsaRes_NC_onco <- runGSA(NC_stat, fc, gsc=onco, geneSetStat="reporter", 
                         signifMethod="nullDist", nPerm=1000, verbose=TRUE)
GSAsummaryTable(gsaRes_NC_onco, save=TRUE, file="../results/Limma/NC_onco_Limma_piano.xls")