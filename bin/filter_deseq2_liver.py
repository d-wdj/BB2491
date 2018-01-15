

### NC
NC_old = open("../results/DESeq2/NC_DESeq2_genesymb_FDR5%"+".tsv", 'r')
NC_filtered = open("../results/DESeq2/NC_DESeq2_genesymb_FDR5%_filtered.tsv", 'w')

NC = NC_old.readlines()
NC_filtered.write(NC[0])

for lines in NC[1:]:
    lines = lines.split('\t')
    if "liver" in lines[2]:
        NC_filtered.write('\t'.join(lines))

NC_old.close()
NC_filtered.close()


### SN
SN_old = open("../results/DESeq2/SN_DESeq2_genesymb_FDR5%"+".tsv", 'r')
SN_filtered = open("../results/DESeq2/SN_DESeq2_genesymb_FDR5%_filtered.tsv", 'w')

SN = SN_old.readlines()
SN_filtered.write(SN[0])

for lines in SN[1:]:
    lines = lines.split('\t')
    if "liver" in lines[2]:
        SN_filtered.write('\t'.join(lines))

SN_old.close()
SN_filtered.close()


### HC
HC_old = open("../results/DESeq2/HS_DESeq2_genesymb_FDR5%"+".tsv", 'r')
HC_filtered = open("../results/DESeq2/HC_DESeq2_genesymb_FDR5%_filtered.tsv", 'w')

HC = HC_old.readlines()
HC_filtered.write(HC[0])

for lines in HC[1:]:
    lines = lines.split('\t')
    if "liver" in lines[2]:
        HC_filtered.write('\t'.join(lines))

HC_old.close()
HC_filtered.close()

print ("Done!")

###ARCHIVED 2018.01.15
