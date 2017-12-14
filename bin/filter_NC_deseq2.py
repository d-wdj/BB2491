
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

print ("Done!")
