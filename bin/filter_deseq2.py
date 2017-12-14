import glob

### Step 1: Create the dictionary of the localisation of expression in the
### tissue. Database: Human Protein Atlas
## Expect --> gene_symb == tissue location(s)
loc = dict()
with open("../reference/proteinatlas.tsv", 'r') as pt:
    pt = pt.readlines()
    countall = 0
    countdup = 0
    for line in pt[1:]: #Because line[0] is just the header
        countall += 1
        line = line.split('\t')
        # There are many genes with unknown localisation so we will discard them
        # The column 'RNA Tissue Category' might be interesting so will save it.
        # The reasoning is that we may not have time to investigate all the
        # genes especially if they are housekeeping genes.
        if line[17] == "":
            line[17] = "NA"
        loc[line[0]] = [line[14], line[17]]
print ("Dictionary generated.")

### Step 2: Open the DESeq2 output file, filter for the adj-p-val FDR 5% then
### filter for tissue localisation.
dir = "../results/DESeq2/"
pattern = "_DESeq2_genesymb.tsv"
newname = "_DESeq2_genesymb_FDR5%"+".tsv"
for dsq in glob.glob(dir+"??_DESeq2_genesymb.tsv"):
    print ("Processing: {}...".format(dsq))
    check = dsq.replace(dir,"").replace(pattern, "")
    newfile = dir+check+newname
    newfile = open(newfile, 'w')
    with open(dsq, 'r') as deseq2:
        deseq2 = deseq2.readlines()
        header = deseq2[0].strip('\n').split()
        h1 = header[0]+'\t'+'\"RNA_tissue_category\"'+'\t'
        h1 += '\"RNA_TS_TPM\"'+'\t'
        h2 = '\t'.join(header[1:])+'\n'
        newfile.write(h1+h2)
        for line in deseq2[1:]:
            line = line.split()
            if line[-1] == "NA" or line[-1] == "":
                continue
            elif float(line[-1]) < 0.05:
                genesymb = line[0]
                if genesymb in loc:
                    newline = line[0]+'\t'+loc[genesymb][0]
                    newline += '\t'+loc[genesymb][1]
                    newline += '\t'+'\t'.join(line[1:])+'\n'
                    newfile.write(newline)
    newfile.close()






print ("Done.")
