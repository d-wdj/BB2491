
## Downloaded protein coding genes from Protein Atlas
## https://www.proteinatlas.org/about/download
### Step 1: Load the dictionary
gene_sym = dict()
with open('../reference/proteinatlas.tsv', 'r') as en:
    en = en.readlines()
    print ("Generating dictionary...")
    for line in en[1:100]:
        line = line.split('\t')
        key, val = line[2], line[0]
        gene_sym[key] = val
    print ("Done.")
