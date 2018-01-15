import sys
import pandas as pd
import numpy as np

### Step 1: Load the dictionary
tID_dict = dict()
with open('../reference/ensembl_ref_dict.txt', 'r') as en:
    en = en.readlines()
    print ("Generating dictionary...")
    for line in en[1:]:
        val = line.split()
        tID_dict[val[2]] = val
    print ("Done.")

## Example dictionary key:value
## {'ENST00000601199': ['FAM231D', 'ENSG00000268674', 'NA']}
## {transcriptID: [g_Name, g_ID, tsl]}

### Step 2: Map the transcript ID into the gene ID
with open("../data/raw_counts_NASH_code.txt", 'r') as ns:
    print ("Mapping transcript ID to gene ID...")
    RT = open("../data/raw_counts_NASH_geneID.txt", 'w')
    ns = ns.readlines()
    RT.write(ns[0])
    transcript_not_found = 0
    total_transcript = 0
    for line in ns[1:]:
        total_transcript += 1
        # print (line)
        enst = (line.split())
        # print (enst)
        if enst[0] in tID_dict:
            # print (enst, tID_dict[enst[0]][2], enst[1:])
            RT.write(tID_dict[enst[0]][1]+'\t'+'\t'.join(enst[1:])+'\n')
        else:
            transcript_not_found += 1
            # print ("Transcript ID: {} not found.".format(enst[0]))
    RT.close()
    print ("{}/{} transcripts not found.".format(transcript_not_found,
                                                total_transcript))
    print ("Done!")

### Step 3: Load the file into a pandas dataframe. Then do (1) remove low-Q
### samples (i.e. NASH code X, U); (2) merge sample replicates; (3) sum values
### for the same gene IDs.
print ("Loading the data...")
data = pd.read_table("../data/raw_counts_NASH_geneID.txt", index_col=0, header=0)
data.drop([col for col in data.columns if 'U' in col or 'X' in col],
        axis=1,inplace=True)
print ("""Removing low-quality samples...
Merging and summing replicates...""")
data = data.groupby(data.index, sort=True).sum()
data = data.T.groupby([s.split('.')[0] for s in data.T.index.values]).sum().T
print ("Writing into file...")
data.to_csv(path_or_buf="../data/processed_raw_counts.tsv", sep='\t')
print ("Done!")

###ARCHIVED 2018.01.15
