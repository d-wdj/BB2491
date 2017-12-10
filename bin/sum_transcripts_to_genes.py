import sys
import pandas as pd
import numpy as np

### Step 1: Load the dictionary
# tID_dict = dict()
# with open('../reference/ensembl_ref_dict.txt', 'r') as en:
#     en = en.readlines()
#     print ("Generating dictionary...")
#     for line in en[1:]:
#         val = line.split()
#         tID_dict[line[2]] = val
#     print ("Done.")

## Example dictionary key:value
## {'ENST00000601199': ['FAM231D', 'ENSG00000268674', 'NA']}
## {transcriptID: [g_Name, g_ID, tsl]}

### Step 2: Map the transcript ID into the gene ID

# with open("../data/raw_counts_NASH_code.txt", 'r') as ns:
#     print ("Mapping transcript ID to gene ID...")
#     RT = open("../data/raw_counts_NASH_geneID.txt", 'w')
#     ns = ns.readlines()
#     RT.write(ns[0])
#     transcript_not_found = 0
#     total_transcript = 0
#     for line in ns[1:]:
#         total_transcript += 1
#         # print (line)
#         enst = (line.split())
#         # print (enst)
#         if enst[0] in tID_dict:
#             # print (enst, tID_dict[enst[0]][2], enst[1:])
#             RT.write(tID_dict[enst[0]][1]+'\t'+'\t'.join(enst[1:])+'\n')
#         else:
#             transcript_not_found += 1
#             # print ("Transcript ID: {} not found.".format(enst[0]))
#     RT.close()
#     print ("{}/{} transcripts not found.".format(transcript_not_found,
#                                                 total_transcript))
#     print ("Done!")

### Step 3: Load the file into a pandas dataframe. Then do (1) remove low-Q
### samples (i.e. NASH code X, U); (2) merge sample replicates; (3) sum values
### for the same gene IDs.
data = pd.read_table("../data/1k_NASH_gene_id.txt", index_col=0, header=0)
data = data.groupby(data.index, sort=True).sum()
# print (data.groupby(data.index, sort=True).sum())



#
# ### Step 3: Possibly use awk(sed?) to replace the first column, i.e. the tID
# ### into the gene ID?
# # key itemgetter groupby? to group by gene and sum them tgt
# ## Trying pandas dataframe...
