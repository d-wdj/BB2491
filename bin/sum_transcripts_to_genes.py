import sys


### Step 1: Load the dictionary
tID_dict = dict()
with open('../data/ensembl_ref_dict.txt', 'r') as en:
    en = en.readlines()
    for line in en[1:]:
        val = []
        line = line.split()
        val.extend(line[:3])
        val.extend(line[-1])
        tID_dict[line[3]] = val
# print (tID_dict)

## Example dictionary key:value
## {'ENST00000601199': [['FAM231D', '510', 'ENSG00000268674'], 'NA']}
## {transcriptID: [g_Name, g_length, g_ID, tsl]}

### Step 2: Map the transcript ID into the gene ID
with open("../data/raw_counts_NASH_code.txt", 'r') as ns:
    RT = open("../data/raw_NASH_geneID.txt", 'w')
    ns = ns.readlines()
    RT.write(ns[0]+'\n')
    transcript_not_found = 0
    total_transcript = 0
    for line in ns[1:]:
        total_transcript += 1
        # print (line)
        enst = (line.split())
        # print (enst)
        if enst[0] in tID_dict:
            # print (enst, tID_dict[enst[0]][2], enst[1:])
            RT.write(tID_dict[enst[0]][2] + '\t'.join(enst[1:]) + '\n')
        else:
            transcript_not_found += 1
            # print ("Transcript ID: {} not found.".format(enst[0]))
    RT.close()
    print ("{}/{} transcripts not found.".format(transcript_not_found,
                                                total_transcript))
    print ("Done!")


### Step 3: Possibly use awk(sed?) to replace the first column, i.e. the tID
### into the gene ID?
# key itemgetter groupby? to group by gene and sum them tgt
