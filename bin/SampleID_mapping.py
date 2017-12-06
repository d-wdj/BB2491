import fileinput
import sys, os
## First load the metadata and save it as a dictionary:
metadata = dict()
with open("../../Data/metadata_patient_code.txt", 'r') as m:
    m = m.readlines()
    for line in m[1:]:
        line = line.rstrip('\n')
        line = line.split('\t')
        metadata[line[0]] = line[1]

## Second, load the raw count files and try to map the each sampleID into the
## correct cohort
# os.system("""head -n 50000 ../Data/raw_counts.txt > 50k_raw_counts.txt""")
with open("../../Data/raw_counts.txt", 'r') as r:
    with open ("../data/raw_counts_NASH_code.txt", 'w') as N:
        r = r.readlines()
        for z, line in enumerate(r):
            sampleID = list()
            if z == 0:
                line = line.strip('\n').split()
                sampleID.append('NASH_ID')
                for nash in line[1:]:
                    if "U" not in (metadata[nash]) and "X" not in (metadata[nash]):
                        sampleID.append(metadata[nash])
                    # else:
                    #     print (nash)
                sampleID = '\t'.join(sampleID)
                N.write(sampleID+'\n')
            else:
                N.write(line)


        # r = r.readlines()
        # header = r[0].strip('\n').split()
        # for tID in header[1:]:
        #     if "U" not in (metadata[tID]) and "X" not in (metadata[tID]):
        #         sampleID.append(metadata[tID])
        # sampleID.insert(0, header[0])
        # sampleID = "\t".join(sampleID)
        # with open ("raw_counts_NASH_code.txt", 'w') as N:
        #     N.write(sampleID+'\n')
        #     for line in r[1:]:
        #         N.write(line+'\n')
        #


#     r = r.read().split()
#     sampleID.append(r[0])
# # print(sampleID)
#     for tID in r[1:]:
# ## Know from metadata that sample that starts with U, X are low-qual, so exclude
#

#     with open("raw_counts_NASH_code.txt", 'w') as n:

## To replace the first line - sampleID with the ones you generated
## https://stackoverflow.com/questions/13438095/replace-the-first-line-in-a-text-file-by-a-string/13438118#13438118
# filename = "1k_raw_counts.txt"

# print(os.system("echo 'Splicing data...'"))
# command = """sed -i '' -e "1,$ s/.*/{$sampleID}/" '1k_raw_counts.txt'"""
# print(os.system(command))

# with open('1k_raw_counts', 'r+') as f:
    # repl_str = sampleID
    # lines = f.readlines()
    # f.seek(0)
    # f.truncate()
    # for l in lines:
    #     f.write(repl_str if l.startswith('NATOMS') else l)
