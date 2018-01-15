import matplotlib.pyplot as plt
import numpy as np
from math import log10# as log10

def collect_pval(filepath):
    print ("Processing {}...".format(filepath))
    qval = list()
    with open(filepath, 'r') as op:
        op = op.readlines()
        for line in op[1:]:
            line = line.split()
            if "NA" not in line[-1]:
                qval.append(log10(float(line[-1])))
    return qval


HS = collect_pval("../results/DESeq2/HS_DESeq2_genesymb"+".tsv")
SN = collect_pval("../results/DESeq2/SN_DESeq2_genesymb"+".tsv")
NC = collect_pval("../results/DESeq2/NC_DESeq2_genesymb"+".tsv")
# NCL = collect_pval("../results/DESeq2/NC_Deseq2_genesymb_FDR5%_filtered.tsv")
# print (HS, SN, NC[:20], NCL[:20])

plt.figure()
plt.hist(HS, histtype='step',color='b', label='H-S')
plt.hist(SN, histtype='step',color='r', label='S-N')
# bins = np.arange(min(NC), max(NC) + 10, 10)
plt.hist(NC, histtype='step',color='c',  label='N-C')
# plt.hist(NCL, histtype='step',color='m',  label='N-C (liver only)')

# plt.xscale('log')
plt.xlabel("log10(q-value)")
# plt.xlim(plt.xlim()[::-1])
plt.ylabel("Frequency")
plt.legend(loc='upper left')
plt.title("Histogram of log10 of Adjusted p-values from DESeq2 Results")
plt.savefig("../pictures/deseq2_qval_histogram_all.png", bbox_inches='tight')
plt.show()

###ARCHIVED 2018.01.15
