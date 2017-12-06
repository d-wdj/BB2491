
def slice_reference():
    with open("../data/transcript.gtf", 'r') as t:
        t = t.readlines()
        # for z, line in enumerate(t):
        ref = open("../data/ensembl_ref_dict.txt", 'w')
        genes = "gene_name\tgene_length\tgene_id"
        transcript = "transcript_id\ttsl"
        ref.write(genes+'\t'+transcript+'\n')
        for line in t:
            line = line.split()
            for z, feature in enumerate(line):
                gene_length = str(int(line[4])+1-int(line[3]))
                if feature == 'gene_id':
                    gene_id = line[z+1].replace(";","").replace("\"", "")
                if feature == 'gene_name':
                    gene_name = line[z+1].replace(";","").replace("\"", "")
                if feature == 'transcript_id':
                    transcript_id = line[z+1].replace(";","").replace("\"", "")
                if feature == 'transcript_support_level':
                    tsl = line[z+1].replace(";","").replace("\"", "")
            ref.write(gene_name+'\t'+gene_length+'\t'+gene_id+'\t')
            ref.write(transcript_id+'\t'+tsl+'\n')
            # print (gene_id, gene_name, gene_length)
            # print (transcript_id, tsl)
        ref.close()
        print ("Done!")
    return True
# slice_reference()
