# List of Scripts and What They Do  
1. Data Pre-processing  
* script1  
* script2  
* scriptn  

2. DESeq2
* *DESeq2.Rmd*: Takes in the processed raw_count files, appropriate the input for DESeq2 in R then generate .tsv files containing the statistical analyses across the different disease states.  

3. Post-DESeq2  
* prep_piano.py: Converts the Ensembl gene IDs to HPA gene symbols. Also filters down the results to only the protein-coding genes.  
