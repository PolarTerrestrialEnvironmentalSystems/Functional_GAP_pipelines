# Author: Ugur Cabuk
# Description: This script calculates NGC (Normalized Gene Counts based on catalog size) values from raw counts data and gene lengths.

import pandas as pd

#Input_file for the pipeline
input_file= "prokGAP_all_lake_lama_gene_quant.raw_count.fixed.tsv""
# the number should be changed based on the gene catalog size.
scaling_factor=6568483
df = pd.read_csv(input_file, sep='\t')

#NGC calculation formula
def ngc(counts, length):
    x = counts / length
    return (x.T * scaling_factor / x.sum()).T

#extract the gene length
gene_lengths = df['gene_length']

for col in df.columns[1:-1]:  #Exclude the first and last column
    print(col)  # Print the current sample column name to double check.
    df[col] = ngc(df[col], gene_lengths)
    
# Save it.
output_file = "gene_abundance_cpm_values_prokgap.tsv"
df.to_csv(output_file, sep='\t', index=False)
