# Author: Ugur Cabuk
#!/bin/python3

#Extract the second column from the len_samples file
gene_lengths = []
with open("len_samples", 'r') as file:
    next(file)  # Skip the first line
    for line in file:
        gene_length = line.split()[1]
        gene_lengths.append(gene_length)

#"numreads_samples.txt", append gene lengths to each line
with open("numreads_samples", 'r') as input_file, open("numreads_sample_fixed.tsv", 'w') as output_file:
    header = input_file.readline().strip().split('\t') 

    #gene_length to the header
    header.append('gene_length')

    #Write the modified header
    output_file.write('\t'.join(header) + '\n')

    #gene lengths to each line
    for line, length in zip(input_file, gene_lengths):
        fields = line.strip().split('\t')

        #the gene length to the fields
        fields.append(length)
        output_file.write('\t'.join(fields) + '\n')
