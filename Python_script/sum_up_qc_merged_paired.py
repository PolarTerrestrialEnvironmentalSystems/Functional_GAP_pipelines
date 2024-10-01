# Author: Ugur Cabuk
#!/bin/python3
import os

def sum_up_num_reads(file1_path, file2_path, output_path):
    # Read for the first file
    data1 = {}
    with open(file1_path, 'r') as file1:
        next(file1)  # skip the header
        for line in file1:
            fields = line.strip().split('\t')
            name = fields[0]
            num_reads = float(fields[4])
            data1[name] = num_reads
    
    with open(file2_path, 'r') as file2: #second file
        with open(output_path, 'w') as output_file:
            header = next(file2).strip()
            output_file.write(header + '\n') 
            for line in file2:
                fields = line.strip().split('\t')
                name = fields[0]
                num_reads = float(fields[4])
                if name in data1:
                    sum_num_reads = data1[name] + num_reads
                    fields[4] = str(sum_num_reads)
                output_file.write('\t'.join(fields) + '\n')

def process_directory(directory):
    for file_name in os.listdir(directory):
        if file_name.endswith("_merged"):
            merged_file_path = os.path.join(directory, file_name, "quant.sf")
            paired_file_path = os.path.join(directory, file_name.replace("_merged", "_paired"), "quant.sf")
            output_file_path = os.path.join(directory, file_name.replace("_merged", "_merged_paired.quant.sf"))
            sum_up_num_reads(merged_file_path, paired_file_path, output_file_path)

directory_path = input("Enter the path to the directory containing merged and paired files: ")
process_directory(directory_path)
