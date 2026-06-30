#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:13:02 2025

@author: mozhganoroujlu
"""

import os

def get_bin(value, bins):
    for i, (start, end) in enumerate(bins):
        if start <= value < end:
            return i
    return None

# Hardcoded folder paths (replace these with your actual paths)
folder1 = 'normalized_contacts/bandnorm/bandnorm_txt/'  # update it 
folder2 = '/normalized_contacts/bandnorm/scKTLD_TADs/'  
output_folder = '/normalized_contacts/bandnorm/filtered_tads_contacts/'  

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# List files in folder1
for file_name in os.listdir(folder1):
    if file_name.endswith('.txt'):
        base_name = file_name.replace('.txt', '')
        file1_path = os.path.join(folder1, file_name)
        
        file2_name = base_name + '_tads.txt'
        file2_path = os.path.join(folder2, file2_name)
        
        if os.path.exists(file2_path):
            # Read bins from file2
            bins = []
            with open(file2_path, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:  # Skip header
                    parts = line.strip().split()
                    if len(parts) == 2:
                        start = int(parts[0])
                        end = int(parts[1])
                        bins.append((start, end))
            
            # Process file1 and write to output
            output_file_name = base_name + '_filtered.txt'
            output_path = os.path.join(output_folder, output_file_name)
            
            with open(file1_path, 'r') as f_in, open(output_path, 'w') as f_out:
                for line in f_in:
                    parts = line.strip().split('\t')
                    if len(parts) == 3:
                        col1 = int(parts[0])
                        col2 = int(parts[1])
                        value = parts[2]
                        if col1 == col2:
                            continue
                        bin1 = get_bin(col1, bins)
                        bin2 = get_bin(col2, bins)
                        if bin1 is not None and bin1 == bin2:
                            f_out.write(f"{col1}\t{col2}\t{value}\n")
