#!/usr/bin/env python

import argparse
import pandas as pd

# Argument parser to accept input files and output file
parser = argparse.ArgumentParser(description="Extract top two rows from MASH results and save as TSV.")
parser.add_argument("-i", "--input", nargs="+", required=True, help="List of MASH results files (space-separated)")
parser.add_argument("-o", "--output", required=True, help="Output TSV file")
args = parser.parse_args()

# Prepare a list to store extracted data
summary_data = []

# Process each input file
for file_path in args.input:
    sample_name = file_path.split(".")[0]  # Extract sample name (before first dot)

    # Read the first two lines from the file
    with open(file_path, "r") as f:
        lines = [line.strip().split("\t") for line in f.readlines()[:2]]  # Read & split by tab

    # Extract relevant columns (first two + last column)
    first_row = True
    for line in lines:
        if len(line) >= 3:  # Ensure the line has enough columns
            identity, coverage, organism = line[0], line[1], line[-1]
            if first_row:
                summary_data.append([sample_name, identity, coverage, organism])
                first_row = False
            else:
                summary_data.append(["", identity, coverage, organism])  # Blank sample name for subsequent rows
    
    # Add an extra blank line after each sample block
    summary_data.append(["", "", "", ""])

# Convert to DataFrame and save as TSV
df = pd.DataFrame(summary_data, columns=["Sample", "Identity", "Coverage", "Organism Found"])
df.to_csv(args.output, sep="\t", index=False, header=True)

print(f"Summary saved to {args.output}")
