#!/usr/bin/env python

import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Extract top N rows from MASH results and save as TSV.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="List of MASH results files")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument("-n", "--num_lines", type=int, default=3, 
                       help="Number of lines to output per file (default: 2)")
    
    args = parser.parse_args()
    summary_data = []

    for file_path in args.input:
        try:
            sample_name = Path(file_path).stem  # More robust filename extraction
            
            with open(file_path, "r") as f:
                # Read only needed lines (memory efficient)
                lines = []
                for _ in range(args.num_lines):
                    line = next(f, None)
                    if line is None:
                        break
                    lines.append(line.strip().split("\t"))
                
                # Process lines if they have enough columns
                for i, line in enumerate(lines):
                    if len(line) >= 3:
                        row = [
                            sample_name if i == 0 else "",  # Only show sample name in first row
                            line[0],  # identity
                            line[1],  # coverage
                            line[-1]  # organism
                        ]
                        summary_data.append(row)
            
            # Add separator only if we found data
            if lines:
                summary_data.append(["", "", "", ""])
                
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            continue

    if summary_data:
        df = pd.DataFrame(
            summary_data,
            columns=["Sample", "Identity", "Coverage", "Organism Found"]
        )
        df.to_csv(args.output, sep="\t", index=False)
        print(f"Successfully saved summary to {args.output}")
    else:
        print("No valid data found to save.")

if __name__ == "__main__":
    main()