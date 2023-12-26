#!/usr/bin/env python
# Author: Joseph JOSEPHIDES
# Institut Curie, Paris

import os
import argparse
import pandas as pd

def create_subfiles(output_dir, bed_files, csv_files):
    for bed_file in bed_files:
        bed_basename = os.path.basename(bed_file)
        bed_prefix = os.path.splitext(bed_basename)[0]

        bed_df = pd.read_csv(bed_file, sep='\t', dtype=str)

        for csv_file in csv_files:
            csv_basename = os.path.basename(csv_file)
            csv_prefix = os.path.splitext(csv_basename)[0]

            csv_df = pd.read_csv(csv_file, dtype=str)

            common_cells = set(bed_df['Cell']).intersection(csv_df['Cell'])

            filtered_bed_df = bed_df[bed_df['Cell'].isin(common_cells)]
            filtered_csv_df = csv_df[csv_df['Cell'].isin(common_cells)]

            output_bed_filename = os.path.join(output_dir, f"{bed_prefix}_{csv_prefix}.bed")
            output_csv_filename = os.path.join(output_dir, f"{bed_prefix}_{csv_prefix}.csv")

            # Reorder columns with 'Cell' column as the first column
            filtered_bed_df = filtered_bed_df[['Cell'] + [col for col in filtered_bed_df.columns if col != 'Cell']]
            filtered_csv_df = filtered_csv_df[['Cell'] + [col for col in filtered_csv_df.columns if col != 'Cell']]
            if 'Subpopulation' not in filtered_csv_df.columns:
                filtered_csv_df['Subpopulation'] = 0

            # Save BED file with header
            filtered_bed_df.to_csv(output_bed_filename, sep='\t', header=True, index=False)

            # Save CSV file
            filtered_csv_df.to_csv(output_csv_filename, index=False, na_rep="NA")

def main():
    parser = argparse.ArgumentParser(description="Filter bed and csv files based on common indices.")
    parser.add_argument("--output-dir", required=True, help="Output directory for sub-bed and sub-csv files")
    parser.add_argument("--bed-files", required=True, nargs='+', help="List of bed files")
    parser.add_argument("--csv-files", required=True, nargs='+', help="List of CSV files")

    args = parser.parse_args()

    output_dir = args.output_dir
    bed_files = args.bed_files
    csv_files = args.csv_files

    create_subfiles(output_dir, bed_files, csv_files)

if __name__ == "__main__":
    main()
