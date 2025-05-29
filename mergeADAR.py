import argparse
import pandas as pd

pd.set_option('display.max_columns', None)
def filter_and_merge(file_add_4to13, file_add_12or13, output_file):
    df_add_4to13 = [pd.read_table(file, sep="\t", header=0) for file in file_add_4to13]
    df_add_12or13 = [pd.read_table(file, sep="\t", header=0) for file in file_add_12or13]
    
    df_add_4to13_filtered = [df[(df['ref'] == 'A') | (df['ref'] == 'T')] for df in df_add_4to13]
    df_add_12or13_filtered = [df[(df['ref'] == 'A') | (df['ref'] == 'T')] for df in df_add_12or13]

    merged_df = df_add_4to13_filtered[0]
    for i, df in enumerate(df_add_4to13_filtered[1:], start=1):
        suffix = f"_file{i+1}"
        merged_df = pd.merge(merged_df, df, on=["chrom", "position", "ref"], suffixes=("", suffix))
    old = len(df_add_4to13_filtered)

    for i, df in enumerate(df_add_12or13_filtered):
        suffix = f"_file{i+1+old}"
        df_selected = df[["chrom", "position", "ref", "C_ref", "G_ref"]] 
        merged_with_df = pd.merge(merged_df, df_selected, on=["chrom", "position", "ref"], how="left", suffixes=("", f"{suffix}"))
    
        merged_with_df[f'file{i+1+old}_EVfre_column'] = merged_with_df.apply(
        lambda row: row[f'G_ref{suffix}'] if row['ref'] == 'A' else 
                    row[f'C_ref{suffix}'] if row['ref'] == 'T' else None, 
        axis=1
    )
    
        merged_df=merged_with_df

    merged_df.to_csv(output_file, sep="\t", index=False, header=True)

def main():
    parser = argparse.ArgumentParser(description="Filter and merge genomic data from multiple files based on A/T reference columns.")
    parser.add_argument('--add_4to13', nargs='+', required=True, help="Files that add columns 4-13")
    parser.add_argument('--add_12or13', nargs='+', required=True, help="Files that add columns 12 or 13")
    parser.add_argument('--output', required=True, help="Output file path")

    args = parser.parse_args()

    filter_and_merge(args.add_4to13, args.add_12or13, args.output)

if __name__ == "__main__":
    main()