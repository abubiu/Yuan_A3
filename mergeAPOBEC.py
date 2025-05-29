import argparse
import pandas as pd

pd.set_option('display.max_columns', None)
def filter_and_merge(file_add_4to13, file_add_10or11, output_file):
    df_add_4to13 = [pd.read_table(file, sep="\t", header=0) for file in file_add_4to13]
    df_add_10or11 = [pd.read_table(file, sep="\t", header=0) for file in file_add_10or11]
    
    df_add_4to13_filtered = [df[(df['ref'] == 'C') | (df['ref'] == 'G')] for df in df_add_4to13]
    df_add_10or11_filtered = [df[(df['ref'] == 'C') | (df['ref'] == 'G')] for df in df_add_10or11]

    merged_df = df_add_4to13_filtered[0]
    for i, df in enumerate(df_add_4to13_filtered[1:], start=1):
        suffix = f"_file{i+1}"
        merged_df = pd.merge(merged_df, df, on=["chrom", "position", "ref"], suffixes=("", suffix))
    old = len(df_add_4to13_filtered)

    for i, df in enumerate(df_add_10or11_filtered):
        suffix = f"_file{i+1+old}"
        df_selected = df[["chrom", "position", "ref", "A_fre", "T_ref"]]  
        merged_with_df = pd.merge(merged_df, df_selected, on=["chrom", "position", "ref"], how="left", suffixes=("", f"{suffix}"))
    
        merged_with_df[f'file{i+1+old}_EVfre_column'] = merged_with_df.apply(
        lambda row: row[f'T_ref{suffix}'] if row['ref'] == 'C' else 
                    row[f'A_fre{suffix}'] if row['ref'] == 'G' else None, 
        axis=1
    )
    
        merged_df=merged_with_df
    merged_df.to_csv(output_file, sep="\t", index=False, header=True)

def main():
    parser = argparse.ArgumentParser(description="Filter and merge genomic data from multiple files based on C/G reference columns.")
    parser.add_argument('--add_4to13', nargs='+', required=True, help="Files that add columns 4-13")
    parser.add_argument('--add_10or11', nargs='+', required=True, help="Files that add columns 10 or 11")
    parser.add_argument('--output', required=True, help="Output file path")
    args = parser.parse_args()

    filter_and_merge(args.add_4to13, args.add_10or11, args.output)

if __name__ == "__main__":
    main()