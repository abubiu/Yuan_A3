import pandas as pd
import argparse

def percent_to_float(value):
    return float(value.rstrip('%'))

def filter_data(input_file, output_file):
    data = pd.read_csv(input_file, sep="\t")
    
    filtered_data = data[
        (data["G"] >= 2) &
        (data["G_ref"].apply(percent_to_float) >= 3) &
        (data["G_file2"] >= 2) &
        (data["G_ref_file2"].apply(percent_to_float) >= 3) &
        (data["file3_EVfre_column"].apply(percent_to_float) == 0) &
        (data["file4_EVfre_column"].apply(percent_to_float) == 0)
    ]
    
    filtered_data.to_csv(output_file, sep="\t", index=False)
    print(f"filter finished，save to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="filter your data")
    parser.add_argument("-i", "--input", required=True, help="inputfile")
    parser.add_argument("-o", "--output", required=True, help="outfile")
    args = parser.parse_args()

    filter_data(args.input, args.output)

if __name__ == "__main__":
    main()