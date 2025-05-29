import pandas as pd
from Bio import SeqIO
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Extract neighboring bases from FASTA based on position and create base pairs")
    parser.add_argument("fasta_file", help="Path to the FASTA genome file")
    parser.add_argument("data_file", help="Path to the input file with chromosome and position information (CSV format)")
    parser.add_argument("output_file", help="Path to the output file (CSV format)")
    return parser.parse_args()

def extract_neighboring_bases(fasta_file, data_file, output_file):
    df = pd.read_csv(data_file, sep="\t") 
    
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_dict[record.id] = record.seq  
    
    pairing_info = []
    extended_pairing_info = []

    for index, row in df.iterrows():
        chrom = row['chrom']
        position = row['position']
        ref_base = row['ref']

        if chrom in seq_dict:
            seq = seq_dict[chrom]
            
            if 2 < position < len(seq) - 1: 
                if ref_base == 'C':
                    neighboring_base = seq[position - 2]  
                    pair = neighboring_base + ref_base            
  
                    extended_neighboring_base = seq[position - 3]  
                    extended_pair = extended_neighboring_base + neighboring_base + ref_base  
                elif ref_base == 'G':
                    neighboring_base = seq[position]  
                    pair = ref_base + neighboring_base  
                    
                    extended_neighboring_base = seq[position + 1]  
                    extended_pair = ref_base + neighboring_base + extended_neighboring_base  
                else:
                    print("******Error！Not C or G******")
                    pair = None
                    extended_pair = None
                pairing_info.append(pair)
                extended_pairing_info.append(extended_pair)
            else:
                print("Position too close to start or end of the sequence")
                pairing_info.append(None)  
                extended_pairing_info.append(None)  
        else:
            print("Chromosome not found in FASTA")
            pairing_info.append(None)  
            extended_pairing_info.append(None)
   
    df['Neighboring_Pair'] = pairing_info
    df['Extended_Neighboring_Pair'] = extended_pairing_info

    df.to_csv(output_file, sep="\t", index=False) 

    print(f"Completed. Output file saved as '{output_file}'")

if __name__ == "__main__":
    args = parse_args()
    extract_neighboring_bases(args.fasta_file, args.data_file, args.output_file)