import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def load_genome(fasta_file):
    genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome[record.id] = record.seq
    return genome

def get_flanking_sequence(chrom, position, genome, flanking_length, strand):
    seq = genome.get(chrom)
    if seq:
        start = max(0, position - flanking_length - 1)
        end = min(len(seq), position + flanking_length)
        extracted_seq = str(seq[start:position]) + str(seq[position:end])
        print(f"Extracted sequence: {extracted_seq}")
        print(f"Strand parameter: '{strand}'")
        
        if strand == '-':
            print("Performing reverse complement")
            extracted_seq = str(Seq(extracted_seq).reverse_complement())
            print(f"Reverse complemented sequence: {extracted_seq}")
        
        return extracted_seq
    return ""

def main(args):
    genome = load_genome(args.genome_fasta)

    df = pd.read_csv(args.input_file, sep="\t") 
    
    flanking_sequences = []
    for _, row in df.iterrows():
        chrom = row['chrom']
        position = row['position']
        
        strand = row.get('strand', args.strand) 
        
        flanking_seq = get_flanking_sequence(chrom, position, genome, args.flanking_length, strand)
        flanking_sequences.append(flanking_seq)

    df['flanking_sequence'] = flanking_sequences

    df.to_csv(args.output_file, index=False, sep="\t")  
    print(f"结果已保存到 {args.output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch extract 6bp flanking sequences around genomic positions")
    
    parser.add_argument("genome_fasta", help="Path to the genome FASTA file")
    parser.add_argument("input_file", help="Path to the input file (CSV/TSV format)")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("--flanking_length", type=int, default=6, help="Number of bases to extract upstream and downstream (default: 6)")
    parser.add_argument("--strand", choices=['+', '-'], default='+', help="Strand direction to use (default: +")
    
    args = parser.parse_args()
    main(args)