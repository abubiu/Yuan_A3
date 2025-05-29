import argparse
import pandas as pd
from pybedtools import BedTool
import tempfile
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Extract gene information based on genomic positions")
    parser.add_argument('input_file', type=str, help='Path to the input file')
    parser.add_argument('gtf_file', type=str, help='Path to the GTF annotation file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    parser.add_argument('--preferred_strand', type=str, choices=['+', '-'], default='+', help='Preferred strand to return')
    return parser.parse_args()

def extract_gene_info(gtf_field):
    try:
        match_gene_id = re.search(r'gene_id "([^"]+)"', gtf_field)
        gene_id = match_gene_id.group(1) if match_gene_id else None
        
        match_gene_biotype = re.search(r'gene_biotype "([^"]+)"', gtf_field)
        gene_biotype = match_gene_biotype.group(1) if match_gene_biotype else None


        return gene_id, gene_biotype
    except Exception as e:
        print(f"Error occurred: {e}")
        return None, None


def main():
    args = parse_args()

    df = pd.read_csv(args.input_file, sep='\t')

    gtf = BedTool(args.gtf_file)

    gene_info = []

    with tempfile.NamedTemporaryFile(mode='w+t', delete=False) as temp:
        for _, row in df.iterrows():
            chrom = row['chrom']  
            position = row['position'] 
            start = position - 1  
            end = position
            temp.write(f"{chrom}\t{start}\t{end}\n")
        temp.flush()  

        query = BedTool(temp.name)

    overlapping_features = query.intersect(gtf, wa=True, wb=True)

    gene_dict = {}
    geneid_dict = {}
    for feature in overlapping_features:
        chrom_pos = f"{feature.chrom}:{feature.start}-{feature.end}"
        gtf_field = feature.fields[11]
        feature_type=feature.fields[5] 
        strand = feature.fields[9]
        gene_id, gene_biotype = extract_gene_info(gtf_field)
        if chrom_pos not in gene_dict:
            gene_dict[chrom_pos] = {'geneid_dict': {}}
        if gene_id not in geneid_dict:
            geneid_dict[gene_id] = {'strand':[],'biotypes': [], "feature_type": []}
        geneid_dict[gene_id]["strand"].append(strand)
        geneid_dict[gene_id]["strand"]=list(set(geneid_dict[gene_id]["strand"]))
        geneid_dict[gene_id]["biotypes"].append(gene_biotype)
        geneid_dict[gene_id]["biotypes"]=list(set(geneid_dict[gene_id]["biotypes"]))
        geneid_dict[gene_id]["feature_type"].append(feature_type) 
        geneid_dict[gene_id]['feature_type'] = list(set(geneid_dict[gene_id]['feature_type']))
          
        gene_dict[chrom_pos]['geneid_dict'][gene_id]=geneid_dict[gene_id]
        
    info=[]
    for _, row in df.iterrows():
        chrom_pos = f"{row['chrom']}:{row['position']-1}-{row['position']}"
        gene_info_list = []
        region_info_list = [] 
        gene_data = gene_dict.get(chrom_pos, {'geneid_dict':{}})
        if gene_data=={'geneid_dict': {}}:
            gene_info_list.append("intergenic")
            region_info_list.append("intergenic region")
        else:
            biotypes = []
            feature_types = []
            region_types =[]
            gene_ids=list(gene_data['geneid_dict'].keys())
       
            strands = []

            for gene_id in gene_ids:
                gene_strands = gene_data['geneid_dict'][gene_id]['strand']
                strands.extend(gene_strands)  
            flat_strands=list(set(strands))
       
            if len(flat_strands) == 1:
                gene_ids=gene_ids
            else:
                preferred_genes = []
                for gene_id in gene_ids:
                    if args.preferred_strand in gene_data['geneid_dict'][gene_id]['strand']:
                        preferred_genes.append(gene_id)
                    else:
                        continue
                if len(preferred_genes)==0:
                    print ("No genes on the preferred strand")
                gene_ids=preferred_genes   
            for gene_id in gene_ids:
                gene_info = gene_data['geneid_dict'].get(gene_id, {'strand':["No_strand"],'biotypes': ["No_biotype"], "feature_type": ["No_feature_type"]})
                gene_biotype = gene_info['biotypes']
                gene_feature_type = gene_info['feature_type']  
                if "protein_coding" in gene_biotype:
                    if "CDS" in gene_feature_type:
                        region_type = "Coding RNA--exonic"
                    elif "CDS" not in gene_feature_type and "exon" in gene_feature_type:
                        region_type = "Coding RNA--untranslated region"
                    else:
                        region_type = "Coding RNA--intron"
                elif gene_biotype:
                    if "exon" in gene_feature_type:
                        region_type="Non-coding RNA--Exonic"
                    else:
                        region_type = "Non-coding RNA--intron"
                else:
                    region_type = " untranscribed intergenic"            
                gene_info_list.append(gene_id)
                region_info_list.append(region_type)
        info.append({
        "Gene": "; ".join(gene_info_list),  
         "Region": "; ".join(region_info_list)  
            })
    info_df = pd.DataFrame(info)
    print("------Work Done------")
    f = pd.concat([df, info_df], axis=1)

    f.to_csv(args.output_file, index=False, sep='\t')

if __name__ == "__main__":
    main()