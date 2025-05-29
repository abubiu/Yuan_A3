####################
#Part1
#1.mapping
hisat2 -p 10 ~/database/GCF_000001405.40_GRCh38.p14_index -1 A3_R1.fq.gz -2 A3_R2.fq.gz -S A3_rep1_hisat2_aln.sam --rna-strandness RF 

#2.samtools
samtools view -Sb A3_rep1_hisat2_aln.sam >A3_rep1_hisat2_aln.bam && samtools sort -@ 5 -O bam -o A3_rep1_hisat2_aln_sorted.bam A3_rep1_hisat2_aln.bam  

#3.samtools
samtools addreplacerg -r "@RG\tID:RG1\tSM:A3\tPL:Illumina\tLB:Library.fa" -o A3_rep1_hisat2_aln_sorted_tag.bam A3_rep1_hisat2_aln_sorted.bam

#4.duplicate
time java -jar ~/miniconda3/envs/backup/share/picard-3.2.0-0/picard.jar MarkDuplicates -I A3_rep1_hisat2_aln_sorted_tag.bam -O A3_rep1_hisat2_aln_sorted_tag_markdup.bam -M A3_rep1_hisat2_aln_sorted_tag_markdup_metrix.txt -REMOVE_DUPLICATES true

#5.clipOverlap
bam clipOverlap --in A3_rep1_hisat2_aln_sorted_tag_markdup.bam --out A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip.bam --poolSize 10000000

#6.plus or minus
samtools view A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip.bam | egrep "XS:A:-" >A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_minus.sam
samtools view A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip.bam | egrep "XS:A:\+" >A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_plus.sam

#7.sam to bam
samtools view -bt ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna.fai -o A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_minus.bam A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_minus.sam
samtools view -bt ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna.fai -o A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_plus.bam A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_plus.sam

#8.mpileup.txt generation
samtools mpileup -f ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna -AB A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_minus.bam -d 10000 -q 1 -Q 20 >A3_rep1_hisat2_minus_mpileup.txt
samtools mpileup -f ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna -AB A3_rep1_hisat2_aln_sorted_tag_markdup_bamclip_plus.bam -d 10000 -q 1 -Q 20 >A3_rep1_hisat2_plus_mpileup.txt

#9.filter sites with total < 6
split -l 1000000 A3_rep1_hisat2_minus_mpileup.txt part_ && ls part_* | parallel -j 36 'awk "\$4 < 6 {print \$1, \$2}" {} > {}.out' && cat part_*.out > A3_rep1_hisat2_minus_mpileup_depth6_ID.txt && rm part_*
split -l 1000000 A3_rep1_hisat2_plus_mpileup.txt part_ && ls part_* | parallel -j 36 'awk "\$4 < 6 {print \$1, \$2}" {} > {}.out' && cat part_*.out > A3_rep1_hisat2_plus_mpileup_depth6_ID.txt && rm part_*

#10.merge and remove duplicate sites between EV and A3s.
#A3 minus
cat 1st_editingdata/A3/A3_rep1_hisat2_minus_mpileup_depth6_ID.txt 2nd_editingdata/A3/A3_rep2_hisat2_minus_mpileup_depth6_ID.txt >A3_rep1rep2_hisat2_minus_mpileup_depth6_ID.txt
awk '!seen[$0]++' A3_rep1rep2_hisat2_minus_mpileup_depth6_ID.txt >A3_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt
#A3 plus
cat 1st_editingdata/A3/A3_rep1_hisat2_plus_mpileup_depth6_ID.txt 2nd_editingdata/A3/A3_rep2_hisat2_plus_mpileup_depth6_ID.txt >A3_rep1rep2_hisat2_plus_mpileup_depth6_ID.txt
awk '!seen[$0]++' A3_rep1rep2_hisat2_plus_mpileup_depth6_ID.txt >A3_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt
#merge minus
cat EV_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt A3_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt >A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_ID.txt
awk '!seen[$0]++ ' A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_ID.txt >A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt
#merge plus
cat EV_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt A3_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt >A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_ID.txt
awk '!seen[$0]++ ' A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_ID.txt >A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt

####################################################
#Part2
#filtering
#1. >=6 total calls in every sample(including N)
#minus
#A3 two duplicate
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt 1st_editingdata/A3/A3_rep1_hisat2_minus_mpileup.txt >filter/EV_HsA3A/1st_A3_rep1_hisat2_minus_mpileup.txt
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt 2nd_editingdata/A3/A3_rep2_hisat2_minus_mpileup.txt >filter/EV_HsA3A/2nd_A3_rep2_hisat2_minus_mpileup.txt
#EV two duplicate
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt 1st_editingdata/EV/EV_rep1_hisat2_minus_mpileup.txt >filter/EV_HsA3A/1st_EV_rep1_hisat2_minus_mpileup.txt
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_minus_mpileup_depth6_uniqueID.txt 2nd_editingdata/EV/EV_rep2_hisat2_minus_mpileup.txt >filter/EV_HsA3A/2nd_EV_rep2_hisat2_minus_mpileup.txt
#plus
#A3 two duplicate
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt 1st_editingdata/A3/A3_rep1_hisat2_plus_mpileup.txt >filter/EV_HsA3A/1st_A3_rep1_hisat2_plus_mpileup.txt
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt 2nd_editingdata/A3/A3_rep2_hisat2_plus_mpileup.txt >filter/EV_HsA3A/2nd_A3_rep2_hisat2_plus_mpileup.txt
#EV two duplicate
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt 1st_editingdata/EV/EV_rep1_hisat2_plus_mpileup.txt >filter/EV_HsA3A/1st_EV_rep1_hisat2_plus_mpileup.txt
awk 'NR==FNR{a[$1" "$2]; next} !($1" "$2 in a)' A3+EV_rep1rep2_hisat2_plus_mpileup_depth6_uniqueID.txt 2nd_editingdata/EV/EV_rep2_hisat2_plus_mpileup.txt >filter/EV_HsA3A/2nd_EV_rep2_hisat2_plus_mpileup.txt

#2. mpileup to baserequency.txt
python basefrequency.py 1st_EV_rep1_hisat2_minus_mpileup.txt 1st_EV_rep1_hisat2_minus_base.txt

#3.
#total A/T/C/G base call >6 in every sample (no N)
#minus
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 1st_EV_rep1_hisat2_minus_base.txt > 1st_EV_rep1_hisat2_minus_basedepth6.txt
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 1st_A3_rep1_hisat2_minus_base.txt >1st_A3_rep1_hisat2_minus_basedepth6.txt
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 2nd_EV_rep2_hisat2_minus_base.txt > 2nd_EV_rep2_hisat2_minus_basedepth6.txt
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 2nd_A3_rep1_hisat2_minus_base.txt >2nd_A3_rep2_hisat2_minus_basedepth6.txt
#plus
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 1st_EV_rep1_hisat2_plus_base.txt >1st_EV_rep1_hisat2_plus_basedepth6.txt
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 1st_A3_rep1_hisat2_plus_base.txt >1st_A3_rep1_hisat2_plus_basedepth6.txt
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 2nd_EV_rep2_hisat2_plus_base.txt >2nd_EV_rep2_hisat2_plus_basedepth6.txt
awk 'BEGIN {FS=OFS="\t"} NR == 1 || $4 >= 6' 2nd_A3_rep2_hisat2_plus_base.txt >2nd_A3_rep2_hisat2_pluss_basedepth6.txt

#4.：
# filtering >6 position between EV and A3
# minus
awk -F'\t' '{print $1, $2}' 1st_EV_rep1_hisat2_minus_basedepth6.txt > Mposition1.txt  
awk -F'\t' '{print $1, $2}' 1st_A3_rep1_hisat2_minus_basedepth6.txt> Mposition2.txt
awk -F'\t' '{print $1, $2}' 2nd_EV_rep2_hisat2_minus_basedepth6.txt> Mposition3.txt
awk -F'\t' '{print $1, $2}' 2nd_A3_rep2_hisat2_minus_basedepth6.txt > Mposition4.txt
#intersection position
python find_commonline.py Mposition1.txt Mposition2.txt Mposition3.txt Mposition4.txt minus_common_positions.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' minus_common_positions.txt 1st_EV_rep1_hisat2_minus_basedepth6.txt> 1st_EV_rep1_hisat2_minus_basedepth6common.txt
aawk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' minus_common_positions.txt 2nd_EV_rep2_hisat2_minus_basedepth6.txt >2nd_EV_rep2_hisat2_minus_basedepth6common.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' minus_common_positions.txt 1st_A3_rep1_hisat2_minus_basedepth6.txt >1st_A3_rep1_hisat2_minus_basedepth6common.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' minus_common_positions.txt 2nd_A3_rep2_hisat2_minus_basedepth6.txt >2nd_A3_rep2_hisat2_minus_basedepth6common.txt
#########################################

# plus
awk -F'\t' '{print $1, $2}' 1st_EV_rep1_hisat2_plus_basedepth6.txt > Pposition1.txt
awk -F'\t' '{print $1, $2}' 1st_A3_rep1_hisat2_plus_basedepth6.txt> Pposition2.txt
awk -F'\t' '{print $1, $2}' 2nd_EV_rep2_hisat2_plus_basedepth6.txt > Pposition3.txt
awk -F'\t' '{print $1, $2}' 2nd_A3_rep2_hisat2_plus_basedepth6.txt > Pposition4.txt
#intersection position
python find_commonline.py Pposition1.txt Pposition2.txt Pposition3.txt Pposition4.txt plus_common_position.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' plus_common_183181450positions.txt 1st_EV_rep1_hisat2_plus_basedepth6.txt> 1st_EV_rep1_hisat2_plus_basedepth6common.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' plus_common_183181450positions.txt 2nd_EV_rep2_hisat2_plus_basedepth6.txt >2nd_EV_rep2_hisat2_plus_basedepth6common.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' plus_common_183181450positions.txt 1st_A3_rep1_hisat2_plus_basedepth6.txt >1st_A3_rep1_hisat2_plus_basedepth6common.txt
awk -F'[[:space:]]+' 'NR==FNR {positions[$1, $2]; next} NR==1 {print; next} ($1, $2) in positions' plus_common_183181450positions.txt 2nd_A3_rep2_hisat2_plus_basedepth6.txt >2nd_A3_rep2_hisat2_plus_basedepth6common.txt

#5.
# filter C or G
#minus
python merge_APOBEC.py --add_4to13 1st_A3_rep1_hisat2_minus_basedepth6common.txt 2nd_A3_rep2_hisat2_minus_basedepth6common.txt --add_10or11 1st_EV_rep1_hisat2_minus_basedepth6common.txt 2nd_EV_rep2_hisat2_minus_basedepth6common.txt --output EV+A3_hisat2_minus_depth6_baseoutdepth6common_filter.txt
#plus
python merge_APOBEC.py --add_4to13 1st_A3_rep1_hisat2_plus_basedepth6common.txt 2nd_A3_rep2_hisat2_plus_basedepth6common.txt --add_10or11 1st_EV_rep1_hisat2_plus_basedepth6common.txt 2nd_EV_rep2_hisat2_plus_basedepth6common.txt --output EV+A3_hisat2_plus_depths6_baseoutdepth6common_filter.txt

#ADAR
#minus
python merge_ADAR.py --add_4to13 ADAR_rep1_hisat2_minus_basedepth6common.txt ADAR_rep2_hisat2_minus_basedepth6common.txt --add_12or13 ADAR_rep1_hisat2_minus_basedepth6common.txt ADAR_rep2_hisat2_minus_basedepth6common.txt --output EV+ADAR_hisat2_minus_depth6_baseoutdepth6common_filter.txt
#plus
python merge_ADAR.py --add_4to13 ADAR_rep1_hisat2_plus_basedepth6common.txt ADAR_rep2_hisat2_plus_basedepth6common.txt --add_12or13 ADAR_rep1_hisat2_plus_basedepth6common.txt ADAR_rep2_hisat2_plus_basedepth6common.txt --output EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filter.txt


#6. filter ref C or G
#A3
#minus
awk 'NR == 1 || $3 == "C"' EV+A3_hisat2_minus_depth6_baseoutdepth6common_filter.txt >EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterC.txt
awk 'NR == 1 || $3 == "G"' EV+A3_hisat2_minus_depth6_baseoutdepth6common_filter.txt >EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterG.txt
#plus
awk 'NR == 1 || $3 == "C"' EV+A3_hisat2_plus_depth6_baseoutdepth6common_filter.txt >EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterC.txt
awk 'NR == 1 || $3 == "G"' EV+A3_hisat2_plus_depth6_baseoutdepth6common_filter.txt >EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterG.txt

#ADAR
#minus
awk 'NR == 1 || $3 == "A"' EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filter.txt >EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filterA.txt
awk 'NR == 1 || $3 == "T"' EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filter.txt >EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filterT.txt
#plus
awk 'NR == 1 || $3 == "A"' EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filter.txt >EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterA.txt
awk 'NR == 1 || $3 == "T"' EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filter.txt >EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterT.txt

#7. fitler site with level >3%
#A3
#minus
python filterC.py -i EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterC.txt -o EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterC_2or3%.txt
python filterG.py -i EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterG.txt -o EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterG_2or3%.txt
#plus
python filterC.py -i EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterC.txt -o EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterC_2or3%.txt
python filterG.py -i EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterG.txt -o EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterG_2or3%.txt

#ADAR
#minus
python filterA.py -i EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filterA.txt -o EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filterA_2or3%.txt
python filterT.py -i EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filterT.txt -o EV+ADAR_hisat2_minus_depths6_baseoutdepth6common_filterT_2or3%.txt
#plus
python filterA.py -i EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterA.txt -o EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterA_2or3%.txt
python filterT.py -i EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterT.txt -o EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterT_2or3%.txt

#8. mutation preference
#APOBEC
python basepreferenceAPOBEC3.py ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna EV+A3_hisat2_plus_depth6_baseoutdepth6common_filterC_2or3%.txt EV+A3_hisat2_plus_depth6_baseoutdepth6common_filterC_2or3%_prference3.txt
#ADAR
python basepreference3ADAR.py ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterA_2or3%.txt EV+ADAR_hisat2_plus_depths6_baseoutdepth6common_filterA_2or3%_preference3.txt

#9. gene name generation
#plus
python gene.py EV+A3_hisat2_plus_depths6_baseoutdepth6common_filterC_2or3%_preference3.txt ~/database/GCF_000001405.40_GRCh38.p14_genomic_notitle.gtf A3plus_results.txt --preferred_strand +
#minus
python gene.py EV+A3_hisat2_minus_depths6_baseoutdepth6common_filterC_2or3%_preference3.txt ~/database/GCF_000001405.40_GRCh38.p14_genomic_notitle.gtf A3_minus_results.txt --preferred_strand -

#10. flanking sequence
#flanking mutation site upstream and downstream 6bp total 13
#minus
python positionmotif.py ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna A3_minus_results.txt A3S_minus_motif6_sequence.txt --strand - --flanking_length 6
#plus
python positionmotif.py ~/database/GCF_000001405.40_GRCh38.p14_genomic.fna A3_plus_results.txt A3_plus_motif6_sequence.txt --strand + --flanking_length 6