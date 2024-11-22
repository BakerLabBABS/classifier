#!/bin/bash

echo "Please enter the relative path of the Jackhmmer directory: Example : /srv/scratch/baker/Jamiema/Testing_pipeline/example/Hits"
read -r DIR


# Check if the provided directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory $DIR does not exist."
    exit 1
fi


# Navigate to the specified directory
cd "$DIR" || exit

mkdir -p HMMER_OUT

### concatenate all HMMER output
cat *_seqHitsTable.txt > HMMER_OUT/All_HMMER_patric.out
cd HMMER_OUT || exit

## change the output to tab-delimited and solve the issues in some of the lines ###
sed -r 's/ {2,}/\t/g' All_HMMER_patric.out |grep -v "#"|sed 's/ -/\t-/g' > All_HMMER_tab_sep.out

### Filter at different evalues
awk -F "\t" '{if($5<1e-10){print}}' All_HMMER_tab_sep.out > All_HMMER_tab_sep_1e-10.out

### remove redundant hits ##
sort -t$'\t' -k1,1 -k6,6nr -gk5,5n All_HMMER_tab_sep_1e-10.out |sort -u -k1,1 > All_HMMER_sorted_1e-10.tsv

### Extract genes and species and remove redundant lines ###
awk -F '[\t[]' '{print $3"\t"$(NF)}' All_HMMER_sorted_1e-10.tsv |awk -F '[|\t]' '{print $1"\t"$(NF-1)}'|sort|uniq > tmp_1e-10.out

## Aggregate genes of same species together like key-value pairs ###
awk -F "\t" '{a[$2] = a[$2]"|"$1} END{for (i in a) print i "\t" a[i]}' tmp_1e-10.out |awk -F '[|]' '{print $0"\t"NF-1}' > species_genes_1e-10.tsv

### Remove intermediate files ##
rm  tmp_1e-10.out

### Generate input for presence-absence matrix ####
cut -f1-2 species_genes_1e-10.tsv |sed 's/\t|/\t/g' |sed 's/|/,/g' |sed  -i '1i species\tgenes' > species_genes_list.tsv
