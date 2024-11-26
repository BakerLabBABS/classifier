#!/bin/sh
  
name="flrC"
BASE="/srv/scratch/baker/Jamiema/New_analysis_set4/HMMER/PATRIC_hmmer/TOTAL_OUT/NewV1/subset_proteins/redo"
dir="/srv/scratch/baker/Jamiema/New_analysis_set5"
foldseek=/srv/scratch/baker/Jamiema/Structural_analysis/foldseek/bin/foldseek
 
cd "$dir"/"$name" || exit
mkdir -p pdb
mv output_"$name"/*rank_001*.pdb pdb/
 
cp -r "$BASE"/"$name"/"$name"_KEGG_alphafold "$dir"/"$name"/
 
### foldseek analysis ##
$foldseek easy-search "$name"_KEGG_alphafold pdb foldseek_out tmp --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,lddt,alntmscore
 
cat foldseek_out |awk -F "\t" '{if($13==1.000){print}}'|cut -f2|sort|uniq|cut -d '_' -f1-2
