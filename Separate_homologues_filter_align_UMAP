#!/bin/sh

echo "Enter the relative path of the directory containing Jackhmmer results: Example: /srv/scratch/baker/Jamiema/Testing_pipeline/example/Hits/HMMER_OUT"
read -r BASEDIR

# Check if the provided directory exists
if [ ! -d "$BASEDIR" ]; then
    echo "Directory $BASEDIR does not exist."
    exit 1
fi

echo "Enter the name of the flagellar protein to be analyzed containing : Example : fliG"
read -r name

# Navigate to the specified directory
cd "$BASEDIR" || exit

# Create a directory for the protein 
mkdir -p "$name"

cd "$name" || exit

# Query directory path - Replace the path with your path of query dir ##
QUERYDIR="/srv/scratch/baker/Jamiema/New_analysis_set4/HMMER/PATRIC_hmmer"

# Load modules
module load seqkit
module load seqtk
module load cdhit
module load mafft

# Extract unique homologs of specific protein
grep "$name|" "$BASEDIR"/All_HMMER_sorted_1e-10.tsv | \
    awk -F "\t" -v OFS="\t" '{print $1,$3,$5,$6,$NF}' | \
    sort -t$'\t' -k4,4nr -k3,3n | \
    sort -t$'\t' -k5,5 -u > "$name"_uniq_homologs.txt

# Extract only the IDs
cut -f1 "$name"_uniq_homologs.txt > "$name"_uniq_homologs.tsv

# Extract long IDs
grep -w -Ff "$name"_uniq_homologs.tsv /srv/scratch/baker/Jamiema/PATRIC_DB/GENOME_set*.faa | \
    cut -d '>' -f2 > "$name"_uniq_homologs1.tsv

# Remove temp files
rm "$name"_uniq_homologs.txt "$name"_uniq_homologs.tsv

# Extract fasta sequences
for i in /srv/scratch/baker/Jamiema/PATRIC_DB/GENOME_*.faa; do
    seqtk subseq "$i" "$name"_uniq_homologs1.tsv
done > "$name"_uniq_homologs1.fa

# Get the length of sequences and extract hits based on annotations
MAX_LENGTH=$(seqkit fx2tab --length --name --header-line "$name"_uniq_homologs1.fa | grep -i "$name" | sort -t$'\t' -k2,2n | cut -f2 | tail -1)
MIN_LENGTH=$(seqkit fx2tab --length --name --header-line "$name"_uniq_homologs1.fa | grep -i "$name" | sort -t$'\t' -k2,2nr | cut -f2 | tail -1)

echo "$MAX_LENGTH" > length_filter.txt
echo "$MIN_LENGTH" >> length_filter.txt

# Use that length and filter hits which have no annotations
seqkit fx2tab --length --name --header-line "$name"_uniq_homologs1.fa | \
    awk -v min_length="$MIN_LENGTH" -v max_length="$MAX_LENGTH" -F "\t" '{if(($2 >= min_length) && ($2 <= max_length)) {print $1}}' > length_filtered_hits.txt

seqtk subseq "$name"_uniq_homologs1.fa length_filtered_hits.txt > length_filtered_hits.fa

# Cluster at 90% identity
cd-hit -i length_filtered_hits.fa -o length_filtered_hits90.fa -c 0.9 -n 5

# Get the query sequences of the specific protein
grep -w "$name" "$QUERYDIR"/flagellar_proteins_renamed1.fa | sed 's/>//g' > kegg_"$name"
seqtk subseq "$QUERYDIR"/flagellar_proteins_renamed1.fa kegg_"$name" > kegg_"$name".fa

# Add the query sequences to the length-filtered and 90% clustered homologs
cat kegg_"$name".fa length_filtered_hits90.fa > "$name"_length_filtered_hits90.fa

# Remove intermediate files
rm kegg_"$name".fa

# Align the sequences
mafft "$name"_length_filtered_hits90.fa > "$name"_length_filtered_hits90.aln

# Run UMAP_DBSCAN script
mkdir -p UMAP_DBSCAN_output
python /srv/scratch/baker/Jamiema/Testing_pipeline/scripts/UMAP_DBSCAN_final.py "$name"_length_filtered_hits90.aln kegg_"$name" UMAP_DBSCAN_output
