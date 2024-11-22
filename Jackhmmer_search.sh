#!/bin/bash

### Load modules ##
module load hmmer/3.2.1

echo "Please enter the relative path of the directory containing fasta files: Example : /srv/scratch/baker/Jamiema/Testing_pipeline/example"
read -r DIR

echo "Please enter the filename with list of query fasta file"
read -r FILE

# Check if the provided directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory $DIR does not exist."
    exit 1
fi

# Check if the provided file exists
if [ ! -f "$FILE" ]; then
    echo "File $FILE does not exist"
    exit 1
fi


# Navigate to the specified directory
cd "$DIR" || exit

### create directories for storing output files ##
mkdir -p Alignment Hits output logs

### Run Jackhmmer on fasta files ###
for filename in *.faa; do
    # Extract the full species name from the filename
    species=$(basename "$filename" .PATRIC.faa)
    
    # Ensure the species-specific fasta file exists
    if [ ! -f "$species.PATRIC.faa" ]; then
        echo "Species fasta file $species.PATRIC.faa does not exist."
        continue
    fi
    
    nohup jackhmmer --noali -A Alignment/"$species"_aln  --tblout Hits/"$species"_seqHitsTable.txt -o output/"$species"_out "$FILE"  "$species".PATRIC.faa > logs/"$species"_log.txt 2>&1 &
done

echo "Jackhmmer analysis started with nohup. Check the logs directory for progress."
