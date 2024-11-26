#!/bin/bash
accession_list="${PWD}/$1"
output_directory="${PWD}/$2"
mkdir -p "$output_directory"

# Loop through each accession in the file and download the corresponding PDB file
while IFS= read -r accession; do
    # Construct the URL for downloading the PDB file
    pdb_url="https://alphafold.ebi.ac.uk/files/AF-${accession}-F1-model_v4.pdb"
    
    # Use wget to download the PDB file and save it to the output directory
    wget -q "$pdb_url" -P "$output_directory"
    
    # Check if the download was successful
    if [ $? -eq 0 ]; then
        echo "Downloaded $accession.pdb"
    else
        echo "Failed to download $accession.pdb"
    fi
done < "$accession_list"

echo "Download completed."
