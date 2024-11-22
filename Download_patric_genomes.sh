#!/bin/bash

echo "Please enter the relative path of the directory: Example : /srv/scratch/baker/Jamiema/pipeline/example"
read -r DIR

echo "Please enter the filename with list of IDS for download:"
read -r FILE

echo "ID file format is a list of ids in separate lines:
    325777.4
    325775.3
    637910.3"

# Check if the provided directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory $DIR does not exist."
    exit 1
fi

# Check if the provided file exists
if [ ! -f "$DIR/$FILE" ]; then
    echo "File $FILE does not exist in directory $DIR."
    exit 1
fi

# Navigate to the specified directory
cd "$DIR" || exit


# Read each line (ID) from the file and use it in the wget command
while IFS= read -r LINE; do
    wget -N "ftp://ftp.bvbrc.org/genomes/$LINE/$LINE.PATRIC.faa"
done < "$FILE"
