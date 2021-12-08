#! /bin/bash

# downloads latest refseq assembly summary
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt -O assembly_summary_refseq.txt

# Uses $1 and $2 are separate parsing variables for organism names
# gets link from refseq db
link=$(grep $1 assembly_summary_refseq.txt | grep $2 | head -n 1 | cut -f 20)

# exits program if no matching organism is found
if [[ $link == '' ]]
then
    echo No matches found for $1 $2
    exit 1
fi

# generates G_species dir
g=$(expr substr $(echo $1) 1 1)
new_dir=$(echo $g)_$2
mkdir $new_dir

# gets file name
file=$(echo $link | awk -F'/' '{print $NF}')_genomic.fna.gz

# gets fna download link
fna=$(echo $link)/$(echo $file)

# downloads, moves, and unzips file
wget $fna
mv $file $new_dir
gzip -d $new_dir/$file