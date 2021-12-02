#! /bin/bash

# gets link from refseq db
link=$(cat $1 | grep $2 | grep $3 | head -n 1 | cut -f 20)

# generates G_species dir
g=$(expr substr $(echo $2) 1 1)
new_dir=$(echo $g)_$3
mkdir $new_dir

# gets file name
file=$(echo $link | awk -F'/' '{print $NF}')_genomic.fna.gz

# gets fna download link
fna=$(echo $link)/$(echo $file)

# downloads, moves, and unzips file
wget $fna
mv $file $new_dir
gzip -d $new_dir/$file