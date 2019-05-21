#!/bin/bash


# https://deweylab.github.io/RSEM/rsem-calculate-expression.html

# indir with sam files
cd $1;

module load perl RSEM

for i in $(ls *unique.map); do 
	
	echo convert-sam-for-rsem -p 10 "$i" "$i.rsem";
	convert-sam-for-rsem -p 10 "$i" "$i.rsem";
	
done


