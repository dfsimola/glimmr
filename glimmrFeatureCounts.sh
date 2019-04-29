#!/bin/bash

cd /home/dsimola/harpybrd/Tcell_IBD.DnaseI.data.23jan19.glimmr/mapped_reads_best/;

INDEX="/home/dsimola/harpybrd/Altius_DHS_index.22jan19.SAF";

mkdir ../featureCounts/

module load subread

for i in $(ls *unique.map); do 
	# echo ../featureCounts/featureCounts.$i.txt;
	
	# original
	# featureCounts -p -T 10 -F SAF -a $INDEX -o ../featureCounts/featureCounts.$i.txt $i;
	
	featureCounts -p -s 0 -d 0 -D 1000 -T 10 -F SAF -a $INDEX -o ../featureCounts.d0-1000/featureCounts.$i.txt $i;
	
done


