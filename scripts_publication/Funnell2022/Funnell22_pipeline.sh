#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

#Funnell 2022 data
mkdir -p Funnell22/Metadata
cd Funnell22
# data from Zenodo

curl -s https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg19.chrom.sizes | sed 's/chr//g' > Metadata/hg19chr.sizes.tsv

datasetcount="$(ls signatures_dataset_folder/DLP/CNA/persample/*gz | wc -l)"
COUNTER=1
for sample in signatures_dataset_folder/DLP/CNA/persample/*gz ; do
	echo MnM on $sample : $COUNTER of $datasetcount
  trimmedname=$(basename "$sample" "_hscn.csv.gz")
	MnM -i $sample -o MnM/$trimmedname -n $trimmedname -g Metadata/hg19chr.sizes.tsv -w 500000 -s --Cellcol cell_id --CNcol state --sep , --seed 18671107
	COUNTER=$[$COUNTER +1]
done
