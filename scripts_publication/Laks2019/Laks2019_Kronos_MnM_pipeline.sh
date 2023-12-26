#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

###############################################################
################# DEFINE FILE/DIRECTORY PATHS #################
###############################################################

# Base Files
OutputDir=Laks2019/Kronos
# Software
TrimGalore=$(which trim_galore)
Cutadapt=$(which cutadapt)
Picard=#path/picard.jar
Java=$(which java)
# Genome files
hg38Ref=#path/WholeGenomeFasta/genome.fa
hg38Blist=#path/hg38-blacklist.v2.bed
hg38BT2=#path/Bowtie2Index/genome

cores=8

mkdir $OutputDir

cd Laks2019

###########
# Binning #
###########
#25k
Kronos binning -R $hg38Ref -c $cores -o ../Binned_Genomes/hg38/25kPE -i $hg38BT2 --bin_size 25000 --paired_ends -B $hg38Blist --chr_range 1:22,X,Y

#######
# CNV #
#######
cd Laks2019
datasetcount="$(ls ./BAMs_hg38 | wc -l)"
cat Metadata/Laks2019_cell_descriptions_with_ploidy.csv | awk -F',' '{ if (system(" test -f Kronos/CNV/"$4"_"$5"_"$6"_cnv_calls.bed ")) { if ($10 == "M") {chr = "X,Y"} else {chr = "X"} {if ($9 == "Homo_sapiens") {system(" date ; echo "$4"_"$5"_"$1": file "NR" of '$datasetcount' ; echo "$2" cells to process. ; Kronos CNV -D BAMs_hg38/"$4"_"$5"_"$1" -B ../Binned_Genomes/hg38/25kPE/*.tsv -c '$cores' -o Kronos/CNV -m "$11" -M "$12" -e "$4"_"$5"_"$6" --chr_range 1:22,"chr)}} } else {print $4,$5,$6,"DONE..Exiting..."} }'

## Cell stats
sh calculate_stats.sh > Metadata/Laks_read_stats.csv

##############
# Diagnostic # For MnM
##############
# remove GM18507_SA928_A90648B and GM18507_SA928_A90685
cd Laks2019
mkdir Kronos/Diagnostic
datasetcount="$(ls ./BAMs_hg38 | wc -l)"
tail -n+2 Metadata/Laks2019_cell_descriptions_with_ploidy.csv | grep -v A90648B | grep -v A90685 | awk -F','  '{
	if ($9 == "Homo_sapiens" && length($13) != 0 && $7=="NA") {
		system(" date ; echo "$4"_"$5"_"$6": dataset "NR" of '$datasetcount' ; Kronos diagnostic -F Kronos/CNV/"$4"_"$5"_"$6"_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Dev/"$4"_"$5" -b "$4"_"$5"_"$6" -m "$13" -d ")
	} else if ($7!="NA") {
        system(" date ; echo "$4"_"$5"_"$6": dataset "NR" of '$datasetcount' ; Kronos diagnostic -F Kronos/CNV/"$4"_"$5"_"$6"_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Dev/"$4"_"$7" -b "$4"_"$5"_"$6" -m "$13" -d ")
    } }'

#################
##### M n M #####
#################
mkdir Kronos/CNV/merged
# Metadata (cell name, cell type)
for bed in Kronos/CNV/*_per_Cell_summary_metrics.csv ; do
	name=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' )
	tail -n+2 Kronos/CNV/"$name"*_per_Cell_summary_metrics.csv | awk -F',' -v n="$name" 'OFS="\t"{print $1,n}'
done > Metadata/All_laks_CellsAndTypes.tsv
# FILTER CELLS
tail -n+2 Metadata/Laks2019_cell_descriptions_with_ploidy.csv | awk -F',' '{if ($9 == "Homo_sapiens" && length($13) != 0) { system(" echo "$4"_"$5" ;  head -n1 Kronos/CNV/"$4"_"$5"_"$6"_cnv_calls.bed > Kronos/CNV/merged/"$4"_"$5".bed ") } }'
for bed in $(ls Kronos/CNV/*_cnv_calls.bed | grep -v tmp ) ; do
	celltype=$(echo $bed | awk -F'_' '{print $1 }' | sed 's/Kronos\/CNV\///' )
	name=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_cnv_calls.bed//' | awk -F'_A' '{print $1 }')
	fullname=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_cnv_calls.bed//')
	echo $celltype $name $fullname
	awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 Kronos/diagnostic/Dev/$celltype*/"$fullname"*_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') Kronos/CNV/"$fullname"*_cnv_calls.bed >> Kronos/CNV/merged/"$name".bed
done

## MnM ##
for bed in $(ls Kronos/CNV/merged/* ) ; do
	celltype=$(echo $bed | sed 's/Kronos\/CNV\/merged\///' | sed 's/.bed//' )
	echo $celltype
	MnM -i $bed -o MnM/$celltype -n $celltype -r -s -b --seed 18671107 --groups Metadata/All_laks_CellsAndTypes.tsv
done
# fix the ones that failed (probably technical failure:
MnM -i Kronos/CNV/merged/184-hTERT_SA906.bed -o MnM/184-hTERT_SA906 -n 184-hTERT_SA906 -r -b --seed 18671107 --groups Metadata/All_laks_CellsAndTypes.tsv
MnM -i Kronos/CNV/merged/ERpos-PDX_SA995X5XB01910.bed -o MnM/ERpos-PDX_SA995X5XB01910 -n ERpos-PDX_SA995X5XB01910 -r -b --seed 18671107 --groups Metadata/All_laks_CellsAndTypes.tsv

##############
# Diagnostic # For Kronos this time
##############
# WhoIsWho
for csv in Kronos/CNV/*_per_Cell_summary_metrics.csv ; do
	name=$(echo $csv | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' | awk -F'_A' '{print $1 }')
	echo $name
	Kronos WhoIsWho -F $csv -o Kronos/Diagnostic/Whoiswho -W MnM/$name/*_metadata.tsv
done

# Diagnostic for Kronos
tail -n+2 Metadata/Laks2019_cell_descriptions_with_ploidy.csv | awk -F',' '{if ($9 == "Homo_sapiens" && length($13) != 0) {system(" date ; echo "$4"_"$5"_"$6" ; Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_"$4"_"$5"_"$6"*.csv -o Kronos/Diagnostic  -b "$4"_"$5"_"$6" -m "$13" -C ") } }'

##############
# Prepare RT #
##############

# Prepare all subpopulations/subfiles
tail -n+2 Metadata/Laks2019_cell_descriptions_with_ploidy.csv | awk -F',' '{if ($9 == "Homo_sapiens") {system(" mkdir -p Kronos/subpopulation_files_split/"$4"_"$5"") } }'
for dataset in MnM/* ; do
	name=$(echo $dataset | sed 's/MnM\///' )
	echo $name
	python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/$name --csv-files $(ls Kronos/Diagnostic/Whoiswho/phased_"$name"*_per_Cell_summary_metrics.csv) --bed-files $(ls $dataset/*.BED.gz)
done

#### remove failed files (technical failure)
rm Kronos/subpopulation_files_split/GM18507_SA928/*SA928_A90685*
rm Kronos/subpopulation_files_split/GM18507_SA928/*SA928_A90648B*
rm Kronos/subpopulation_files_split/GM18507_SA928/*SA928_A96226B*
#failed RT of 1 GM cell sample. Diagnostic file should exist for RT, regardless of the parameters (S phase not concerned here). Copying it here for this reason
cp Kronos/Diagnostic/GM18507_SA928_A96172B_settings.txt Kronos/Diagnostic/GM18507_SA928_A90648B_settings.txt

######
# RT #
######
cores=4
for dataset in $(ls MnM ) ; do
	name=$(echo $dataset | sed 's/MnM\///' )
	echo $name
	data=$(ls Kronos/subpopulation_files_split/$name/*.csv | sort | awk -F'_phased_' '{print $2}' | sed 's/_per_Cell_summary_metrics.csv//g')
	cell_type=$(echo $data | awk -F'_' '{print $1}')
	sample_id=$(echo $data | awk -F'_' '{print $2}')
	library_id=$(echo $data | awk -F'_' '{print $3}')
	subpop=$(ls Kronos/subpopulation_files_split/$name/*.csv | sort | awk -F'.BED' '{print $1}' | sed 's/Kronos\/subpopulation_files_split\/'$name'\/'$name'//g' | sed 's/_/S/g')

	Kronos RT \
	-F $(ls Kronos/subpopulation_files_split/$name/*.csv | sort | tr '\n' ',' | sed 's/,$//') \
	-T $(ls Kronos/subpopulation_files_split/$name/*.bed | sort | tr '\n' ',' | sed 's/,$//') \
	-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
	-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
	-o Kronos/RT/$name \
	-b $(paste -d '_' <(echo "$library_id") <(echo "$subpop") | tr '\n' ',' | sed 's/,$//' ) \
	-f 200kb \
	-g $(paste -d '_' <(echo "$cell_type") <(echo "$sample_id") <(echo "$subpop") | tr '\n' ',' | sed 's/,$//' ) \
	-S $(paste -d '\0' <(yes "Kronos/Diagnostic/" | head -n $(echo $data | wc -l)) <(echo "$data") <(yes "_settings.txt" | head -n $(echo $data | wc -l))  | tr '\n' ',' | sed 's/,$//') \
	-B 200000 \
	-c $cores \
	-N 5 \
	-p \
	--extract_G1_G2_cells \
	--chr_range 1:22,X
done

#### DRED #####
for dataset in $(ls MnM ) ; do
	name=$(echo $dataset | sed 's/MnM\///' )
	echo $name
	Kronos DRed -C Kronos/RT/$name/200kb_single_cells_CNV_200Kb.tsv -o Kronos/DRed/$name -f $name -U --CNV_values B -c 6 -s 18671107
done
