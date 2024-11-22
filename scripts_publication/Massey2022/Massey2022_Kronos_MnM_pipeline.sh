#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

#########
# PATHS #
#########
# Kronos
cores=6

TrimGalore=$(which trim_galore)
Cutadapt=$(which cutadapt)
Picard=#path/picard.jar
Java=$(which java)

hg38BT2=#path/Bowtie2Index/genome
hg38Ref=#path/WholeGenomeFasta/genome.fa
hg38Blist=#path/hg38-blacklist.v2.bed

BAMDir=Massey2022/BAM
RefRTs=Massey2022/Kronos/Reference_RT

mkdir -p $RefRTs
cd Massey2022/Kronos

# rename BAM directories originating from BAM files (= Not from SRA)
cd Massey2022/BAM
mv GM12878_G1_library1_replicate1 SRR16352033_GM12878-G1
mv GM12878_G2_library1_replicate1 SRR16352032_GM12878-G2
mv GM12878_S_library1_replicate1 SRR16352030_GM12878-S
mv GM12878_earlyS_library1_replicate1 SRR16352031_GM12878-earlyS
mv GM12878_lateS_library1_replicate1 SRR16352029_GM12878-lateS
mv GM12878_unsorted_library1_replicate1 SRR16352028_GM12878-unsorted1

####################################
# Binning Reference Genome in 25kb #
####################################
Kronos binning -R $hg38Ref -c $cores -o ../Binned_Genomes/hg38/25kPE -i $hg38BT2 --bin_size 25000 --paired_ends -B $hg38Blist --chr_range 1:22,X,Y

#######
# CNV #
#######
cd Massey2022
# if CNV file exists, skip the dataset. Otherwise, if M(ale), then add both X & Y chromosomes, else: X only
datasetcount="$(ls ./BAM | wc -l)"
tail -n+2 Metadata/PRJNA770772_v2.txt | awk -F'\t' '{ if (system(" test -f Kronos/CNV/"$18"_"$1"_cnv_calls.bed ")) { if ($26 == "male") {chr = "X,Y"} else {chr = "X"} {if ($8 == "Human") {system(" date ; echo "$18"_"$1": dataset "NR" of '$datasetcount' ; Kronos CNV -D BAM/"$1"_"$18" -B ../Binned_Genomes/hg38/25kPE/*.tsv -c '$cores' -o 'Kronos'/CNV -m "$31" -M "$32" -e "$18"_"$1" -n 50000 --chr_range 1:22,"chr)}} } else {print $18,$1,"DONE..Exiting..."} }'

## Cell stats
sh calculate_stats.sh > Metadata/MK_read_stats.csv

# Create cell phase correspondance files
cd Massey2022/Kronos
echo "Cell\tPhase" > ./Metadata/GM12878_ALL_Cell_Phases.tsv
cat CNV/*-G1*.csv | awk -F',' 'OFS="\t" {print $1,"G1"}' | sort | tail -n+6 >> ./Metadata/GM12878_ALL_Cell_Phases.tsv
cat CNV/*-G2*.csv | awk -F',' 'OFS="\t" {print $1,"G1"}' | sort | tail -n+6 >> ./Metadata/GM12878_ALL_Cell_Phases.tsv
cat CNV/*-S_*.csv | awk -F',' 'OFS="\t" {print $1,"S"}' | sort | tail -n+6 >> ./Metadata/GM12878_ALL_Cell_Phases.tsv
cat CNV/*-lateS_*.csv | awk -F',' 'OFS="\t" {print $1,"S"}' | sort | tail -n+6 >> ./Metadata/GM12878_ALL_Cell_Phases.tsv
cat CNV/*-earlyS_*.csv | awk -F',' 'OFS="\t" {print $1,"S"}' | sort | tail -n+6 >> ./Metadata/GM12878_ALL_Cell_Phases.tsv

#Merge CNV of each cell-type:
cat CNV/*GM12878*.csv | sort | uniq > CNV/GM12878_ALL_per_Cell_summary_metrics.csv
cat CNV/*GM12878*.bed | sort | uniq  > CNV/GM12878_ALL_calls.bed
cd Massey2022

##############
# Diagnostic #
##############
# FOR MnM
datasetcount="$(ls ./BAM | wc -l)"
tail -n+2 Metadata/PRJNA770772_v2.txt | head -n5 | awk -F'\t' '{ system(" date ; echo "$18"_"$1": dataset "NR" of '$datasetcount' ; Kronos diagnostic -F Kronos/CNV/"$18"*csv -c 3 -o Kronos/Diagnostic/Dev -b "$18" -m "$33" -G "$34" -S "$35" -d ") }'
# To filter other data for MnM
tail -n+7 Metadata/PRJNA770772_v2.txt | awk -F'\t' '{ system(" date ; echo "$18"_"$1": dataset "NR" of '$datasetcount' ; Kronos diagnostic -F Kronos/CNV/"$18"*csv -c 3 -o Kronos/Diagnostic/Dev -b "$18" -m "$33" -d") }'

#################
##### M n M #####
#################
mkdir Kronos/CNV/merged
# Metadata (cell name, cell type)
for bed in Kronos/CNV/*_per_Cell_summary_metrics.csv ; do
	name=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' | awk -F'_SRR' '{print $1 }')
	tail -n+2 Kronos/CNV/"$name"*_per_Cell_summary_metrics.csv | awk -F',' -v n="$name" 'OFS="\t"{print $1,n}'
done > Metadata/All_MK_CellsAndTypes.tsv

## filter CNV data
tail -n+2 Metadata/PRJNA770772_v2.txt | awk -F'\t' '{ system(" date ; echo "$18" ;  head -n1 Kronos/CNV/"$18"*_cnv_calls.bed > Kronos/CNV/merged/"$30".bed ") }'

for bed in $(ls Kronos/CNV/*_cnv_calls.bed | grep -v G2) ; do
	celltype=$(echo $bed | awk -F'-' '{print $1 }' | sed 's/Kronos\/CNV\///' | sed 's/MCF7/MCF-7/' | sed 's/HCT116/HCT-116/')
	name=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_cnv_calls.bed//' | awk -F'_SRR' '{print $1 }')
	echo $celltype $name
	awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 Kronos/diagnostic/Dev/"$name"_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') Kronos/CNV/"$name"*_cnv_calls.bed >> Kronos/CNV/merged/"$celltype".bed
done

## MnM ##
for bed in $(ls Kronos/CNV/merged/* ) ; do
	celltype=$(echo $bed | sed 's/Kronos\/CNV\/merged\///' | sed 's/.bed//' )
	echo $celltype
	MnM -i $bed -o MnM/$celltype -n $celltype -r -s -b --seed 18671107 --groups Metadata/All_MK_CellsAndTypes.tsv --cpu 4
done

##############
# Diagnostic # For Kronos this time
##############
# WhoIsWho
for bed in Kronos/CNV/*_per_Cell_summary_metrics.csv ; do
	name=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' | awk -F'_SRR' '{print $1 }')
	celltype=$(echo $bed | awk -F'-' '{print $1 }' | sed 's/Kronos\/CNV\///' | sed 's/MCF7/MCF-7/' | sed 's/HCT116/HCT-116/')
	echo $name
	Kronos WhoIsWho -F $bed -o Kronos/Diagnostic/Whoiswho -W MnM/$celltype/*_metadata.tsv
done

# Diagnostic for Kronos
tail -n+2 Metadata/PRJNA770772_v2.txt | awk -F'\t' '{ system(" echo "$18" ; Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_"$18"*.csv -c 3 -o Kronos/Diagnostic -b "$18" -m "$33" -C ") }'
# fix H1 that didn't work
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep2-lane3_SRR16362019_per_Cell_summary_metrics.csv -c 3 -o Kronos/Diagnostic -b H1-unsorted-rep2-lane3 -m 15 -C -f 1 -s 0.5

######
# RT #
######
# GM12878
mkdir -p Kronos/subpopulation_files_split/GM12878
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/GM12878 --csv-files Kronos/Diagnostic/Whoiswho/phased_GM12878-G1_SRR16352033_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12878-S_SRR16352030_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12878-earlyS_SRR16352031_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12878-lateS_SRR16352029_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12878-unsorted1_SRR16352028_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12878-unsorted2_SRR16352027_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12878-unsorted3_SRR16352026_per_Cell_summary_metrics.csv --bed-files MnM/GM12878/GM12878.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_GM12878-G1_SRR16352033_per_Cell_summary_metrics.csv,\
Kronos/Diagnostic/Whoiswho/phased_GM12878-S_SRR16352030_per_Cell_summary_metrics.csv,\
Kronos/Diagnostic/Whoiswho/phased_GM12878-earlyS_SRR16352031_per_Cell_summary_metrics.csv,\
Kronos/Diagnostic/Whoiswho/phased_GM12878-lateS_SRR16352029_per_Cell_summary_metrics.csv,\
Kronos/Diagnostic/Whoiswho/phased_GM12878-unsorted1_SRR16352028_per_Cell_summary_metrics.csv,\
Kronos/Diagnostic/Whoiswho/phased_GM12878-unsorted2_SRR16352027_per_Cell_summary_metrics.csv,\
Kronos/Diagnostic/Whoiswho/phased_GM12878-unsorted3_SRR16352026_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-G1_SRR16352033_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-S_SRR16352030_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-earlyS_SRR16352031_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-lateS_SRR16352029_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-unsorted1_SRR16352028_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-unsorted2_SRR16352027_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/GM12878/GM12878.BED_phased_GM12878-unsorted3_SRR16352026_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/GM12878 \
-b G1,S,earlyS,lateS,unsorted1,unsorted2,unsorted3 \
-f 200kb \
-g GM12878,GM12878,GM12878,GM12878,GM12878,GM12878,GM12878 \
-S Kronos/Diagnostic/GM12878-G1_settings.txt,\
Kronos/Diagnostic/GM12878-S_settings.txt,\
Kronos/Diagnostic/GM12878-earlyS_settings.txt,\
Kronos/Diagnostic/GM12878-lateS_settings.txt,\
Kronos/Diagnostic/GM12878-unsorted1_settings.txt,\
Kronos/Diagnostic/GM12878-unsorted2_settings.txt,\
Kronos/Diagnostic/GM12878-unsorted3_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_1/GSM923451_hg19_wgEncodeUwRepliSeqGm12878WaveSignalRep1_hg38.BedGraph \
--ref_name GM12878\ Bulk\ RT

# GM12891
mkdir Kronos/subpopulation_files_split/GM12891
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/GM12891 --csv-files Kronos/Diagnostic/Whoiswho/phased_GM12891-unsorted-lane1_SRR16328252_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12891-unsorted-lane2_SRR16328251_per_Cell_summary_metrics.csv --bed-files MnM/GM12891/GM12891.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_GM12891-unsorted-lane1_SRR16328252_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_GM12891-unsorted-lane2_SRR16328251_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/GM12891/GM12891.BED_phased_GM12891-unsorted-lane1_SRR16328252_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/GM12891/GM12891.BED_phased_GM12891-unsorted-lane2_SRR16328251_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/GM12891 \
-b lane1,lane2 \
-f 200kb \
-g GM12891,GM12891 \
-S Kronos/Diagnostic/GM12891-unsorted-lane1_settings.txt,Kronos/Diagnostic/GM12891-unsorted-lane2_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_1/GSM923451_hg19_wgEncodeUwRepliSeqGm12878WaveSignalRep1_hg38.BedGraph \
--ref_name GM12878\ Bulk\ RT

# GM12892
mkdir Kronos/subpopulation_files_split/GM12892
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/GM12892 --csv-files Kronos/Diagnostic/Whoiswho/phased_GM12892-unsorted-lane1_SRR16328249_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_GM12892-unsorted-lane2_SRR16328248_per_Cell_summary_metrics.csv --bed-files MnM/GM12892/GM12892.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_GM12892-unsorted-lane1_SRR16328249_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_GM12892-unsorted-lane2_SRR16328248_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/GM12892/GM12892.BED_phased_GM12892-unsorted-lane1_SRR16328249_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/GM12892/GM12892.BED_phased_GM12892-unsorted-lane2_SRR16328248_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/GM12892 \
-b lane1,lane2 \
-f 200kb \
-g GM12892,GM12892 \
-S Kronos/Diagnostic/GM12892-unsorted-lane1_settings.txt,Kronos/Diagnostic/GM12892-unsorted-lane2_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_1/GSM923451_hg19_wgEncodeUwRepliSeqGm12878WaveSignalRep1_hg38.BedGraph \
--ref_name GM12878\ Bulk\ RT

# H1
mkdir Kronos/subpopulation_files_split/H1
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/H1 --csv-files Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep1_SRR16362022_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep2-lane2_SRR16362020_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep2-lane3_SRR16362019_per_Cell_summary_metrics.csv --bed-files MnM/H1/H1.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep1_SRR16362022_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep2-lane2_SRR16362020_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_H1-unsorted-rep2-lane3_SRR16362019_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/H1/H1.BED_phased_H1-unsorted-rep1_SRR16362022_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/H1/H1.BED_phased_H1-unsorted-rep2-lane2_SRR16362020_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/H1/H1.BED_phased_H1-unsorted-rep2-lane3_SRR16362019_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/H1 \
-b rep1,rep2-lane2,rep2-lane3 \
-f 200kb \
-g H1,H1,H1 \
-S Kronos/Diagnostic/H1-unsorted-rep1_settings.txt,Kronos/Diagnostic/H1-unsorted-rep2-lane2_settings.txt,Kronos/Diagnostic/H1-unsorted-rep2-lane3_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_2/wgEncodeFsuRepliChipH1hescWaveSignalRep1_hg38.BedGraph \
--ref_name H1\ Bulk\ RT

# H7
mkdir Kronos/subpopulation_files_split/H7
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/H7 --csv-files Kronos/Diagnostic/Whoiswho/phased_H7-unsorted-lane1_SRR16362018_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_H7-unsorted-lane2_SRR16362017_per_Cell_summary_metrics.csv --bed-files MnM/H7/H7.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_H7-unsorted-lane1_SRR16362018_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_H7-unsorted-lane2_SRR16362017_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/H7/H7.BED_phased_H7-unsorted-lane1_SRR16362018_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/H7/H7.BED_phased_H7-unsorted-lane2_SRR16362017_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/H7 \
-b lane1,lane2 \
-f 200kb \
-g H7,H7 \
-S Kronos/Diagnostic/H7-unsorted-lane1_settings.txt,Kronos/Diagnostic/H7-unsorted-lane2_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_2/wgEncodeFsuRepliChipH7esWaveSignalRep1_hg38.BedGraph \
--ref_name H7\ Bulk\ RT

# H9
mkdir Kronos/subpopulation_files_split/H9
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/H9 --csv-files Kronos/Diagnostic/Whoiswho/phased_H9-unsorted_SRR16362016_per_Cell_summary_metrics.csv --bed-files MnM/H9/H9.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_H9-unsorted_SRR16362016_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/H9/H9.BED_phased_H9-unsorted_SRR16362016_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/H9 \
-b H9 \
-f 200kb \
-g H9 \
-S Kronos/Diagnostic/H9-unsorted_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_2/wgEncodeFsuRepliChipH9esWaveSignalRep1_hg38.BedGraph \
--ref_name H9\ Bulk\ RT

# HCT-116
mkdir Kronos/subpopulation_files_split/HCT-116
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/HCT-116 --csv-files Kronos/Diagnostic/Whoiswho/phased_HCT116-unsorted-lane1_SRR16328247_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_HCT116-unsorted-lane2_SRR16328246_per_Cell_summary_metrics.csv --bed-files MnM/HCT-116/HCT-116_1.BED.gz MnM/HCT-116/HCT-116_2.BED.gz
Kronos RT \
-F Kronos/subpopulation_files_split/HCT-116/HCT-116_1.BED_phased_HCT116-unsorted-lane1_SRR16328247_per_Cell_summary_metrics.csv,Kronos/subpopulation_files_split/HCT-116/HCT-116_1.BED_phased_HCT116-unsorted-lane2_SRR16328246_per_Cell_summary_metrics.csv,Kronos/subpopulation_files_split/HCT-116/HCT-116_2.BED_phased_HCT116-unsorted-lane1_SRR16328247_per_Cell_summary_metrics.csv,Kronos/subpopulation_files_split/HCT-116/HCT-116_2.BED_phased_HCT116-unsorted-lane2_SRR16328246_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/HCT-116/HCT-116_1.BED_phased_HCT116-unsorted-lane1_SRR16328247_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/HCT-116/HCT-116_1.BED_phased_HCT116-unsorted-lane2_SRR16328246_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/HCT-116/HCT-116_2.BED_phased_HCT116-unsorted-lane1_SRR16328247_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/HCT-116/HCT-116_2.BED_phased_HCT116-unsorted-lane2_SRR16328246_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/HCT-116 \
-b lane1\ S1,lane2\ S1,lane1\ S2,lane2\ S2 \
-f 200kb \
-g HCT-116_1,HCT-116_1,HCT-116_2,HCT-116_2 \
-S Kronos/Diagnostic/HCT116-unsorted-lane1_settings.txt,Kronos/Diagnostic/HCT116-unsorted-lane2_settings.txt,Kronos/Diagnostic/HCT116-unsorted-lane1_settings.txt,Kronos/Diagnostic/HCT116-unsorted-lane2_settings.txt \
-B 200000 \
-c $cores \
-N 3 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/Du2021/GSE158008_HCT116_WA_hg38.BedGraph \
--ref_name HCT-116\ Bulk\ RT \

# MCF-7
mkdir Kronos/subpopulation_files_split/MCF-7
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/MCF-7 --csv-files Kronos/Diagnostic/Whoiswho/phased_MCF7-unsorted-lane2_SRR16328242_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_MCF7-unsorted-lane3_SRR16328250_per_Cell_summary_metrics.csv --bed-files MnM/MCF-7/MCF-7.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_MCF7-unsorted-lane2_SRR16328242_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_MCF7-unsorted-lane3_SRR16328250_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/MCF-7/MCF-7.BED_phased_MCF7-unsorted-lane2_SRR16328242_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/MCF-7/MCF-7.BED_phased_MCF7-unsorted-lane3_SRR16328250_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/MCF-7 \
-b lane2,lane3 \
-f 200kb \
-g MCF-7,MCF-7 \
-S Kronos/Diagnostic/MCF7-unsorted-lane2_settings.txt,Kronos/Diagnostic/MCF7-unsorted-lane3_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38/ENCODE_1/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1_hg38.BedGraph \
--ref_name MCF-7\ Bulk\ RT

# RKO
mkdir Kronos/subpopulation_files_split/RKO
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/RKO --csv-files Kronos/Diagnostic/Whoiswho/phased_RKO-unsorted-lane1_SRR16328245_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_RKO-unsorted-lane2_SRR16328244_per_Cell_summary_metrics.csv --bed-files MnM/RKO/RKO_1.BED.gz MnM/RKO/RKO_2.BED.gz
Kronos RT \
-F Kronos/Diagnostic/Whoiswho/phased_RKO-unsorted-lane1_SRR16328245_per_Cell_summary_metrics.csv,Kronos/Diagnostic/Whoiswho/phased_RKO-unsorted-lane2_SRR16328244_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/RKO/RKO_1.BED_phased_RKO-unsorted-lane1_SRR16328245_per_Cell_summary_metrics.bed,Kronos/subpopulation_files_split/RKO/RKO_1.BED_phased_RKO-unsorted-lane2_SRR16328244_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/RKO \
-b lane1,lane2 \
-f 200kb \
-g RKO,RKO \
-S Kronos/Diagnostic/RKO-unsorted-lane1_settings.txt,Kronos/Diagnostic/RKO-unsorted-lane2_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X
