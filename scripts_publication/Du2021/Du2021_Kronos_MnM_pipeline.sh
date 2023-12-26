#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

#########
# PATHS #
#########
OutputDir=Du2021/Kronos
cores=8

TrimGalore=$(which trim_galore)
Cutadapt=$(which cutadapt)
Java=$(which java)
Picard=#path/picard.jar

hg38BT2=#path/Bowtie2Index/genome
hg38Ref=#path/WholeGenomeFasta/genome.fa
hg38Blist=#path/hg38-blacklist.v2.bed

BAMDir=Du2021/BAM
RefRTs=Du2021/Bulk_RepliSeq

cd $OutputDir

####################################
# Binning Reference Genome in 20kb #
####################################
echo "Binning the Reference Genome"
Kronos binning -R $hg38Ref -c $cores -o $OutputDir/Binning -i $hg38BT2 --bin_size 20000 --paired_ends -B $hg38Blist

#######
# CNV #
#######
cd Du2021/Kronos
echo "Estimating Copy-Number Variation"
for BAMs in BAM/* ; do
  echo $(basename "${BAMs}")
  Kronos CNV -D $BAMs -B ../Binned_Genomes/hg38/25kPE/*.tsv -c $cores -o Kronos/CNV -e $(basename "${BAMs}") -m 1.25 -M 4 -n 150000 --chr_range 1:22,X,Y
done

## Cell stats
cd Du2021
sh calculate_stats.sh > Metadata/Du_read_stats.csv

# Create cell phase correspondance files
cd Du2021/Kronos
echo "Cell\tPhase" > ./Metadata/HCT116_ALL_Cell_Phases.tsv
cat CNV/*G1*.csv | awk -F',' 'OFS="\t" {print $1,"G1"}' | sort | tail -n+3 >> ./Metadata/HCT116_ALL_Cell_Phases.tsv
cat CNV/*S*.csv | awk -F',' 'OFS="\t" {print $1,"S"}' | sort | tail -n+3 >> ./Metadata/HCT116_ALL_Cell_Phases.tsv

#Merge CNV of each cell-type:
# merge G and S phases for WT and DKO1 and then all together.
cd Du2021
mkdir Kronos/CNV/merged
# BED files
cat Kronos/CNV/HCT116_WT_scRepli-Seq_G1_cnv_calls.bed > Kronos/CNV/merged/HCT116_WT_cnv_calls.bed ; tail -n+2 Kronos/CNV/HCT116_WT_scRepli-Seq_S_cnv_calls.bed >> Kronos/CNV/merged/HCT116_WT_cnv_calls.bed
cat Kronos/CNV/HCT116_DKO1_scRepli-Seq_G1_cnv_calls.bed > Kronos/CNV/merged/HCT116_DKO1_cnv_calls.bed ; tail -n+2 Kronos/CNV/HCT116_DKO1_scRepli-Seq_S_cnv_calls.bed >> Kronos/CNV/merged/HCT116_DKO1_cnv_calls.bed
cat Kronos/merged/CNV/HCT116_WT_cnv_calls.bed > Kronos/CNV/merged/HCT116_ALL_cnv_calls.bed ; tail -n+2 Kronos/CNV/merged/HCT116_DKO1_cnv_calls.bed >> Kronos/CNV/merged/HCT116_ALL_cnv_calls.bed
# CSV files
cat Kronos/CNV/HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv > Kronos/CNV/merged/HCT116_WT_per_Cell_summary_metrics.csv ; tail -n+2 Kronos/CNV/HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv >> Kronos/CNV/merged/HCT116_WT_per_Cell_summary_metrics.csv
cat Kronos/CNV/HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv > Kronos/CNV/merged/HCT116_DKO1_per_Cell_summary_metrics.csv ; tail -n+2 Kronos/CNV/HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv >> Kronos/CNV/merged/HCT116_DKO1_per_Cell_summary_metrics.csv
cat Kronos/CNV/HCT116_WT_per_Cell_summary_metrics.csv > Kronos/CNV/merged/HCT116_ALL_per_Cell_summary_metrics.csv ; tail -n+2 Kronos/CNV/HCT116_DKO1_per_Cell_summary_metrics.csv >> Kronos/CNV/merged/HCT116_ALL_per_Cell_summary_metrics.csv

##############
# Diagnostic #
##############
## for MnM
Kronos diagnostic -F Kronos/CNV/HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Dev -b HCT116_G1_WT -c 3 -m 400 -S 1.15 -d
Kronos diagnostic -F Kronos/CNV/HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Dev -b HCT116_S_WT -c 3 -m 400 -S 0.5 -d
Kronos diagnostic -F Kronos/CNV/HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Dev -b HCT116_G1_DKO1 -c 3 -m 400 -S 1.05 -d
Kronos diagnostic -F Kronos/CNV/HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Dev -b HCT116_S_DKO1 -c 3 -m 400 -G .68 -S 0.75 -d

##### MNM #######
cd Du2021

# Organise metadata
tail -n+2 Kronos/CNV/HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"HCT116_WT_G1"}' > Metadata/DuCellsAndTypes.tsv
tail -n+2 Kronos/CNV/HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"HCT116_WT_S"}' >> Metadata/DuCellsAndTypes.tsv
tail -n+2 Kronos/CNV/HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"HCT116_DKO1_G1"}' >> Metadata/DuCellsAndTypes.tsv
tail -n+2 Kronos/CNV/HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"HCT116_DKO1_S"}' >> Metadata/DuCellsAndTypes.tsv

# MnM
MnM -i Kronos/CNV/merged/HCT116_WT_cnv_calls.bed -o MnM/WT -n HCT116\ WT --seed 18671107 -r -s -b --groups Metadata/DuCellsAndTypes.tsv
MnM -i Kronos/CNV/merged/HCT116_DKO1_cnv_calls.bed -o MnM/DKO1 -n HCT116\ DKO1 --seed 18671107 -r -s -b --groups Metadata/DuCellsAndTypes.tsv
MnM -i Kronos/CNV/merged/HCT116_ALL_cnv_calls.bed -o MnM/ALL -n HCT116-Du2021 --seed 18671107 -r -s -b --groups Metadata/DuCellsAndTypes.tsv

##############
# Diagnostic # For Kronos this time
##############
# WhoIsWho
Kronos WhoIsWho -F Kronos/CNV/HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Whoiswho -W MnM/WT/HCT116\ WT_metadata.tsv
Kronos WhoIsWho -F Kronos/CNV/HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Whoiswho -W MnM/WT/HCT116\ WT_metadata.tsv
Kronos WhoIsWho -F Kronos/CNV/HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Whoiswho -W MnM/DKO1/HCT116\ DKO1_metadata.tsv
Kronos WhoIsWho -F Kronos/CNV/HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv -o Kronos/Diagnostic/Whoiswho -W MnM/DKO1/HCT116\ DKO1_metadata.tsv

# Diagnostic for Kronos
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b HCT-116_WT_G1 -m 400 -C
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b HCT-116_WT_S -m 400 -C
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b HCT-116_DKO1_G1 -m 400 -C -f 0.95 -s 0.55
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b HCT-116_DKO1_S -m 400 -C

######
# RT #
######
mkdir -p Kronos/subpopulation_files_split/DKO1 Kronos/subpopulation_files_split/WT
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/WT --csv-files Kronos/Diagnostic/Whoiswho/phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv --bed-files MnM/WT/HCT116\ WT_1.BED.gz MnM/WT/HCT116\ WT_2.BED.gz
python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/DKO1 --csv-files Kronos/Diagnostic/Whoiswho/phased_HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv Kronos/Diagnostic/Whoiswho/phased_HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv --bed-files MnM/DKO1/HCT116\ DKO1.BED.gz

# DKO1
Kronos RT \
-F Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/DKO1 \
-b G1,S \
-f 200kb \
-g HCT-116\ DKO1,HCT-116\ DKO1 \
-S Kronos/Diagnostic/HCT-116_DKO1_G1_settings.txt,Kronos/Diagnostic/HCT-116_DKO1_S_settings.txt \
-B 200000 \
-c 8 \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R Bulk_RepliSeq/GSE158008_DKO1_WA_hg38.bed \
--ref_name HCT-116\ DKO1\ Bulk\ RT

# WT
Kronos RT \
-F Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/WT \
-b G1\ 1,S\ 1,\G1\ 2,S\ 2 \
-f 200kb \
-g HCT-116\ WT\ S1,HCT-116\ WT\ S1,HCT-116\ WT\ S2,HCT-116\ WT\ S2 \
-S Kronos/Diagnostic/HCT-116_WT_G1_settings.txt,Kronos/Diagnostic/HCT-116_WT_S_settings.txt,Kronos/Diagnostic/HCT-116_WT_G1_settings.txt,Kronos/Diagnostic/HCT-116_WT_S_settings.txt \
-B 200000 \
-c 4 \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R Bulk_RepliSeq/GSE158008_HCT116_WA_hg38.bed \
--ref_name HCT-116\ Bulk\ RT

# Both
Kronos RT \
-F Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.csv,\
Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.csv \
-T Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_1.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_G1_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/WT/HCT116\ WT_2.BED_phased_HCT116_WT_scRepli-Seq_S_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_G1_per_Cell_summary_metrics.bed,\
Kronos/subpopulation_files_split/DKO1/HCT116\ DKO1.BED_phased_HCT116_DKO1_scRepli-Seq_S_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r ../Kronos_data_analyses/Metadata/regions_to_plot.bed \
-o Kronos/RT/WT_and_DKO1 \
-b G1\ 1,S\ 1,\G1\ 2,S\ 2,G1,S \
-f 200kb \
-g HCT-116\ WT\ S1,HCT-116\ WT\ S1,HCT-116\ WT\ S2,HCT-116\ WT\ S2,HCT-116\ DKO1,HCT-116\ DKO1 \
-S Kronos/Diagnostic/HCT-116_WT_G1_settings.txt,Kronos/Diagnostic/HCT-116_WT_S_settings.txt,Kronos/Diagnostic/HCT-116_WT_G1_settings.txt,Kronos/Diagnostic/HCT-116_WT_S_settings.txt,Kronos/Diagnostic/HCT-116_DKO1_G1_settings.txt,Kronos/Diagnostic/HCT-116_DKO1_S_settings.txt \
-B 200000 \
-c 6 \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R Bulk_RepliSeq/GSE158008_HCT116_WA_hg38.bed \
--ref_name HCT-116\ Bulk\ RT



###### DRED ########
Kronos DRed -C Kronos/RT/WT/200kb_single_cells_CNV_200Kb.tsv -o Kronos/DRed -f HCT-116 -U --CNV_values B -c 6 -s 18671107
