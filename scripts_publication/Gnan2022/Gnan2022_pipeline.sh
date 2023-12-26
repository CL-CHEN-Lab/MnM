#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

#################
# Declare Paths #
#################
OutputDir=Kronos_data_analyses
cores=6

Picard=#path/picard.jar
hg38BT2=#path/Bowtie2Index/genome
hg38Ref=#path/WholeGenomeFasta/genome.fa
hg38Blist=#path/hg38-blacklist.v2.bed

BAMDir=#path
RefRTs=#path
Binning_dir=../Binned_Genomes/hg38/20kPE
cd $OutputDir

###########
# Binning #
###########
Kronos binning -R $hg38Ref -c $cores -o $Binning_dir -i $hg38BT2 --bin_size 20000 --paired_ends -B $hg38Blist --chr_range 1:22,X,Y

#######
# CNV #
#######
date ; echo JEFF_1 ; Kronos CNV -D $BAMDir/Hela1 -B ../Binned_Genomes/hg38/20kPE/*.tsv -c $cores -o $OutputDir/CNV/20kb_2023 -e JEFF_1 --chr_range 1:22,X,Y -m 1.25 -M 4
date ; echo JEFF_2 ; Kronos CNV -D $BAMDir/Hela2 -B ../Binned_Genomes/hg38/20kPE/*.tsv -c $cores -o $OutputDir/CNV/20kb_2023 -e JEFF_2 --chr_range 1:22,X,Y -m 1.25 -M 4
date ; echo HeLa_1 ; Kronos CNV -D $BAMDir/Jeff1 -B ../Binned_Genomes/hg38/20kPE/*.tsv -c $cores -o $OutputDir/CNV/20kb_2023 -e HeLa_1 --chr_range 1:22,X -m 2 -M 6
date ; echo HeLa_2 ; Kronos CNV -D $BAMDir/Jeff2 -B ../Binned_Genomes/hg38/20kPE/*.tsv -c $cores -o $OutputDir/CNV/20kb_2023 -e HeLa_2 --chr_range 1:22,X -m 2 -M 6
date ; echo MCF7_unsorted ; Kronos CNV -D $BAMDir/First_exp -B ../Binned_Genomes/hg38/20kPE/*.tsv -c $cores -o $OutputDir/CNV/20kb_2023 -e MCF7_unsorted --chr_range 1:22,X -m 2 -M 6.5
date ; echo MCF7_Normal ; Kronos CNV -D $BAMDir/Normal -B ../Binned_Genomes/hg38/20kPE/*.tsv -c $cores -o $OutputDir/CNV/20kb_2023 -e MCF7_Normal --chr_range 1:22,X -m 2 -M 6.5

## Cell stats
sh calculate_stats.sh CNV/20kb_2023 > Metadata/Kronos_read_stats.csv

##############
# Diagnostic # For replication state predictor (dev mode) -d
##############
Kronos diagnostic -F CNV/20kb_2023/HeLa_1_per_Cell_summary_metrics.csv -o diagnostic/Dev -b HeLa_1 -m 115 -c 3 -G 0.8 -S 0.85 -d
Kronos diagnostic -F CNV/20kb_2023/HeLa_2_per_Cell_summary_metrics.csv -o diagnostic/Dev -b HeLa_2 -m 115 -c 3 -G 0.65 -S 0.7 -d
Kronos diagnostic -F CNV/20kb_2023/JEFF_1_per_Cell_summary_metrics.csv -o diagnostic/Dev -b JEFF_1 -m 115 -c 3 -G 0.7 -S 0.75 -d
Kronos diagnostic -F CNV/20kb_2023/JEFF_2_per_Cell_summary_metrics.csv -o diagnostic/Dev -b JEFF_2 -m 115 -c 3 -G 0.6 -S 0.65 -d
Kronos diagnostic -F CNV/20kb_2023/MCF7_unsorted_per_Cell_summary_metrics.csv -o diagnostic/Dev -b MCF7_unsorted -m 115 -c 3 -G 1.025 -S 1.1 -d
Kronos diagnostic -F CNV/20kb_2023/MCF7_Normal_per_Cell_summary_metrics.csv -o diagnostic/Dev -b MCF7_Normal -m 115 -c 3 -G .7 -S .75 -d

#################
##### M n M #####
#################
cd Kronos_data_analyses/
mkdir CNV/20kb_2023/merged
mkdir Metadata

#Group R1 and R2 toghether
## HeLa
head -n1 CNV/20kb_2023/HeLa_1_cnv_calls.bed > CNV/20kb_2023/merged/HeLa.bed
awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 diagnostic/Dev/HeLa_1_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') CNV/20kb_2023/HeLa_1_cnv_calls.bed >> CNV/20kb_2023/merged/HeLa.bed
awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 diagnostic/Dev/HeLa_2_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') CNV/20kb_2023/HeLa_2_cnv_calls.bed >> CNV/20kb_2023/merged/HeLa.bed
## JEFF
head -n1 CNV/20kb_2023/JEFF_1_cnv_calls.bed > CNV/20kb_2023/merged/JEFF.bed
awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 diagnostic/Dev/JEFF_1_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') CNV/20kb_2023/JEFF_1_cnv_calls.bed >> CNV/20kb_2023/merged/JEFF.bed
awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 diagnostic/Dev/JEFF_2_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') CNV/20kb_2023/JEFF_2_cnv_calls.bed >> CNV/20kb_2023/merged/JEFF.bed
## MCF7
head -n1 CNV/20kb_2023/MCF7_Normal_cnv_calls.bed > CNV/20kb_2023/merged/MCF7.bed
awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 diagnostic/Dev/MCF7_Normal_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') CNV/20kb_2023/MCF7_Normal_cnv_calls.bed >> CNV/20kb_2023/merged/MCF7.bed
awk 'FNR==NR { patterns[$0]; next } ($1 in patterns)' <(tail -n+2 diagnostic/Dev/MCF7_unsorted_phases_from_diagnostic.tsv | grep -v 'Low Coverage' | awk -F'\t' '{print $1}') CNV/20kb_2023/MCF7_unsorted_cnv_calls.bed >> CNV/20kb_2023/merged/MCF7.bed
## HeLa + JEFF
(cat CNV/20kb_2023/merged/JEFF.bed ; tail -n+2 CNV/20kb_2023/merged/HeLa.bed ) > CNV/20kb_2023/merged/Mix_HeLa_JEFF.bed

## cell phases file
tail -n+2 CNV/20kb_2023/HeLa_1_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"HeLa_R1"}' > Metadata/AllKronosCellsAndTypes.tsv
tail -n+2 CNV/20kb_2023/HeLa_2_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"HeLa_R2"}' >> Metadata/AllKronosCellsAndTypes.tsv
tail -n+2 CNV/20kb_2023/JEFF_1_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"JEFF_R1"}' >> Metadata/AllKronosCellsAndTypes.tsv
tail -n+2 CNV/20kb_2023/JEFF_2_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"JEFF_R2"}' >> Metadata/AllKronosCellsAndTypes.tsv
tail -n+2 CNV/20kb_2023/MCF7_Normal_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"MCF7"}' >> Metadata/AllKronosCellsAndTypes.tsv
tail -n+2 CNV/20kb_2023/MCF7_unsorted_per_Cell_summary_metrics.csv | awk -F',' 'OFS="\t" {print $1,"MCF7_unsorted"}' >> Metadata/AllKronosCellsAndTypes.tsv

#### MnM
echo ; echo HeLa
MnM -i CNV/20kb_2023/merged/HeLa.bed -o MnM/HeLa -n HeLa -r -s --groups Metadata/AllKronosCellsAndTypes.tsv -b --seed 18671107
echo ; echo JEFF
MnM -i CNV/20kb_2023/merged/JEFF.bed -o MnM/JEFF -n JEFF -r -s --groups Metadata/AllKronosCellsAndTypes.tsv -b --seed 18671107
echo ; echo MCF7
MnM -i CNV/20kb_2023/merged/MCF7.bed -o MnM/MCF7 -n MCF-7 -r -s --groups Metadata/AllKronosCellsAndTypes.tsv -b --seed 18671107
echo ; echo HeLa JEFF mix
MnM -i <(cat CNV/20kb_2023/merged/Mix_HeLa_JEFF.bed | grep -v chrY) -o MnM/Mix_HeLa_JEFF -n Mix\ HeLa+JEFF -r -s --groups Metadata/AllKronosCellsAndTypes.tsv --seed 18671107

##############
# Diagnostic # For Kronos this time
##############

# WhoIsWho module from MnM Replication State Predictor
Kronos WhoIsWho -F CNV/20kb_2023/HeLa_1_per_Cell_summary_metrics.csv -o diagnostic/Whoiswho -W MnM/HeLa/HeLa_metadata.tsv
Kronos WhoIsWho -F CNV/20kb_2023/HeLa_2_per_Cell_summary_metrics.csv -o diagnostic/Whoiswho -W MnM/HeLa/HeLa_metadata.tsv
Kronos WhoIsWho -F CNV/20kb_2023/JEFF_1_per_Cell_summary_metrics.csv -o diagnostic/Whoiswho -W MnM/JEFF/JEFF_metadata.tsv
Kronos WhoIsWho -F CNV/20kb_2023/JEFF_2_per_Cell_summary_metrics.csv -o diagnostic/Whoiswho -W MnM/JEFF/JEFF_metadata.tsv
Kronos WhoIsWho -F CNV/20kb_2023/MCF7_unsorted_per_Cell_summary_metrics.csv -o diagnostic/Whoiswho -W MnM/MCF7/MCF-7_metadata.tsv
Kronos WhoIsWho -F CNV/20kb_2023/MCF7_Normal_per_Cell_summary_metrics.csv -o diagnostic/Whoiswho -W MnM/MCF7/MCF-7_metadata.tsv

# Diagnostic for Kronos
Kronos diagnostic -F diagnostic/Whoiswho/phased_HeLa_1_per_Cell_summary_metrics.csv -o diagnostic -b HeLa_1 -m 115 -c 3 -C
Kronos diagnostic -F diagnostic/Whoiswho/phased_HeLa_2_per_Cell_summary_metrics.csv -o diagnostic -b HeLa_2 -m 115 -c 3 -C
Kronos diagnostic -F diagnostic/Whoiswho/phased_JEFF_1_per_Cell_summary_metrics.csv -o diagnostic -b JEFF_1 -m 115 -c 3 -C
Kronos diagnostic -F diagnostic/Whoiswho/phased_JEFF_2_per_Cell_summary_metrics.csv -o diagnostic -b JEFF_2 -m 115 -c 3 -C
Kronos diagnostic -F diagnostic/Whoiswho/phased_MCF7_unsorted_per_Cell_summary_metrics.csv -o diagnostic -b MCF7_unsorted -m 115 -c 3 -C
Kronos diagnostic -F diagnostic/Whoiswho/phased_MCF7_Normal_per_Cell_summary_metrics.csv -o diagnostic -b MCF7_Normal -m 115 -c 3 -C

######################
# Replication Timing # RT from Kronos
######################
#Create chromosome size file
#retain 2 first columns (chr number and length)
# (echo "chr\tsize" ; cut -f1,2 #path/WholeGenomeFasta/genome.fa.fai) > Metadata/Chr_size.txt

### RT ###

# HeLa
mkdir -p subpopulation_files_split/HeLa
python split_subpopulations_csv_bed.py --output-dir subpopulation_files_split/HeLa --csv-files diagnostic/Whoiswho/phased_HeLa_1_per_Cell_summary_metrics.csv diagnostic/Whoiswho/phased_HeLa_2_per_Cell_summary_metrics.csv --bed-files MnM/HeLa/HeLa.BED.gz
Kronos RT \
-F diagnostic/Whoiswho/phased_HeLa_1_per_Cell_summary_metrics.csv,diagnostic/Whoiswho/phased_HeLa_2_per_Cell_summary_metrics.csv \
-T subpopulation_files_split/HeLa/HeLa.BED_phased_HeLa_1_per_Cell_summary_metrics.bed,subpopulation_files_split/HeLa/HeLa.BED_phased_HeLa_2_per_Cell_summary_metrics.bed \
-C Metadata/Chr_size.txt \
-r Metadata/regions_to_plot.bed \
-o RT/HeLa \
-b R1,R2 \
-f 200kb \
-g HeLa,HeLa \
-S diagnostic/HeLa_1_settings.txt,diagnostic/HeLa_2_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38//HeLaS3_hg38.bg \
--ref_name HeLa\ Bulk\ RT

# JEFF
mkdir subpopulation_files_split/JEFF
python split_subpopulations_csv_bed.py --output-dir subpopulation_files_split/JEFF --csv-files diagnostic/Whoiswho/phased_JEFF_1_per_Cell_summary_metrics.csv diagnostic/Whoiswho/phased_JEFF_2_per_Cell_summary_metrics.csv --bed-files MnM/JEFF/JEFF.BED.gz
Kronos RT \
-F diagnostic/Whoiswho/phased_JEFF_1_per_Cell_summary_metrics.csv,diagnostic/Whoiswho/phased_JEFF_2_per_Cell_summary_metrics.csv \
-T subpopulation_files_split/JEFF/JEFF.BED_phased_JEFF_1_per_Cell_summary_metrics.bed,subpopulation_files_split/JEFF/JEFF.BED_phased_JEFF_2_per_Cell_summary_metrics.bed \
-C Metadata/Chr_size.txt \
-r Metadata/regions_to_plot.bed \
-o RT/JEFF \
-b R1,R2 \
-f 200kb \
-g JEFF,JEFF \
-S diagnostic/JEFF_1_settings.txt,diagnostic/JEFF_2_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R /RT_references/hg38//GM12878_hg38.bg \
--ref_name GM12878\ Bulk\ RT

# MCF-7
mkdir subpopulation_files_split/MCF-7
python split_subpopulations_csv_bed.py --output-dir subpopulation_files_split/MCF-7 --csv-files diagnostic/Whoiswho/phased_MCF7_Normal_per_Cell_summary_metrics.csv diagnostic/Whoiswho/phased_MCF7_unsorted_per_Cell_summary_metrics.csv --bed-files MnM/MCF7/MCF-7_1.BED.gz MnM/MCF7/MCF-7_2.BED.gz
Kronos RT \
-F subpopulation_files_split/MCF-7/MCF-7_1.BED_phased_MCF7_Normal_per_Cell_summary_metrics.csv,subpopulation_files_split/MCF-7/MCF-7_2.BED_phased_MCF7_Normal_per_Cell_summary_metrics.csv,subpopulation_files_split/MCF-7/MCF-7_1.BED_phased_MCF7_unsorted_per_Cell_summary_metrics.csv,subpopulation_files_split/MCF-7/MCF-7_2.BED_phased_MCF7_unsorted_per_Cell_summary_metrics.csv \
-T subpopulation_files_split/MCF-7/MCF-7_1.BED_phased_MCF7_Normal_per_Cell_summary_metrics.bed,subpopulation_files_split/MCF-7/MCF-7_2.BED_phased_MCF7_Normal_per_Cell_summary_metrics.bed,subpopulation_files_split/MCF-7/MCF-7_1.BED_phased_MCF7_unsorted_per_Cell_summary_metrics.bed,subpopulation_files_split/MCF-7/MCF-7_2.BED_phased_MCF7_unsorted_per_Cell_summary_metrics.bed \
-C ../Kronos_data_analyses/Metadata/Chr_size.txt \
-r Metadata/regions_to_plot.bed \
-o RT/MCF-7_Normal_Unsorted \
-b Normal\ 1,Normal\ 2,Unsorted\ 1,Unsorted\ 2 \
-f 200kb \
-g MCF-7\ \S1,MCF-7\ S2,MCF-7\ S1,MCF-7\ S2 \
-S diagnostic/MCF7_Normal_settings.txt,diagnostic/MCF7_unsorted_settings.txt,diagnostic/MCF7_Normal_settings.txt,diagnostic/MCF7_unsorted_settings.txt \
-B 200000 \
-c $cores \
-N 5 \
-p \
--extract_G1_G2_cells \
--chr_range 1:22,X \
-R  /RT_references/hg38//MCF7_hg38.bg \
--ref_name MCF-7\ Bulk\ RT


###### DRED ########
Kronos DRed -C RT/HeLa/200kb_single_cells_CNV_200Kb.tsv -o DRed -f HeLa -U --CNV_values B -c 6 -s 18671107
Kronos DRed -C RT/JEFF -o DRed -f JEFF -U --CNV_values B -c 6 -s 18671107
Kronos DRed -C RT/MCF-7_Normal_Unsorted/200kb_single_cells_CNV_200Kb.tsv -o DRed -f MCF-7 -U --CNV_values B -c 6 -s 18671107
Kronos DRed -C RT/MCF-7_Normal_Unsorted/200kb_single_cells_CNV_200Kb.tsv -o DRed -f MCF-7 -U --CNV_values B -c 6 -s 18671107 --per_Chr
