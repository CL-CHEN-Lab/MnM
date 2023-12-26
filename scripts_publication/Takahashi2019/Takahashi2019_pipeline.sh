#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

mkdir -p Takahashi2019/HUMAN/Metadata

#############################
# Download single cell data #
#############################
cd Takahashi2019/HUMAN/fastq
# download
awk -F ";" 'NR>1 { system(" echo "$1" ; fasterq-dump -3 -p --skip-technical -O ./ "$1" ") }' ./Metadata/SRP130912_SraRunTable.csv
# merge paired-end files
awk -F ";" 'NR>1 {if ($19 == "PAIRED") { system("for R in 1 2 ; do cat "$1"_$R.fastq >> "$17"_$R.fastq && rm "$1"_$R.fastq ; done ")} }' ./Metadata/SRP130912_SraRunTable.csv
# rename single-end files to match paired-end files
awk -F ";" 'NR>1 {if ($19 == "SINGLE") { system("cat "$1".fastq >> "$17".fastq && rm "$1".fastq ") } }' ./Metadata/SRP130912_SraRunTable.csv

################
# fastq to bam #
################
hg38BT2=#Path/Bowtie2Index/genome
hg38Ref=/#Path/WholeGenomeFasta/genome.fa
hg38Blist=#Path/hg38-blacklist.v2.bed
Picard=#Path/picard.jar
cores=6
cd Takahashi2019/HUMAN
awk -F";" 'NR>1{print $1,$19}' | sort | uniq | awk '{if ($19 == "SINGLE") {system("Kronos fastqtoBAM --one=fastq/"$1".fastq -i '$hg38Ref' -c '$cores' -o BAM/single --path_to_trim_galore=$(which trim_galore) --path_to_cutadapt=$(which cutadapt) --path_to_java=$(which java) --path_to_picard='$Picard' && pigz fastq/"$1"*.fastq ")} else {system("Kronos fastqtoBAM --one=fastq/"$1"_1.fastq --two=fastq/"$1"_2.fastq -i '$hg38Ref' -c '$cores' -o BAM/paired/ --path_to_trim_galore=$(which trim_galore) --path_to_cutadapt=$(which cutadapt) --path_to_java=$(which java) --path_to_picard='$Picard' && pigz fastq/"$1"*.fastq ")} }' ./Metadata/SRP130912_SraRunTable.csv

###########
# Binning #
###########
## Paired-end
Kronos binning -R $hg38Ref -c $cores -o ../Binned_Genomes/hg38/25kPE -i $hg38BT2 --bin_size 25000 --paired_ends -B $hg38Blist --chr_range 1:22,X,Y
## Single-end
Kronos binning -R $hg38Ref -c $cores -o ../Binned_Genomes/hg38/25kSE -i $hg38BT2 --bin_size 25000 -B $hg38Blist --chr_range 1:22,X,Y

#######
# CNV #
#######
cd Takahashi2019/HUMAN
## Paired-end
Kronos CNV -D BAM/paired/delDup -B ../Binned_Genomes/hg38/25kPE/*.tsv -c $cores -o Kronos/CNV -e P_hTERT-RPE1 -n 75000 --chr_range 1:22,X
## Single-end
Kronos CNV -D BAM/single/delDup -B ../Binned_Genomes/hg38/25kSE/*.tsv -c $cores -o Kronos/CNV -e S_hTERT-RPE1 -n 75000 --chr_range 1:22,X
# Merge CNV files from single and paired
cat Kronos/CNV/*.bed | sort | uniq > Kronos/CNV/hTERT-RPE1_cnv_calls.bed
cat Kronos/CNV/*.csv | sort | uniq > Kronos/CNV/hTERT-RPE1_per_Cell_summary_metrics.csv

############
# WhoIsWho #
############
Kronos WhoIsWho -F Kronos/CNV/hTERT-RPE1_per_Cell_summary_metrics.csv -W Metadata/Phase_Cell_human_correspondance.tsv -o Kronos/WhoIsWho

##############
# Diagnostic #
##############
#identify appropriate thresholds for Kronos RT
Kronos diagnostic -F Kronos/WhoIsWho/phased_hTERT-RPE1_per_Cell_summary_metrics.csv -o Kronos/diagnostic -b hTERT-RPE1 -c 6 -m 70 -f 0.95 -s 0.5 -C

#######
# MnM #
#######
MnM -i Kronos/CNV/hTERT-RPE1_cnv_calls.bed -o MnM -n hTERT-RPE1 -r -b --seed 18671107

######
# RT #
######
#liftover Reference RT from hg19
liftOver Metadata/Reference_RT/GSE108556_RPE1_rpm_w200ks40k_BrdUIP_Percent_q0.05.bedGraph Metadata/Reference_RT/hg19ToHg38.over.chain Metadata/Reference_RT/hTERT-RPE1_hg38_Reference.bed Metadata/Reference_RT/unMapped_DKO1.bed

#scRT
Kronos RT -F Kronos/WhoIsWho/phased_hTERT-RPE1_per_Cell_summary_metrics.csv -T Kronos/CNV/hTERT-RPE1_cnv_calls.bed -C ../genome_assemblies/HUMAN/name_length_hg38.txt -o Kronos/RT -b hTERT-RPE1\ Takahashi2019 -f 200kb -S Kronos/diagnostic/hTERT-RPE1_settings.txt -c 6 -p -r chr2:30Mb-60Mb -B 200000 -N 3 -R Metadata/Reference_RT/hTERT-RPE1_hg38_Reference.bed
