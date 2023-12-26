#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

# DATA: https://www-ncbi-nlm-nih-gov.insb.bib.cnrs.fr/Traces/study/?acc=PRJNA629885&f=assay_type_s%3An%3Awgs%3Blibrarysource_s%3An%3Agenomic%3Bcell_line_sam_ss_dpl110_ss%3An%3Acell%2520line%2520157%2Ccell%2520line%2520bt20%2Ccell%2520line%2520453%2Ccellline231%2Cnot%2520collected%3Ac%3Blibrarylayout_s%3An%3Blibraryselection_s%3An%3Bbytes_l%3An%3Bbiomaterial_provider_sam_ss%3An%3Bbases_l%3An%3Bavgspotlen_l%3An&o=biosample_s%3Aa

mkdir -p Minussi2021/Metadata
mkdir -p Minussi2021/SRA
mkdir -p Minussi2021/fastq
mkdir -p Minussi2021/BAM_hg38
cd Minussi2021

################################
#### DOWNLOAD DATA FROM SRA ####
################################
awk -F';' 'NR>1{  system(" echo "$1" ; fasterq-dump "$1" -O fastq/"$28" -m 750MB -e 4 -p -3  ") }' Metadata/Minussi2021_data_summary.csv

# move 10X data to another folder
mkdir -p fastq/TN3_10x fastq/TN1_10x
mv fastq/TN1/SRR14132810* fastq/TN1_10x/
mv fastq/TN3/SRR14132809* fastq/TN3_10x/

# Demultiplex 10X files.
mkdir Minussi2021/fastq/TN1_10x/SCs
mkdir Minussi2021/fastq/TN3_10x/SCs
demultiplex demux -r -e 16 -p fastq/TN1_10x/SCs/ MasseyKoren2022/Metadata/10x_barcode_whitelist.tsv fastq/TN1_10x/TN1_S2_L009_corrected_R1.fastq fastq/TN1_10x/TN1_S2_L009_corrected_R2.fastq && pigz fastq/TN1_10x/TN1_S2_L009_corrected_R1.fastq fastq/TN1_10x/TN1_S2_L009_corrected_R2.fastq
demultiplex demux -r -e 16 -p fastq/TN3_10x/SCs/ MasseyKoren2022/Metadata/10x_barcode_whitelist.tsv fastq/TN3_10x/SRR14132809_1.fastq fastq/TN3_10x/SRR14132809_2.fastq && pigz fastq/TN3_10x/SRR14132809_1.fastq fastq/TN3_10x/SRR14132809_2.fastq

# Record Number of lines for 10X datasets
cd Minussi2021
for celltype in fastq/TN*_10x ; do
	for fastq in $celltype/SCs/*.fastq ; do wc -l $fastq ; done
done > Metadata/wc-l_SRR_fastq.txt
# Convert file to make easier to read in R
cat Metadata/wc-l_SRR_fastq.txt | sed 's/  */\t/g' | awk ' OFS="\t" {print $1,$2}' | sed 's/fastq\///' > Metadata/wc-l_SRR_fastq_tab.txt

# Run following script to create files to faciliate cell filtering and mapping
Rscript ../SCutoff_v1.R -f Metadata/wc-l_SRR_fastq_tab.txt -n TN1_10x -o Metadata

# Move filtered cells to new directory
mkdir -p fastq/TN1_10x/Retained_SCs fastq/TN3_10x/Retained_SCs
cd Minussi2021/fastq
awk 'NR>1 { system(" printf "$5"\,\   ; mv "$2" "$6" ") }' ./Metadata/2023.05.10_CellsToRetain_TN_10x.tsv

#remove unused cells
cd Minussi2021/fastq
rm -r TN*10x/SCs

#Trim adaptors from reads
cd Minussi2021
# 10x files (paired end)
filecount="$(cat ./Metadata/2023.05.10_PairedMatches_TN_10x.tsv | wc -l)"
awk ' { system(" echo ; date ; echo "$7": file "NR" of '$filecount' ; trim_galore --paired "$1" "$2" --output_dir "$3" --clip_R1 16 --fastqc --cores 4 --gzip && echo Compressing "$7" ; pigz "$1" "$2" ") }' ./Metadata/2023.05.10_PairedMatches_TN_10x.tsv
# SRA files (single or paired end accordingly)
filecount="$(cat Metadata/Minussi2021_data_summary.csv | wc -l)"
awk -F';' 'NR>1{ if ($4 != "202") { if ($20 == "SINGLE") {system(" echo "$1": single file "NR" of '$filecount'  ; trim_galore fastq/"$28"/"$1".fastq --output_dir fastq/"$28"/Trimmed_SCs --fastqc --cores 5 --gzip && echo Compressing "$1" && pigz fastq/"$28"/"$1".fastq ")} else {system("echo "$1": paired files "NR" of '$filecount'; trim_galore --paired fastq/"$28"/"$1"_1.fastq fastq/"$28"/"$1"_2.fastq --output_dir fastq/"$28"/Trimmed_SCs --fastqc --cores 5 --gzip && echo Compressing "$1" && pigz fastq/"$28"/"$1"_1.fastq fastq/"$28"/"$1"_2.fastq ")} } }' Metadata/Minussi2021_data_summary.csv

# Quality control & create required directories for mapping + stats
cd Minussi2021/fastq
for dir in * ; do mkdir -p ./Metadata/QC/$dir/Picard_metrics ./Metadata/QC/$dir/Alignment_stats ../BAM/$dir && multiqc -o ./Metadata/QC/$dir $dir/Trimmed_SCs/ ; done

#declare paths
hg38=#Path/genome.fa
picard=#Path/picard.jar
cd Minussi2021
# 10X files
filecount="$(cat Metadata/2023.05.10_PairedMatches_TN_10x.tsv | wc -l)"
awk ' { system(" date ; \
echo Decompress "$7" : file "NR" of '$filecount' ; unpigz fastq/"$4".gz fastq/"$5".gz && \
echo Align "$7" with BWA && bwa mem -M -t 6 '$hg38' fastq/"$4" fastq/"$5" > ./BAM/"$6"/"$7"-1.sam && \
echo Convert "$7" from SAM to BAM && samtools fixmate -@ 3 -O bam ./BAM/"$6"/"$7"-1.sam ./BAM/"$6"/"$7"-1_unsorted.bam && \
echo Sort "$7" by coordinates && samtools sort -@ 3 -O bam -o ./BAM/"$6"/"$7"-1_sorted.bam ./BAM/"$6"/"$7"-1_unsorted.bam && \
echo Remove duplicates of "$7" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAM/"$6"/"$7"-1_sorted.bam OUTPUT=BAM/"$6"/"$7"-1.bam METRICS_FILE=Metadata/QC/"$6"/Picard_metrics/"$7"_metrics.out && \
echo Index "$7" && samtools index -@ 3 -b BAM/"$6"/"$7"-1.bam && \
echo Generate alignment statistics "$7" && samtools flagstat -@ 3 BAM/"$6"/"$7"-1.bam > Metadata/QC/"$6"/Alignment_stats/"$7"_stats.out && \
rm fastq/"$4" fastq/"$5" && rm ./BAM/"$6"/"$7"-1.sam && rm ./BAM/"$6"/"$7"-1_unsorted.bam && rm ./BAM/"$6"/"$7"-1_sorted.bam ") }' Metadata/2023.05.10_PairedMatches_TN_10x.tsv
# SRA files
filecount="$(cat Metadata/Minussi2021_data_summary.csv | wc -l)"
awk -F';' 'NR>1{ if ($4 != "202") { if ($20 == "SINGLE") {system(" date ; \
echo "$1" : file "NR" of '$filecount' ; unpigz fastq/"$28"/Trimmed_SCs/"$1"*.gz && \
echo Align "$1" with BWA && bwa mem -M -t 6 '$hg38' fastq/"$28"/Trimmed_SCs/"$1"_trimmed.fq > ./BAM/"$28"/"$1".sam && \
echo Sort "$1" by coordinates && samtools sort -@ 4 -O bam -o ./BAM/"$28"/"$1"_sorted.bam ./BAM/"$28"/"$1".sam && \
echo Remove duplicates of "$1" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAM/"$28"/"$1"_sorted.bam OUTPUT=BAM/"$28"/"$1".bam METRICS_FILE=Metadata/QC/"$28"/Picard_metrics/"$1"_metrics.out && \
echo Index "$1" && samtools index -@ 4 -b BAM/"$28"/"$1".bam && \
echo Generate alignment statistics "$1" && samtools flagstat -@ 4 BAM/"$28"/"$1".bam > Metadata/QC/"$28"/Alignment_stats/"$1"_stats.out && \
rm fastq/"$28"/Trimmed_SCs/"$1"*.fq && rm ./BAM/"$28"/"$1".sam && rm ./BAM/"$28"/"$1"_sorted.bam ")} else {system(" date ; \
echo "$1" : file "NR" of '$filecount' ; unpigz fastq/"$28"/Trimmed_SCs/"$1"*gz && \
echo Align "$1" with BWA && bwa mem -M -t 6 '$hg38' fastq/"$28"/Trimmed_SCs/"$1"*_1.fq fastq/"$28"/Trimmed_SCs/"$1"*_2.fq > ./BAM/"$28"/"$1".sam && \
echo Convert "$1" from SAM to BAM && samtools fixmate -@ 4 -O bam ./BAM/"$28"/"$1".sam ./BAM/"$28"/"$1"_unsorted.bam && \
echo Sort "$1" by coordinates && samtools sort -@ 4 -O bam -o ./BAM/"$28"/"$1"_sorted.bam ./BAM/"$28"/"$1"_unsorted.bam  && \
echo Remove duplicates of "$1" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAM/"$28"/"$1"_sorted.bam OUTPUT=BAM/"$28"/"$1".bam METRICS_FILE=Metadata/QC/"$28"/Picard_metrics/"$1"_metrics.out && \
echo Index "$1" && samtools index -@ 4 -b BAM/"$28"/"$1".bam && \
echo Generate alignment statistics "$1" && samtools flagstat -@ 4 BAM/"$28"/"$1".bam > Metadata/QC/"$28"/Alignment_stats/"$1"_stats.out && \
rm fastq/"$28"/Trimmed_SCs/"$1"*.fq && rm ./BAM/"$28"/"$1".sam && rm ./BAM/"$28"/"$1"_sorted.bam && rm ./BAM/"$28"/"$1"_unsorted.bam ")} } }' Metadata/Minussi2021_data_summary.csv

#Split single and paired end MDAMB231_popp31 cells
mkdir -p BAM/MDAMB231_popp31_Single BAM/MDAMB231_popp31_Paired
cat Metadata/Minussi2021_data_summary.csv | grep MDAMB231_popp31 | awk -F';' '{ if ($4 != "202") { if ($20 == "SINGLE") {system(" mv BAM/MDAMB231_popp31/"$1"* BAM/MDAMB231_popp31_Single/ ")} else {system(" mv BAM/MDAMB231_popp31/"$1"* BAM/MDAMB231_popp31_Paired/ ")} } }'
rm -r BAM/MDAMB231_popp31

################################################
#################### KRONOS ####################
################################################
# CNV
cd Minussi2021
mkdir -p Kronos/CNV

#Paired
for celltype in $(echo MDAMB231_popp31_Paired TN1_10x TN3_10x) ; do
  date ; echo $celltype
  Kronos CNV -D BAM/$celltype -B ../Binned_Genomes/hg38/25kPE/*.tsv -c 8 -o Kronos/CNV -M 10 -e $celltype -n 100000 --chr_range 1:22,X,M && echo done
done
#Single
for celltype in $(echo BT20 MDAMB231_popp31_Single MDAMB231c28 MDAMB231c8 TN1 TN2 TN3 TN4 TN5 TN6 TN7 TN8 mb157 mb453) ; do
  date ; echo $celltype
  Kronos CNV -D BAM/$celltype -B ../Binned_Genomes/hg38/25kSE/*.tsv -c 6 -o Kronos/CNV -M 10 -e $celltype -n 100000 --chr_range 1:22,X,M && echo done
done

## Cell stats
sh calculate_stats.sh > Metadata/Minussi_read_stats.csv

##############
# Diagnostic # For MnM
##############
datasetcount="$(ls ./BAM | wc -l)"
tail -n+2 Metadata/Minussi_processing_for_Kronos_and_MnM.txt | awk -F'\t' '{ system(" date ; echo "$1": dataset "NR" of '$datasetcount' ; Kronos diagnostic -F Kronos/CNV/"$2" -o Kronos/Diagnostic/Dev -b "$1" -m "$10" -d ") }'

#################
##### M n M #####
#################
mkdir Kronos/CNV/merged
# Metadata (cell name, cell type)
for bed in Kronos/CNV/*_per_Cell_summary_metrics.csv ; do
	name=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' | awk -F'_SRR' '{print $1 }')
	tail -n+2 Kronos/CNV/"$name"_per_Cell_summary_metrics.csv | awk -F',' -v n="$name" 'OFS="\t"{print $1,n}'
done > Metadata/All_Minussi_CellsAndTypes.tsv

# Merge
## MDA-MB-231
(cat Kronos/CNV/MDAMB231_popp31_Paired_cnv_calls.bed ; tail -n+2 Kronos/CNV/MDAMB231_popp31_Single_cnv_calls.bed ; tail -n+2 Kronos/CNV/MDAMB231c28_cnv_calls.bed ; tail -n+2 Kronos/CNV/MDAMB231c8_cnv_calls.bed ) > Kronos/CNV/merged/MDA-MB-231.bed
## TN1
(cat Kronos/CNV/TN1_10x_cnv_calls.bed ; tail -n+2 Kronos/CNV/TN1_cnv_calls.bed ) > Kronos/CNV/merged/TN1.bed
## TN3
(cat Kronos/CNV/TN3_10x_cnv_calls.bed ; tail -n+2 Kronos/CNV/TN3_cnv_calls.bed ) > Kronos/CNV/merged/TN3.bed

## MnM ##
for bed in $(awk -F'\t' '{print $11}' Metadata/Minussi_processing_for_Kronos_and_MnM.txt | tail -n+2 | sort | uniq | grep -E 'TN1|TN3') ; do
	celltype=$(echo $bed | sed 's/Kronos\/CNV\///' | sed 's/.bed//' | sed 's/merged\///' | sed 's/_cnv_calls//' | sed 's/MDAMB231_popp31_Paired/MDA-MB-231/' | sed 's/MDAMB231_popp31_Single/MDA-MB-231/' |  sed 's/MDAMB231c28/MDA-MB-231/' |  sed 's/MDAMB231c8/MDA-MB-231/' | sed 's/mb157/MDA-MB-157/' | sed 's/mb453/MDA-MB-453/' )
	echo $celltype
	MnM -i $bed -o MnM/$celltype -n $celltype -r -s -b --seed 18671107 --groups Metadata/All_Minussi_CellsAndTypes.tsv
done

##############
# Diagnostic # For Kronos this time
##############
# WhoIsWho
for csv in Kronos/CNV/*_per_Cell_summary_metrics.csv ; do
	name=$(echo $csv | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' | awk -F'_SRR' '{print $1 }')
	celltype=$(echo $csv | sed 's/Kronos\/CNV\///' | sed 's/_per_Cell_summary_metrics.csv//g' | sed 's/_10x//' | sed 's/MDAMB231_popp31_Paired/MDA-MB-231/' | sed 's/MDAMB231_popp31_Single/MDA-MB-231/' |  sed 's/MDAMB231c28/MDA-MB-231/' |  sed 's/MDAMB231c8/MDA-MB-231/' | sed 's/mb157/MDA-MB-157/' | sed 's/mb453/MDA-MB-453/' )
	echo $celltype
	Kronos WhoIsWho -F $csv -o Kronos/Diagnostic/Whoiswho -W MnM/$celltype/*_metadata.tsv
done

# Diagnostic for Kronos
tail -n+2 Metadata/Minussi_processing_for_Kronos_and_MnM.txt | awk -F'\t' '{ system(" date ; echo "$1": dataset "NR" ; Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_"$2" -o Kronos/Diagnostic -b "$1" -m "$10" -C && echo ; echo ; echo  ") }'
# fix TN1,TN1_10x,TN2,TN3,TN3_10x,TN6 that didn't work
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_TN1_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b TN1 -m 130 -C -f 0.95 -s 0.55
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_TN1_10x_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b TN1_10x -m 80 -C -f 0.95 -s 0.55
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_TN2_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b TN2 -m 150 -C -f 0.95 -s 0.55
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_TN3_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b TN3 -m 70 -C -f 0.95 -s 0.55
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_TN3_10x_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b TN3_10x -m 70 -C -f 0.95 -s 0.55
Kronos diagnostic -F Kronos/Diagnostic/Whoiswho/phased_TN6_per_Cell_summary_metrics.csv -o Kronos/Diagnostic -b TN6 -m 75 -C -f 0.95 -s 0.55

# Create file to split Subpopulations/Reps
for dataset in $(ls MnM) ; do
	name=$(echo $dataset | sed 's/MnM\///' )
	echo $name
  mkdir -p Kronos/subpopulation_files_split/$name
  number=$(echo $name | sed 's/MnM\///' | sed 's/MDA-MB-//' )
	python split_subpopulations_csv_bed.py --output-dir Kronos/subpopulation_files_split/$name --csv-files $(ls Kronos/Diagnostic/Whoiswho/phased_*"$number"*_per_Cell_summary_metrics.csv) --bed-files $(ls MnM/$dataset/*.BED.gz)
done

######
# RT #
######
cores=6
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
	-b $(paste -d '_' <(echo "$library_id" ) <(echo "$subpop") | sed 's/^_//g' | tr '\n' ',' | sed 's/,$//' ) \
	-f 200kb \
	-g $(paste -d '_' <(echo "$cell_type") <(echo "$sample_id") <(echo "$subpop") | tr '\n' ',' | sed 's/,$//' | sed 's/__/_/g' ) \
	-S $(paste -d '\0' <(yes "Kronos/Diagnostic/" | head -n $(echo $data | wc -l)) <(echo "$data") <(yes "_settings.txt" | head -n $(echo $data | wc -l))  | tr '\n' ',' | sed 's/,$//') \
	-B 200000 \
	-c $cores \
	-N 5 \
	-p \
	--extract_G1_G2_cells \
	--chr_range 1:22,X
done
