#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

#declare paths
hg38=#path/WholeGenomeFasta/genome.fa
picard=#path/picard.jar

#################################################
################# DOWNLOAD DATA #################
#################################################

#Download dataset with pyega3
#awk -F "," 'NR > 1 {system(" pyega3 -c 30 fetch -M 50 --saveto ./BAMs_hg19/"$1"/ "$1" ")}' Metadata/Laks_Datasets.csv

#remove unmerged slices
for dataset in ./EGAD*[0-9^] ; do
  printf "${dataset##*/}\t"
  rm $dataset/*/*.slice
done

#find downloaded dataset sizes
du -shm ./BAMs_hg19/*

#Number of downloaded BAM/CRAMS in datasets from EGA
for dataset in ./BAMs_hg19/* ; do
  printf "${dataset##*/}\t"
  ls $dataset/*/*am | wc -l
done
#Number of downloaded fastq in datasets from EGA
for dataset in ./BAMs_hg19/* ; do
  printf "${dataset##*/}\t"
  ls $dataset/*/*fastq.gz | wc -l
done

########################################################
################# CONVERT HG19 -> HG38 #################
########################################################

# Create FASTQ files from BAMs (2 files per BAM because paired ends)
#create directories
mkdir Laks2019
mkdir Laks2019/fastq
mkdir Laks2019/Metadata
mkdir Laks2019/Trimmed
mkdir Laks2019/BAMs_hg38

#Move fastq files to new destination
mv Laks2019/BAMs_hg19/EGAD00001007756/*/*fastq.gz Laks2019/fastq/EGAD00001007756/
ls -R Laks2019/fastq/EGAD00001007756/ | grep '\.fastq.gz$' | sed -e 's/\.concat.fastq.gz$//' > Laks2019/Metadata/EGAD00001007756_filename_list.txt

#launch BAM to FASTQ conversion with progress information (prints: File n of TOTAL). Input is SORTED BAM file.
datasetcount="$(ls $data3/BAMs_hg19/ | wc -l)"
filecount="$(ls -R $data3/BAMs_hg19/*/ | grep '\.bam$' | wc -l)"
i=1
j=1
for dataset in $data3/BAMs_hg19/* ; do
 echo "${dataset##*/}": DATASET $j of $datasetcount && ((j=j+1))
 mkdir Laks2019/fastq/"$(basename "${dataset%.*}")"
 ls -R $data3/BAMs_hg19/"$(basename "${dataset%.*}")"/ | grep '\.bam$' | sed -e 's/\..*$//' > Laks2019/Metadata/"$(basename "${dataset%.*}")"_filename_list.txt
 for file in $dataset/*/*.bam ; do
   echo "${file##*/}": File $i of $filecount && ((i=i+1))
   samtools fastq -@ 5 -1 Laks2019/fastq/"$(basename "${dataset%.*}")"/"$(basename "${file%.*}")"_1.fastq.gz -2 Laks2019/fastq/"$(basename "${dataset%.*}")"/"$(basename "${file%.*}")"_2.fastq.gz <(samtools sort -n -@ 5 $file)
 done
done

#Trim adaptors from reads
cd Laks2019/
datasetcount="$(ls ./fastq | wc -l)"
j=1
for dataset in ./fastq/* ; do
 echo "${dataset##*/}": DATASET $j of $datasetcount && ((j=j+1))
 mkdir ./Metadata/"${dataset##*/}"_QC
 awk -F " " ' { system("trim_galore --paired ./fastq/'${dataset##*/}'/"$1"_1.fastq.gz ./fastq/'${dataset##*/}'/"$1"_2.fastq.gz --output_dir ./Trimmed/'${dataset##*/}' --fastqc --cores 2 --gzip") }' ./Metadata/"${dataset##*/}"_filename_list.txt && (multiqc -o ./Metadata/"${dataset##*/}"_QC ./Trimmed/"${dataset##*/}" & )
done

#index reference genome (only to be done once)
mkdir ../BWA_Index
mkdir ../BWA_Index/hg38
cd ../BWA_Index/hg38
bwa index "$hg38".fa

# Align each set of FASTQ files against reference genome using BWA
cd Laks2019/
gzip fastq/*/*.fastq
datasetcount="$(ls ./Trimmed | wc -l)"
j=1
for dataset in ./Trimmed/* ; do
 echo "${dataset##*/}": DATASET $j of $datasetcount && ((j=j+1))
 mkdir ./BAMs_hg38/"${dataset##*/}"
 mkdir ./Metadata/"${dataset##*/}"_QC/Picard_metrics
 mkdir ./Metadata/"${dataset##*/}"_QC/Alignment_stats
 awk -F " " ' { system("echo Decompress "$1" && gunzip ./Trimmed/'${dataset##*/}'/"$1"*.gz && echo Align "$1" with BWA && bwa mem -M -t 4 '$hg38' ./Trimmed/'${dataset##*/}'/"$1"_1_val_1.fq ./Trimmed/'${dataset##*/}'/"$1"_2_val_2.fq > ./BAMs_hg38/'${dataset##*/}'/"$1".sam && echo Convert "$1" from SAM to BAM && samtools fixmate -O bam ./BAMs_hg38/'${dataset##*/}'/"$1".sam ./BAMs_hg38/'${dataset##*/}'/"$1"_unsorted.bam && echo Sort "$1" by coordinates && samtools sort -O bam -@ 4 -o ./BAMs_hg38/'${dataset##*/}'/"$1"_sorted.bam ./BAMs_hg38/'${dataset##*/}'/"$1"_unsorted.bam && echo Remove duplicates of "$1" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAMs_hg38/'${dataset##*/}'/"$1"_sorted.bam OUTPUT=BAMs_hg38/'${dataset##*/}'/"$1".bam METRICS_FILE=Metadata/'${dataset##*/}'_QC/Picard_metrics/"$1"_metrics.out && echo Index "$1" && samtools index -b BAMs_hg38/'${dataset##*/}'/"$1".bam && echo Generate alignment statistics "$1" && samtools flagstat BAMs_hg38/'${dataset##*/}'/"$1".bam > Metadata/'${dataset##*/}'_QC/Alignment_stats/"$1"_stats.out && rm ./Trimmed/'${dataset##*/}'/"$1"* && rm ./BAMs_hg38/'${dataset##*/}'/"$1".sam && rm ./BAMs_hg38/'${dataset##*/}'/"$1"_unsorted.bam && rm ./BAMs_hg38/'${dataset##*/}'/"$1"_sorted.bam ") }' ./Metadata/"${dataset##*/}"_filename_list.txt
done

# View fastq read names
awk 'NR%4==1 {print substr($1,2)}' ./Trimmed/EGAD00001004764_test/SA1087-A96150A-R16-C09*fq | head

# Rename BAM Files to give sample names
cd Laks2019/
awk -F';' 'NR>1{ system("echo "$4"_"$5"_"$1" ; mkdir BAMs_hg38/"$4"_"$5"_"$1" ; mv "$3"* BAMs_hg38/"$4"_"$5"_"$1"/") }' Metadata/Laks2019_cell_descriptions.csv
# Verify all EGAD* folders are empty
ls BAMs_hg38/EGAD*/
# Delete all EGAD* folders which are now empty
rm -r BAMs_hg38/EGAD*/
