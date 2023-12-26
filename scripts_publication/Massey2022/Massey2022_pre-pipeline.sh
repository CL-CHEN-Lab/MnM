#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

#################################
########### SRA FILES ###########
#################################

## Contains:
## GM12878 (unsorted), GM12891, GM12892, H1, H7, H9, HCT-116, MCF-7, RKO.
## Note: GM12878 are downloaded as BAM files but also tested as SRA files (21 July 2022).
## Note 2: BAM files preffered for sorted GM12878 files because their respective SRA files no longer contain the barcode in R1 and therefore impossible to demultiplex (30 Aug 2022).

# Make directories
mkdir Massey2022
cd Massey2022
mkdir SRA fastq Metadata

#Download metadata from https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA770772 and save to Metadata directory
cd SRA

#Download data as SRA files
awk -F ";" 'NR>1 { system(" echo "$1" ; prefetch -C yes -f no -T sra -v -p -X 900G -O ./ "$1" &") }' ./Metadata/PRJNA770772.csv

# convert SRA to fastq
# fasterq-dump not working on NAS (maybe file overload). Moving SRA files to HD as a temporary processing location.
mv SRR* SRA/
awk -F ";" 'NR>1 { system(" echo "$1" ; fasterq-dump --include-technical -S -x -p -e 5 -O ../fastq -t SRA SRA/"$1"/"$1".sra && pigz ../fastq/"$1"*.fastq && rm -r SRA/"$1" ") }' ./Metadata/PRJNA770772.csv

# create subdirectories, one for each cell population with descriptive names (e.g. SRR16328242_MCF7-unsorted-lane2)
awk -F ";" 'NR>1 && NR<22 { system(" echo "$1"_"$18" ; mkdir ../fastq/"$1"_"$18" ; mv ../fastq/"$1"*.gz ../fastq/"$1"_"$18"/ ") }' ./Metadata/PRJNA770772.csv

# get 10x barcode whitelist
cd Massey2022/Metadata
wget https://raw.githubusercontent.com/TheKorenLab/Single-cell-replication-timing/main/align/10x_barcode_whitelist.txt
cat 10x_barcode_whitelist.txt | awk '{print $1"-1",$1}' > 10x_barcode_whitelist.tsv

# demultiplex by barcode (16bp from R1; R2 also sorted so carried on from R1 barcode)
cd Massey2022/fastq
awk -F ";" 'NR>1 && NR<22  { system(" echo "$1"_"$18" ; pigz -d "$1"_"$18"/"$1"_*.fastq.gz ; echo Demultiplexing "$1"_"$18" ; demultiplex demux -r -e 16 -p ./"$1"_"$18"/SCs ./Metadata/10x_barcode_whitelist.tsv "$1"_"$18"/"$1"_1.fastq "$1"_"$18"/"$1"_2.fastq && pigz "$1"_"$18"/"$1"_*.fastq ") }' ./Metadata/PRJNA770772.csv


#################################
########### BAM FILES ###########
#################################
#Download BAM files from SRA

# sort BAM files by read name and convert to fastq with cell-specific barcode in header (CB field).
cd Massey2022/BAMs
for bam in *.bam ; do
	echo $(basename "${bam%.*}")
	samtools sort -m 5G -n -@ 5 $bam -o $(basename "${bam%.*}")_sorted.bam && samtools fastq -@ 5 -T CB --barcode-tag CB -1 ../fastq/$(basename "${bam%.*}")_1.fastq.gz -2 ../fastq/$(basename "${bam%.*}")_2.fastq.gz $(basename "${bam%.*}")_sorted.bam && rm $(basename "${bam%.*}")_sorted.bam
done

# get 10x barcode whitelist
cd Massey2022/Metadata
wget https://raw.githubusercontent.com/TheKorenLab/Single-cell-replication-timing/main/align/10x_barcode_whitelist.txt
cat 10x_barcode_whitelist.txt | awk '{print $1"-1",$1"-1"}' > 10x_barcode_whitelist-1.tsv

# move fastq files to directories and demultiplex fastq files
cd Massey2022/BAMs
for bam in *.bam ; do
	echo ; echo $(basename "${bam%.*}")
	mkdir ../fastq/$(basename "${bam%.*}")/ ; mv ../fastq/$(basename "${bam%.*}")*.fastq.gz ../fastq/$(basename "${bam%.*}")/
	pigz -d ../fastq/$(basename "${bam%.*}")/$(basename "${bam%.*}")_*.fastq.gz
	echo Demultiplex $(basename "${bam%.*}")
	demultiplex demux -p ../fastq/$(basename "${bam%.*}")/SCs -m 0 --format x ./Metadata/10x_barcode_whitelist-1.tsv ../fastq/$(basename "${bam%.*}")/$(basename "${bam%.*}")_1.fastq ../fastq/$(basename "${bam%.*}")/$(basename "${bam%.*}")_2.fastq
	echo Done
	pigz ../fastq/$(basename "${bam%.*}")/$(basename "${bam%.*}")_*.fastq
done

cd Massey2022
rm Massey2022/BAMs

#####################################################
########### All Demultiplexed Fastq files ###########
#####################################################
# Remove invalid barcodes
rm */SCs/*UNKNOWN*

# Record number of lines for each Demultiplexed fastq file to remove false positives (valid barcode but not enough reads).
cd Massey2022/fastq
for fastq in GM*/SCs/* ; do wc -l $fastq ; done > ./Metadata/wc-l_GM12878_BAM.txt
for celltype in SRR* ; do
	for fastq in $celltype/SCs/* ; do wc -l $fastq ; done
done > ./Metadata/wc-l_SRR_fastq.txt

# Convert file to make easier to read in R
cat ./Metadata/wc-l_GM12878_BAM.txt | sed 's/  */\t/g' | awk ' OFS="\t" {print $1,$2}'  > ./Metadata/wc-l_GM12878_BAM_tab.txt
cat ./Metadata/wc-l_SRR_fastq.txt | sed 's/  */\t/g' | awk ' OFS="\t" {print $1,$2}'  > ./Metadata/wc-l_SRR_fastq_tab.txt

# Run following script to create files to faciliate cell filtering and mapping
Rscript ../SCutoff_v1.R -f ./Metadata/wc-l_GM12878_BAM_tab.txt -n Massey2022_BAM -o ./Metadata
Rscript ../SCutoff_v1.R -f ./Metadata/wc-l_SRR_fastq_tab.txt -n Massey2022_SRA -o ./Metadata

# Create new directories for filtered cells to use
cd ../fastq
for dir in ./* ; do mkdir $dir/Retained_SCs ; done

# Move filtered cells to new directory
awk 'NR>1 { system(" printf "$5"\,\   ; mv "$2" "$6" ") }' ./Metadata/2022.10.14_CellsToRetain_Massey2022_BAM.tsv
awk 'NR>1 { system(" printf "$5"\,\   ; mv "$2" "$6" ") }' ./Metadata/2022.11.21_CellsToRetain_Massey2022_SRA.tsv

#Trim adaptors from reads
filecount="$(cat ./Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv | wc -l)"
awk ' { system(" echo ; date ; echo "$7": file "NR" of '$filecount' ; trim_galore --paired "$1" "$2" --output_dir "$3" --fastqc --cores 8 --gzip && echo Compressing "$7" ; pigz "$1" "$2" ") }' ./Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv
filecount="$(cat ./Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv | wc -l)"
awk ' { system(" echo ; date ; echo "$7": file "NR" of '$filecount' ; trim_galore --paired "$1" "$2" --output_dir "$3" --clip_R1 16 --fastqc --cores 4 --gzip && echo Compressing "$7" ; pigz "$1" "$2" ") }' ./Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv

#Quality control
for dir in GM* ; do mkdir -p ./Metadata/QC/$dir ; multiqc -o ./Metadata/QC/$dir $dir/Trimmed_SCs/ ; done
for dir in SRR* ; do mkdir -p ./Metadata/QC/$dir ; multiqc -o ./Metadata/QC/$dir $dir/Trimmed_SCs/ ; done

# create required directories for mapping
cd ..
mkdir BAM
mkdir $(cat Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv | awk '{print "BAM/"$6}' | sort | uniq)
mkdir $(cat Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv | awk '{print "BAM/"$6}' | sort | uniq)
mkdir -p $(cat ./Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv | awk '{print "Metadata/QC/"$6"/Picard_metrics"}' | sort | uniq)
mkdir -p $(cat ./Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv | awk '{print "Metadata/QC/"$6"/Alignment_stats"}' | sort | uniq)
mkdir -p $(cat ./Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv | awk '{print "Metadata/QC/"$6"/Alignment_stats"}' | sort | uniq)
mkdir -p $(cat ./Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv | awk '{print "Metadata/QC/"$6"/Picard_metrics"}' | sort | uniq)

#declare paths
hg38=#path/BWAIndex/genome.fa
picard=#path/picard.jar

### Mapping reads to create BAM file ###
cd Massey2022
#gzip fastq/*/Trimmed_SCs/*.fq
filecount="$(cat Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv | wc -l)"
awk ' { system(" date ; echo Decompress "$7" : file "NR" of '$filecount' ; unpigz fastq/"$4".gz fastq/"$5".gz && echo Align "$7" with BWA && bwa mem -M -t 7 '$hg38' fastq/"$4" fastq/"$5" > ./BAM/"$6"/"$7"-1.sam && echo Convert "$7" from SAM to BAM && samtools fixmate -@ 6 -O bam ./BAM/"$6"/"$7"-1.sam ./BAM/"$6"/"$7"-1_unsorted.bam && echo Sort "$7" by coordinates && samtools sort -@ 6 -O bam -o ./BAM/"$6"/"$7"-1_sorted.bam ./BAM/"$6"/"$7"-1_unsorted.bam && echo Remove duplicates of "$7" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAM/"$6"/"$7"-1_sorted.bam OUTPUT=BAM/"$6"/"$7"-1.bam METRICS_FILE=Metadata/QC/"$6"/Picard_metrics/"$7"_metrics.out && echo Index "$7" && samtools index -@ 6 -b BAM/"$6"/"$7"-1.bam && echo Generate alignment statistics "$7" && samtools flagstat -@ 6 BAM/"$6"/"$7"-1.bam > Metadata/QC/"$6"/Alignment_stats/"$7"_stats.out && rm fastq/"$4" fastq/"$5" && rm ./BAM/"$6"/"$7"-1.sam && rm ./BAM/"$6"/"$7"-1_unsorted.bam && rm ./BAM/"$6"/"$7"-1_sorted.bam ") }' Metadata/2022.10.14_PairedMatches_Massey2022_BAM.tsv

# fix file : sed 's/-1.fastq/-1_val_/g' Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv > Metadata/2022.11.21_PairedMatches_Massey2022_SRA_mod.tsv
#gzip fastq/*/Trimmed_SCs/*.fq
filecount="$(cat Metadata/2022.11.21_PairedMatches_Massey2022_SRA.tsv | wc -l)"
awk ' { system(" date ; echo Decompress "$7" : file "NR" of '$filecount' ; unpigz fastq/"$4"1.fq.gz fastq/"$5"2.fq.gz && echo Align "$7" with BWA && bwa mem -M -t 7 '$hg38' fastq/"$4"1.fq fastq/"$5"2.fq > ./BAM/"$6"/"$7"-1.sam && echo Convert "$7" from SAM to BAM && samtools fixmate -@ 6 -O bam ./BAM/"$6"/"$7"-1.sam ./BAM/"$6"/"$7"-1_unsorted.bam && echo Sort "$7" by coordinates && samtools sort -@ 6 -O bam -o ./BAM/"$6"/"$7"-1_sorted.bam ./BAM/"$6"/"$7"-1_unsorted.bam && echo Remove duplicates of "$7" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAM/"$6"/"$7"-1_sorted.bam OUTPUT=BAM/"$6"/"$7"-1.bam METRICS_FILE=Metadata/QC/"$6"/Picard_metrics/"$7"_metrics.out && echo Index "$7" && samtools index -@ 6 -b BAM/"$6"/"$7"-1.bam && echo Generate alignment statistics "$7" && samtools flagstat -@ 6 BAM/"$6"/"$7"-1.bam > Metadata/QC/"$6"/Alignment_stats/"$7"_stats.out && rm fastq/"$4"* fastq/"$5"* && rm ./BAM/"$6"/"$7"-1.sam && rm ./BAM/"$6"/"$7"-1_unsorted.bam && rm ./BAM/"$6"/"$7"-1_sorted.bam ") }' Metadata/2022.11.21_PairedMatches_Massey2022_SRA_mod.tsv

#remove unused cells
cd Massey2022
rm -r fastq/*/SCs
