#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

cd Du2021
mkdir SRA fastq Metadata
cd SRA
# download SRA files with prefetch
awk 'NR>1 { system(" echo ; echo "$1" ; prefetch -C yes -f no -v -p -X 500G -O ./ "$1" &") }' ./Metadata/SraRunTable.txt
# convert SRA files to R1 and R2 fastq files, remove technical R3 file and compress (gz).
awk 'NR>1 { system(" echo ; echo "$1" ; fasterq-dump --include-technical -S -x -p -e 20 -O ../fastq "$1"/"$1".sra && rm ../fastq/"$1"_3* && pigz ../fastq/"$1"*.fastq && rm -r "$1" ") }' ./Metadata/SraRunTable.txt

# move into named directories
cd Du2021/fastq
mkdir HCT116_WT_scRepli-Seq_G1 HCT116_WT_scRepli-Seq_S HCT116_DKO1_scRepli-Seq_G1 HCT116_DKO1_scRepli-Seq_S
mv SRR12646263* HCT116_WT_scRepli-Seq_G1/
mv SRR12646264* HCT116_WT_scRepli-Seq_S/
mv SRR12646265* HCT116_DKO1_scRepli-Seq_G1/
mv SRR12646266* HCT116_DKO1_scRepli-Seq_S/

# get 10x barcode whitelist, add "-1" as name to barcode and make 2nd column the actual sequence of the barcode. save as tsv.
cd ./Metadata
wget https://raw.githubusercontent.com/TheKorenLab/Single-cell-replication-timing/main/align/10x_barcode_whitelist.txt
cat 10x_barcode_whitelist.txt | awk '{print $1"-1",$1}' > 10x_barcode_whitelist.tsv

# demultiplex
cd ../fastq
for DIR in * ; do
	echo Demultiplex $DIR
	demultiplex demux -r -e 16 -p $DIR/SCs ./Metadata/10x_barcode_whitelist.tsv $DIR/*_1.fastq.gz $DIR/*_2.fastq.gz
done

#find number of reads to remove false positives
for fastq in */SCs/* ; do printf $fastq ; gunzip -c $fastq | wc -l  ; done > ./Metadata/wc-l.txt

# Number of lines in each fastq file (1 read = 4 lines)
wc -l ./Metadata/wc-l.txt ; echo ; cat ./Metadata/wc-l.txt | sort -nr -k2 | head
# Convert file to make easier to read in R
cd ./Metadata
cat wc-l.txt | sed 's/  */\t/g' | awk ' OFS="\t" {print $2,$1}'  > wc-l_tab.txt

# Run following script to create files to faciliate cell filtering and mapping
Rscript ../SCutoff_v1.R -f wc-l_tab.txt -n Du2021 -o ./
# Note: In the publication, they removed barcodes with less than ~1M reads (Methods; Table S6). Not the same approach here.

# Create new directories for filtered cells to use
cd ../fastq
for dir in ./* ; do mkdir $dir/Retained_SCs ; done

# Move filtered cells to new directory
awk 'NR>1 { system(" printf "$5"\,\   ; mv "$2" "$6" ") }' ./Metadata/2022.10.17_CellsToRetain_Du2021.tsv

#Trim adaptors from reads
filecount="$(cat ./Metadata/2022.10.17_PairedMatches_Du2021.tsv | wc -l)"
awk ' { system(" echo ; date ; echo "$7": file "NR" of '$filecount' ; trim_galore --paired "$1" "$2" --output_dir "$3" --clip_R1 16 --fastqc --cores 6 --gzip") }' ./Metadata/2022.10.17_PairedMatches_Du2021.tsv

#Quality control
for dir in * ; do mkdir -p ./Metadata/QC/$dir ; multiqc -o ./Metadata/QC/$dir $dir/Trimmed_SCs/ ; done

#remove false barcode reads
rm ./*/SCs/*UNKNOWN*

# create required directories for mapping
cd Du2021
mkdir BAM
mkdir $(cat ./Metadata/2022.10.17_PairedMatches_Du2021.tsv | awk '{print "BAM/"$6}' | sort | uniq)
mkdir -p $(cat ./Metadata/2022.10.17_PairedMatches_Du2021.tsv | awk '{print "Metadata/QC/"$6"/Picard_metrics"}' | sort | uniq)
mkdir -p $(cat ./Metadata/2022.10.17_PairedMatches_Du2021.tsv | awk '{print "Metadata/QC/"$6"/Alignment_stats"}' | sort | uniq)

#declare paths
hg38=../reference_genomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa
picard=#path/picard.jar

### Mapping reads to create BAM file ###
cd Du2021
#gzip fastq/*/Trimmed_SCs/*.fq
filecount="$(cat Metadata/2022.10.17_PairedMatches_Du2021.tsv | wc -l)"
awk ' { system(" date ; echo Decompress "$7" : file "NR" of '$filecount' ; ((i=i+1)) && unpigz fastq/"$4".gz fastq/"$5".gz ; echo Align "$7" with BWA && bwa mem -M -t 6 '$hg38' fastq/"$4" fastq/"$5" > ./BAM/"$6"/"$7"-1.sam && echo Convert "$7" from SAM to BAM && samtools fixmate -@ 6 -O bam ./BAM/"$6"/"$7"-1.sam ./BAM/"$6"/"$7"-1_unsorted.bam && echo Sort "$7" by coordinates && samtools sort -@ 6 -O bam -@ 4 -o ./BAM/"$6"/"$7"-1_sorted.bam ./BAM/"$6"/"$7"-1_unsorted.bam && echo Remove duplicates of "$7" && java -Xmx16g -jar '$picard' MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAM/"$6"/"$7"-1_sorted.bam OUTPUT=BAM/"$6"/"$7"-1.bam METRICS_FILE=Metadata/QC/"$6"/Picard_metrics/"$7"_metrics.out && echo Index "$7" && samtools index -@ 6 -b BAM/"$6"/"$7"-1.bam && echo Generate alignment statistics "$7" && samtools flagstat -@ 6 BAM/"$6"/"$7"-1.bam > Metadata/QC/"$6"/Alignment_stats/"$7"_stats.out && rm fastq/"$4" fastq/"$5" && rm ./BAM/"$6"/"$7"-1.sam && rm ./BAM/"$6"/"$7"-1_unsorted.bam && rm ./BAM/"$6"/"$7"-1_sorted.bam && echo Done "$7" ; echo ; echo ") }' ./Metadata/2022.10.17_PairedMatches_Du2021.tsv
