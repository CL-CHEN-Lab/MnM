#!/bin/sh
# Author: JOSEPH JOSEPHIDES
# Institut Curie, Paris

cd Connolly2022
mkdir Metadata fastq
cd Metadata
curl -O https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10234/E-MTAB-10234.sdrf.txt
cd ../fastq

awk -F"\t" 'NR>1{ system(" wget ht$(cut -c2- <<< "$36") -O "$34" & ") }' ./Metadata/E-MTAB-10234.sdrf.txt

gunzip *.gz

#treated as single-end
###### Trim adaptors from reads
cd Connolly2022
fastqcount="$(ls ./fastq/*.fastq | wc -l)"
j=1
mkdir ./Metadata/QC
for fastq in ./fastq/*.fastq ; do
 echo "${fastq##*/}": fastq $j of $fastqcount && ((j=j+1))
 trim_galore $fastq --output_dir ./Trimmed/ --fastqc --cores 4
done
multiqc -o ./Metadata/QC ./Trimmed/

#### Mapping
hg38=#path/WholeGenomeFasta/genome.fa
picard=#path/picard.jar

cd Connolly2022
fastqcount="$(ls ./Trimmed/*.fq | wc -l)"
j=1
mkdir ./BAMs
mkdir ./Metadata/QC/Picard_metrics
mkdir ./Metadata/QC/Alignment_stats
for fastq in ./Trimmed/*.fq ; do
  echo ; echo
  echo "${fastq##*/}": FASTQ $j of $fastqcount && ((j=j+1))
  echo $fastq && bwa mem -M -t 4 $hg38 $fastq > ./BAMs/$(basename "${fastq%.*}").sam && echo Convert SAM to BAM and sort by coordinates && samtools sort -O bam -@ 4 -o ./BAMs/$(basename "${fastq%.*}")_sorted.bam ./BAMs/$(basename "${fastq%.*}").sam && echo Remove duplicates && java -Xmx16g -jar $picard MarkDuplicates ASSUME_SORT_ORDER=coordinate INPUT=BAMs/$(basename "${fastq%.*}")_sorted.bam OUTPUT=BAMs/$(basename "${fastq%.*}").bam METRICS_FILE=Metadata/QC/Picard_metrics/$(basename "${fastq%.*}")_metrics.out && echo Index && samtools index -b BAMs/$(basename "${fastq%.*}").bam && echo Generate alignment statistics && samtools flagstat BAMs/$(basename "${fastq%.*}").bam > Metadata/QC/Alignment_stats/$(basename "${fastq%.*}")_stats.out && rm ./BAMs/$(basename "${fastq%.*}").sam && rm ./BAMs/$(basename "${fastq%.*}")_sorted.bam && rm Trimmed/$(basename "${fastq%.*}")*
done

rm -r Trimmed
pigz -p 4  fastq/*.fastq

#### ALL are single END
#### KRONOS ####
hg38BT2=#path/Bowtie2Index/genome
hg38Ref=#path/WholeGenomeFasta/genome.fa
hg38Blist=#path/hg38-blacklist.v2.bed
Picard=#path/picard.jar
cores=6
#######
# CNV #
#######
cd Connolly2022
## Single-end
Kronos CNV -D BAMs -B ../Binned_Genomes/hg38/25kSE/*.tsv -c $cores -o Kronos/CNV -e hTERT-RPE1 -n 100000 --chr_range 1:22,X

############
# WhoIsWho #
############
Kronos WhoIsWho -F Kronos/CNV/hTERT-RPE1_per_Cell_summary_metrics.csv -W Metadata/Connolly2022_bam_phase.txt -o Kronos/WhoIsWho

##############
# Diagnostic #
##############
Kronos diagnostic -F Kronos/WhoIsWho/phased_hTERT-RPE1_per_Cell_summary_metrics.csv -o Kronos/diagnostic -b hTERT-RPE1 -c 6 -m 70 -f 0.95 -s 0.5 -C

#######
# MnM #
#######
MnM -i Kronos/CNV/hTERT-RPE1_cnv_calls.bed -o MnM -n hTERT-RPE1 -r -b --seed 18671107

######
# RT #
######
#Calculates scReplication profiles and scRT
Kronos RT -F Kronos/WhoIsWho/phased_hTERT-RPE1_per_Cell_summary_metrics.csv -T Kronos/CNV/hTERT-RPE1_cnv_calls.bed -C ../genome_assemblies/HUMAN/name_length_hg38.txt -o Kronos/RT -b hTERT-RPE1 -f 250kb -S Kronos/diagnostic/hTERT-RPE1_settings.txt -c 6 -p -r chr2:20Mb-70Mb -B 250000 -N 3 -R Takahashi2019/HUMAN/Metadata/Reference_RT/hTERT-RPE1_hg38_Reference.bed

# RT with Connolly2022 + Takahashi2019 data
Kronos RT -F Kronos/WhoIsWho/phased_hTERT-RPE1_per_Cell_summary_metrics.csv,Takahashi2019/HUMAN/Kronos/WhoIsWho/phased_hTERT-RPE1_per_Cell_summary_metrics.csv -T Kronos/CNV/hTERT-RPE1_cnv_calls.bed,Takahashi2019/HUMAN/Kronos/CNV/hTERT-RPE1_cnv_calls.bed -C ../genome_assemblies/HUMAN/name_length_hg38.txt -o Kronos/RT_Connolly_Takahashi -b RPE1_Connolly2022,RPE1_Takahashi2019 -f 250kb -S Kronos/diagnostic/hTERT-RPE1_settings.txt,Takahashi2019/HUMAN/Kronos/diagnostic/hTERT-RPE1_settings.txt -c 6 -p -r chr2:20Mb-70Mb -B 250000 -N 3 -R Takahashi2019/HUMAN/Metadata/Reference_RT/hTERT-RPE1_hg38_Reference.bed
