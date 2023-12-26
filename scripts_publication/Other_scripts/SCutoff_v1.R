#!/usr/bin/env Rscript
# Author: Joseph JOSEPHIDES, Institut Curie, Paris FR
# Created: 2022-09-14

print('Launching SCutoff')
# Check if dependencies are installed and load
pckges1 <- c("stringr", "purrr", "ggplot2", "dplyr", "tidyr","scales","argparser")
new.packages <- pckges1[!(pckges1 %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  print(paste('Installing:', new.packages))
  install.packages(new.packages)
}
pckges2 <- c("devtools","cutoff")
new.packages <- pckges2[!(pckges2 %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  install.packages("devtools")
  devtools::install_github("choisy/cutoff")
}

print("Loading required packages")
invisible(suppressPackageStartupMessages(sapply(c(pckges1,pckges2), require, character.only = TRUE)))

p = arg_parser("SCutoff")
p = add_argument(p, "--file", help="Tabulated line-count file with number of lines and filename.")
p = add_argument(p, "--name", help="Dataset name.", default = "Dataset")
p = add_argument(p, "--output", help="Output location.", default = "./")
p = add_argument(p, "--reads", help="Minimal number of reads per barcode to exlude automatically.", default = 30000, type="integer")
ARGs = parse_args(p)

ARGs$output = str_replace(paste0(ARGs$output , '/'), '//', '/')

# Load data
print('Loading data...')
wc = read.csv(ARGs$file, sep = "\t", header = F, col.names = c('Line_Count', 'File_Path'), )

# Remove any false barcode files
wc = wc[!grepl("UNKNOWN", wc$File_Path),]


# File lists number of lines. 1 read = 4 lines. Make a read count column.
print('Transforming number of lines to number of reads.')
wc$Read_Count = as.numeric(wc$Line_Count / 4)

# Remove certain false positive barcodes.
print('Removing files from list that contain <1000 reads.')
wc2 = wc[wc$Read_Count>=1000,]

# Get a distinct cell type column
print('Extracting cell type(s):')
Cell_Type = unlist(map(strsplit(wc2$File_Path, '/'), 1))
wc2$Cell_Type = Cell_Type
# Make Cell Type column a factor
wc2$Cell_Type = as.factor(wc2$Cell_Type)
print(unique(wc2$Cell_Type))

# Get a distinct cell barcode column
print('Extracting barcodes.')
Barcode = unlist(map(strsplit(wc2$File_Path, '/'), 3))
Barcode = unlist(map(strsplit(Barcode, '-'), 1))
Barcode = str_replace(Barcode, '_R?[12]_', '_')
wc2$Barcode = as.factor(Barcode)

# Sum read count by barcode
print('Merging paired-ends.')
wc3 = wc2 %>%
  group_by(Barcode) %>%
  summarise(Cell_Reads = sum(Read_Count), across()) %>%
  select(Barcode, Cell_Type, Cell_Reads)
wc3 = distinct(wc3)

# Remove barcodes with under 30k reads
print(paste0('Removing barcodes that contain less than ',ARGs$reads,' reads.'))
wc3 = wc3[wc3$Cell_Reads>=(ARGs$reads),]

# Count number of barcodes in each cell type
#table(wc3$Cell_Type)

title = table(wc3$Cell_Type)
title = toString(paste0("n",1:length(title),"=",str_replace_all(format(title, big.mark=','), ' ', '') ))
wc_tmp = wc3
wc_tmp$Cell_Type = str_replace_all(wc_tmp$Cell_Type, '_', '_\n')
boxplot1 = ggplot(wc_tmp, aes(Cell_Type, Cell_Reads)) + geom_boxplot() + theme_linedraw() +
  scale_y_continuous(labels = scales::comma) +
  labs(title = paste0("Reads per cell from ",ARGs$name," data prior filtering"),subtitle = title) +
  xlab("Sample") + ylab("Number of reads")

print('Launching EM algorithm to detect # reads cut-off value for each sample.')

cutoffs = c()
types=c()
wc4=data.frame()
wc3$coff = NA
for (type in sort(unique(wc3$Cell_Type))){
  print(paste0(type,' cutoff estimates:'))
  cell_type_dta = wc3[wc3$Cell_Type == type,]
  em_out <- em(cell_type_dta$Cell_Reads,"log-normal","log-normal")
  cut_off <- cutoff(em_out, distr = 1)
  cutoffs = c(cutoffs, as.integer(cut_off[1]))
  types = c(types, type)
  print(cut_off)
  wc3$coff[wc3$Cell_Type == type] = as.integer(cut_off[1])
  wc4=rbind(wc4,cell_type_dta[cell_type_dta$Cell_Reads>=cut_off[1],])
}
names(cutoffs)=types

histos = ggplot(wc3, aes(x = Cell_Reads)) +
  geom_histogram(aes(color = Cell_Type, fill = Cell_Type),
                 position = "identity", bins = 75, alpha=0.5) +
  facet_wrap(~Cell_Type, ncol = as.integer(length(cutoffs)/2), scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "Read distribution per sample with cut-off values",subtitle = paste0('For barcodes containing ≥', format(ARGs$reads, big.mark=','),' reads') )+
  xlab("Number of reads per cell") + ylab("Number of cells")  + theme_light() + theme(legend.position='none') +
  scale_x_continuous(labels = label_scientific(digits = 0)) + geom_vline(data=wc3 %>% group_by(Cell_Type), aes(xintercept=coff), linetype="longdash")

wc4=subset(wc4, select=-c(coff))

category_n = table(wc4$Cell_Type)
categories = length(category_n)
title2 = toString(paste0("n",1:categories,"=",format(category_n, big.mark=',')))
title2 = paste0(title2, '\nCut-offs: ', toString(paste0("c",1:categories,"=",format(cutoffs, big.mark=',', trim = T)) ))
boxplot2 = ggplot(wc_tmp[wc_tmp$Barcode %in% wc4$Barcode,], aes(Cell_Type, Cell_Reads)) +
  geom_boxplot() + theme_linedraw() + scale_y_continuous(labels = scales::comma) +
  labs(title = paste0("Reads per cell from ",ARGs$name," data after filtering"),subtitle = title2) +
  xlab("Sample") + ylab("Number of reads")

print('Saving plots.')
ggsave(
  paste0(ARGs$output,'/Reads_prior_filtering_',ARGs$name,'.png'),
  plot = boxplot1,
  dpi = 300,
  width = as.integer(length(cutoffs)/2)*4,
  height = ceiling(length(cutoffs)/2)*3,
)
ggsave(
  paste0(ARGs$output,'/Reads_after_filtering_',ARGs$name,'.png'),
  plot = boxplot2,
  dpi = 300,
  width = as.integer(length(cutoffs)/2)*4,
  height = ceiling(length(cutoffs)/2)*3,
)
ggsave(
  paste0(ARGs$output,'/Reads_distribution_with_cutoff_',ARGs$name,'.png'),
  plot = histos,
  dpi = 300,
  width = as.integer(length(cutoffs)/2)*3.5,
  height = ceiling(length(cutoffs)/2)*3,
)

print("Saving files.")

New_Location = str_replace(wc2$File_Path, '/SCs/', '/Retained_SCs/')
wc2$New_Location = New_Location
filetokeep = wc2[wc2$Barcode %in% wc4$Barcode,]

write.table(x = filetokeep, file = paste0(ARGs$output,format(Sys.time(), "%Y.%m.%d"),'_CellsToRetain_',ARGs$name,'.tsv'), sep = '\t', quote = F, col.names = T, row.names = F)

paired_ends = filetokeep[str_detect(filetokeep$File_Path, '_R?1_')  ,]
Pair = str_replace(paired_ends$New_Location, '_1_', '_2_')
Pair = str_replace(Pair, '_R1_', '_R2_')
paired_ends$Pair = Pair
Trimmed_loc = str_replace(paired_ends$New_Location, 'Retained_SCs.*$', 'Trimmed_SCs')

paired_ends$Trimmed_loc = Trimmed_loc

Trim_1 = str_replace(paired_ends$New_Location, 'Retained_SCs', 'Trimmed_SCs')
Trim_1 = str_replace(Trim_1, '.fastq.gz', '_val_1.fq')
Trim_1 = str_replace(Trim_1, '.fastq$', '_val_1.fq')
Trim_2 = str_replace(paired_ends$Pair, 'Retained_SCs', 'Trimmed_SCs')
Trim_2 = str_replace(Trim_2, '.fastq.gz', '_val_2.fq')
Trim_2 = str_replace(Trim_2, '.fastq$', '_val_2.fq')
paired_ends$Trim_1 = Trim_1
paired_ends$Trim_2 = Trim_2
paired_ends = paired_ends %>% 
  arrange(Barcode)

paired_ends = paired_ends[, c('New_Location','Pair','Trimmed_loc','Trim_1','Trim_2','Cell_Type','Barcode')]
write.table(x = paired_ends, file = paste0(ARGs$output,format(Sys.time(), "%Y.%m.%d"),'_PairedMatches_',ARGs$name,'.tsv'), sep = '\t', quote = F, col.names = F, row.names = F)

#export cutoff values
write.table(x = as.data.frame(cutoffs), file = paste0(ARGs$output,format(Sys.time(), "%Y.%m.%d"),'_cutoffvalues_',ARGs$name,'.tsv'), sep = '\t', quote = F, col.names = T, row.names = T)
print(as.data.frame(cutoffs))

print('Succesfully completed cut-off task. Termination.')
