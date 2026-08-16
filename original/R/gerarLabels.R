library(readr)
library(Biostrings)
library(R.matlab)
library(sjmisc)

#leitura do tsv de completo
complete_variant_df1 <- read_tsv("Z:/Artigo_CGR/Dataset_TSV/DataSetGISAIDAll_1_17-01-2022.tsv")

#definindo diretório de trabalho
setwd("Z:/Artigo_CGR/Dataset_TSV/DateFastaFiles")

#abertura do fasta de exemplo
#path = 'C:/Users/AiGenomics/Desktop/AiGenomics/CGR/DataFastaFiles/'
#file = 'Start_01-Mar-2020_End_07-Mar-2020_n=1620'
#files = c("Start_01-Mar-2020_End_07-Mar-2020_n=1620","Start_02-Jan-2022_End_08-Jan-2022_n=2553","Start_10-Jan-2021_End_16-Jan-2021_n=34335","Start_26-Dec-2021_End_01-Jan-2022_n=7918")

#array com os nomes dos arquivos
file = "Start_07-Feb-2021_End_13-Feb-2021_n=38800"
files = paste0(file,'.fasta')

#Leitura do fasta
variant_fasta <- readDNAStringSet(files)

#transformar o fasta em um dataframe com nome e sequência
seq_name <- names(variant_fasta)
sequence <- paste(variant_fasta)
names_sequences <- data.frame(seq_name, sequence)
  
#comparar os nomes do tsv com fasta
concat1 <- complete_variant_df1$VirusName %in% names_sequences$seq_name

#concatenando Virus Name com Variant
data_name_labels <- data.frame(VirusName = c(complete_variant_df1$VirusName[concat1]), Variant = c(complete_variant_df1$Variant[concat1]))

#setwd('Z:\Artigo_CGR\Dataset_TSV\DateFastaFiles\unique_fasta\Experimento_1')
  
#salvar nome do virus e variante em .mat (referente a cada arquivo .fasta)
#writeMat(paste0(file, '_labels.mat'), data_name_labels = data_name_labels)
write.csv(data_name_labels, file = paste0(file, '_labels.csv'), row.names = FALSE)