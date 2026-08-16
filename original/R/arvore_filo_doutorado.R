library("seqinr")
library("stringr")
library(ggtree)
library(treeio)
library("readxl")
library("stringi")
library(ggplot2)
library("plyr")
library(dplyr)
library(Biostrings)
library(stringr)

#Set directory path
setwd('E:/FASTA/Gabriel/Results_Experimento1/Test_2/Filogenia')

#Read tree and put as.treedata 
tree <-read.nexus("fila_align_tempest.fasta.nex")
arvore <-as.treedata(tree)                                 

# File paths
fasta_files <- c("alpha_teste_Start_02_01_2022.fasta",
                 "beta_teste_Start_02_01_2022.fasta",
                 "delta_teste_Start_02_01_2022.fasta",
                 "gamma_teste_Start_02_01_2022.fasta")

# Function to extract variant from file name
get_variant <- function(file_name) {
  gsub("_teste_Start_.*", "", file_name)
}

# Function to read sequences from FASTA files
read_fasta <- function(file_path) {
  fasta <- readDNAStringSet(file_path)
  variant <- get_variant(basename(file_path))
  data.frame(Sequence = names(fasta), Variant = variant)
}

# Read sequences from all FASTA files and combine into a single data frame
sequences_df <- do.call(rbind, lapply(fasta_files, read_fasta))

#First Letter Capital and HCoV to hCoV
sequences_df$Variant <- str_to_title(sequences_df$Variant)
sequences_df$Sequence <- gsub('HCoV', 'hCoV',sequences_df$Sequence)  # Lowercase the first letter

merged_data <- left_join(as_tibble(arvore), sequences_df, by = c('label'= 'Sequence'))
merged_tree <- as.treedata(merged_data)

#Vertical tree
plot_arvore <- ggtree(merged_tree, layout = "circular", branch.length='none') + geom_tippoint(aes(color = Variant, fill=Variant), size=3) + theme_tree(legend.position='left')
plot(plot_arvore)   

ggsave('arvore_teste_1.png')

val_complete_df <- read.csv('val_complete_df.csv')
