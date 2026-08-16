library(seqinr)

val_complete_df <- read.csv("val_complete_df.csv")

file_out <- file("validation_fasta.fasta", "w")

# Write the sequences in fasta format
for (i in 1:nrow(val_complete_df)) {
  # Write the header line
  writeLines(paste(">", val_complete_df$VirusName[i]), file_out)
  # Write the sequence line
  writeLines(val_complete_df$Sequence[i], file_out)
}

# Close the file
close(file_out)

