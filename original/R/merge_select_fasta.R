setwd("E:/FASTA/Gabriel/unique_fasta/")

# Read the filelist.txt
file_list <- readLines("fasta_names.txt")

# Open a connection to merged.fasta for writing
output_file <- file("merged.fasta", "w")

# Loop through each file in the filelist and concatenate to merged.fasta
for (file_name in file_list) {
  cat(readLines(file_name), file = output_file, sep = "\n")
}

# Close the output file connection
close(output_file)