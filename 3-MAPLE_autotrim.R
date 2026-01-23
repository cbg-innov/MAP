##### MAPLE v 1.0 #####

# Get argument from bash script
args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
wd <- args[2]

if (length(sampleid) == 0) {
  stop("Error: Missing input argument (sampleid).\n", call. = FALSE)
}

# Set working directory based on plate's folder, load required libraries
setwd(wd)
library("Biostrings")
library(muscle)
library(foreach)
library(doParallel)

#### Auto-trim # 1 ####

# Input ASV consensus sequences
input.fasta <- readDNAStringSet(paste0(sampleid, "_consensus.fasta"))
num_seqs <- length(input.fasta)
batch_size <- 25

# Split FASTA file into batches of 25 sequences
batches <- split(input.fasta, ceiling(seq_along(input.fasta) / batch_size))

# Function to align and trim a batch
align_and_trim <- function(batch) {
  aligned <- DNAStringSet(muscle(batch))
  
  # Get start and end positions
  column_counts <- colSums(as.matrix(aligned) == "-")
  columns_less_than_5p_gaps <- which(column_counts < 0.05 * length(aligned))
  
  first_col <- as.numeric(columns_less_than_5p_gaps[1])
  last_col <- as.numeric(columns_less_than_5p_gaps[length(columns_less_than_5p_gaps)])
  
  # Trim alignment to start and end positions, then un-align
  trimmed_alignment <- subseq(aligned, start = first_col, end = last_col)
  unaligned <- DNAStringSet(sapply(trimmed_alignment, gsub, pattern = "-", replacement = "", USE.NAMES = FALSE))
  
  return(unaligned)
}

# Use mclapply for parallel processing of batches
output_list <- mclapply(batches, align_and_trim, mc.cores = detectCores() - 1)

# Combine results
output.fasta <- Reduce(c, output_list)

#### Auto-trim # 2 ####
input.fasta <- output.fasta
num_seqs <- length(input.fasta)
batch_size <- 25

# Split FASTA file into batches of 25 sequences
batches <- split(input.fasta, ceiling(seq_along(input.fasta) / batch_size))

# Function to align and trim a batch
align_and_trim <- function(batch) {
  aligned <- DNAStringSet(muscle(batch))
  
  # Get start and end positions
  column_counts <- colSums(as.matrix(aligned) == "-")
  columns_less_than_5p_gaps <- which(column_counts < 0.05 * length(aligned))
  
  first_col <- as.numeric(columns_less_than_5p_gaps[1])
  last_col <- as.numeric(columns_less_than_5p_gaps[length(columns_less_than_5p_gaps)])
  
  # Trim alignment to start and end positions, then un-align
  trimmed_alignment <- subseq(aligned, start = first_col, end = last_col)
  unaligned <- DNAStringSet(sapply(trimmed_alignment, gsub, pattern = "-", replacement = "", USE.NAMES = FALSE))
  
  return(unaligned)
}

# Use mclapply for parallel processing of batches
output_list2 <- mclapply(batches, align_and_trim, mc.cores = detectCores() - 1)

# Combine results
output.fasta2 <- Reduce(c, output_list)

# Output trimmed sequences
writeXStringSet(output.fasta2,sprintf("%s_finalasvs.fas",sampleid), format = "fasta")
