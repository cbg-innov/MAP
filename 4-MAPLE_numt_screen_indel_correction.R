##### MAPLE v 1.0 #####

# Get argument from bash script
args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]
cores <- args[2]
cores <- as.numeric(cores) - 1
ampsize <- args[3]
ampsize <- as.numeric(ampsize)
wkdir <- args[4]

if (length(sampleid) == 0) {
  stop("Error: Missing input argument (sampleid).\n", call. = FALSE)
}
if (length(cores) == 0) {
  stop("Error: Missing input argument (cores).\n", call. = FALSE)
}
if (length(ampsize) == 0) {
  stop("Error: Missing input argument (ampsize).\n", call. = FALSE)
}

# Set working directory
setwd(wkdir)
library("Biostrings")
library("dplyr")
library("coil")
library("parallel")
library(DECIPHER)
library(pwalign)

# Import ASVs and convert to table
asv.seqs <- readDNAStringSet(sprintf("%s_FinalASVs.fast", sampleid))
df <- data.frame("SeqName" = names(asv.seqs),
                 "Sequence" = as.character(asv.seqs, use.names = FALSE))
# Run HMM screen on ASVs
df$trans_table = lapply(1:nrow(df), function(i){
  which_trans_table("Arthropoda")
})

df$trans_table <- as.numeric(df$trans_table)
options(mc.cores = cores)
temp.list = mclapply(1:nrow(df), function(i){
  coi5p_pipe(df$Sequence[i],
             name = df$SeqName[i],
             triple_translate = TRUE,
             indel_threshold = -425.21,
             trans_table = df$trans_table[i])
})

df$INDEL <- unlist(lapply(temp.list, function(x){ x$indel_likely }))
df$STOPCODON <- unlist(lapply(temp.list, function(x){ x$stop_codons }))
df$HMM <- ifelse(df$INDEL == FALSE & df$STOPCODON == FALSE, "OK", "HMM_ISSUE")

# Add sequence length and ambiguous base count to master no-hit table
df$SeqLength <- nchar(df$Sequence)
df$Ns <- nchar(df$Sequence) - nchar(gsub("N", "", df$Sequence))

# Identify ASVs that are obvious NUMTs
df$NUMT <- ifelse(df$INDEL == FALSE & df$STOP == TRUE, TRUE,
                               ifelse(df$SeqLength < (ampsize - 12), TRUE,
                                      ifelse(df$SeqLength > (ampsize + 12), TRUE,
                                             FALSE)))

# Remove unknown ASVs that are NUMTs or have more than 6 Ns
df <- df[!df$NUMT,]
df$Ns <- as.numeric(df$Ns)
df <- df[df$Ns <= 6,]

# Generate new FASTA file
output.fasta <- DNAStringSet(df$Sequence)
output.seq.names <- df$SeqName
names(output.fasta) <- output.seq.names

#################################################################################################
###################### Try to correct polymerase error indels ###################################
#################################################################################################

seqs <- output.fasta
expected_lengths <- c(646, 649, 652, 655, 658, 661, 664, 667, 670)

# Partition sequences
seq_lengths <- width(seqs)
good_seqs <- seqs[seq_lengths %in% expected_lengths]
bad_seqs  <- seqs[!seq_lengths %in% expected_lengths]

n_good <- length(good_seqs)
n_bad  <- length(bad_seqs)
seqs <- c(good_seqs, bad_seqs)
dists <- DistanceMatrix(seqs, method = "overlap", includeTerminalGaps = TRUE, processors = NULL)

# Find nearest neighbor good sequence for each bad sequence
best_refs <- sapply(seq_along(bad_seqs), function(i) {
  bad_index <- n_good + i
  distances_to_good <- dists[bad_index, 1:n_good]
  best_match_index <- which.min(distances_to_good)
  names(good_seqs)[best_match_index]
})

# Align bad sequences with nearest neighbor
aligned <- lapply(seq_along(bad_seqs), function(i) {
  pwalign::pairwiseAlignment(
    pattern = bad_seqs[[i]],
    subject = good_seqs[[best_refs[i]]],
    type = "global",
    gapOpening = 15,
    gapExtension = 2
  )
})

# Function to correct single bp indels in one alignment
correct_single_bp_indels <- function(aln) {
  bad_aln <- as.character(pwalign::pattern(aln))
  ref_aln <- as.character(subject(aln))
  
  stopifnot(nchar(bad_aln) == nchar(ref_aln))
  
  # Count single-bp indels
  in_gap <- FALSE
  gap_length <- 0
  n_indels <- 0
  for (i in seq_len(nchar(bad_aln))) {
    b <- substr(bad_aln, i, i)
    r <- substr(ref_aln, i, i)
    if (b == "-" || r == "-") {
      if (!in_gap) { in_gap <- TRUE; gap_length <- 1 } else { gap_length <- gap_length + 1 }
    } else {
      if (in_gap && gap_length == 1) n_indels <- n_indels + 1
      in_gap <- FALSE
      gap_length <- 0
    }
  }
  if (in_gap && gap_length == 1) n_indels <- n_indels + 1
  
  # If no single-bp indels, return cleaned sequence immediately
  if (n_indels == 0) return(DNAString(gsub("-", "", bad_aln)))
  
  # Otherwise, correct single-bp indels
  bad_chars <- strsplit(bad_aln, "")[[1]]
  ref_chars <- strsplit(ref_aln, "")[[1]]
  corrected <- character(length(bad_chars))
  
  j <- 1
  while (j <= length(bad_chars)) {
    b <- bad_chars[j]
    r <- ref_chars[j]
    
    if (b == "-" || r == "-") {
      # Start of gap block
      gap_start <- j
      while (j <= length(bad_chars) && (bad_chars[j] == "-" || ref_chars[j] == "-")) j <- j + 1
      gap_end <- j - 1
      gap_len <- gap_end - gap_start + 1
      
      if (gap_len == 1) {
        if (b == "-" && r != "-") corrected[gap_start] <- "N"  # single deletion
        # single insertion: skip base, leave corrected as ""
      } else {
        corrected[gap_start:gap_end] <- bad_chars[gap_start:gap_end]
      }
    } else {
      corrected[j] <- b
      j <- j + 1
    }
  }
  
  DNAString(paste0(corrected[corrected != ""], collapse = ""))
}

# Apply correction to all bad sequences
corrected_seqs <- lapply(aligned, correct_single_bp_indels)
names(corrected_seqs) <- names(bad_seqs)

# Final combined fasta: good + corrected bad sequences
output.fasta <- c(good_seqs, DNAStringSet(unlist(corrected_seqs)))
out <- DNAStringSet(gsub("-", "", output.fasta))

# Output trimmed and corrected/uncorrected sequences
writeXStringSet(out, sprintf("%s_FinalASVs.fasta", sampleid))









