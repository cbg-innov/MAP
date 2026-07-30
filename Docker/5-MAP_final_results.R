##### MAP #####

# Get argument from bash script
args <- commandArgs(trailingOnly = TRUE)
runid <- args[1]
marker <- args[2]
min.reads <- args[3]
min.reads <- as.numeric(min.reads)
num.reps <- args[4]
num.reps <- as.numeric(num.reps)
amp.size <- args[5]
amp.size <- as.numeric(amp.size)
wkdir <- args[6]
reflib <- args[7]
componentreads <- args[8]
param_file <- args[9]

setwd(wkdir)
library(Biostrings)
library(plyr)
library(dplyr)
library(openxlsx)
library(tidyr)
library(readxl)
library(stringr)

# Import master tables
df <- read.table("df.txt", sep = "\t", header = T)

#####################################################################################
######################## GENERATE REPLICATE TABLE ###################################
#####################################################################################

# Import metadata from parameters file


# Unique original sample names (may contain "#")
mapping <- read.delim(paste("temp_mapping_", runid, "_", marker, "_", amp.size, ".txt", sep = ""), header = F, check.names = F, fill = T)
sample.list <- mapping %>%
  # Remove trailing _Rep<number>
  mutate(sample.list = str_remove(V3, "_Rep\\d+$")) %>%

  # Keep unique
  distinct(sample.list) %>%
  pull(sample.list)

# Build mapping: sanitized name -> original name
orig <- sample.list
sanitized <- gsub("#", "-", orig)
sample_map <- setNames(orig, sanitized)

# Restore original sample names in df$Sample
df$Sample <- ifelse(df$Sample %in% names(sample_map),
                    sample_map[df$Sample],
                    df$Sample)

# Remove any non-sample rows (i.e. positive and negative controls)
df <- df[df$Sample %in% sample.list,]

# Rename column headers
names(df) <- c("Sample", "Replicate", "Sample_OTU_Name", "Reads", "Sample_OTU_Sequence", 
               "Run_OTU_Name", "Run_OTU_Sequence", "Tax")

# Add sequence length and ambiguous base count to tables
df$Sample_OTU_Length <- nchar(df$Sample_OTU_Sequence)
df$Sample_OTU_Ns <- nchar(df$Sample_OTU_Sequence) - nchar(gsub("N", "", df$Sample_OTU_Sequence))
df$Sample_OTU_Ns <- as.numeric(df$Sample_OTU_Ns)

# Rename run-wide OTU names to more human-readable form
runwide.names <- df$Run_OTU_Name[order(-df$Reads)]
runwide.names <- runwide.names[!duplicated(runwide.names)]
num_digits <- nchar(length(runwide.names))
temp.df <- data.frame(
  "OldName" = runwide.names,
  "NewName" = sprintf(paste0("OTU_%0", num_digits, "d"), 1:length(runwide.names))
)
df <- merge(df, temp.df, by.x = "Run_OTU_Name", by.y = "OldName", all.x = TRUE)
df <- df[,c("Sample", "Replicate", "NewName", "Reads", "Sample_OTU_Name", "Sample_OTU_Length", "Sample_OTU_Ns", "Sample_OTU_Sequence",
            "Run_OTU_Sequence", "Tax")]
names(df)[c(3,9)] <- c("OTU_ID", "Run_OTU_Consensus_Sequence")

# Expand taxonomy into individual columns
nmax <- max(lengths(strsplit(df$Tax, ",")))
rank_names <- c("Kingdom", "Phylum", "Class", "Order",
                "Family", "Genus", "Species")
df <- df %>%
  separate(
    col   = Tax,
    into  = rank_names,
    sep   = ",",
    fill  = "right"
  )

# Remove sintax rank indicators from data
df$Kingdom <- gsub("k:", "", df$Kingdom)
df$Phylum <- gsub("p:", "", df$Phylum)
df$Class <- gsub("c:", "", df$Class)
df$Order <- gsub("o:", "", df$Order)
df$Family <- gsub("f:", "", df$Family)
df$Genus <- gsub("g:", "", df$Genus)
df$Species <- gsub("s:", "", df$Species)

# Backup table for later
df.bkp <- df

# Sort and format final table
df <- df[order(df$Sample, df$Replicate, -df$Reads),]
df[is.na(df)] <- ""

#####################################################################################
########################### GENERATE SAMPLE TABLE ###################################
#####################################################################################
df <- df[order(df$Sample, df$OTU_ID, -df$Reads),]

# Summarize the grouped data
df.summary <- df %>%
  group_by(Sample, OTU_ID) %>%
  summarize(
    Reads = sum(Reads),
    OTU_Consensus_Sequence = first(Run_OTU_Consensus_Sequence),
    Seq_Length = nchar(first(Run_OTU_Consensus_Sequence)),
    Number_N = nchar(first(Run_OTU_Consensus_Sequence)) - nchar(gsub("N", "", first(Run_OTU_Consensus_Sequence))),
    Replicates = n_distinct(Replicate),
    Kingdom = first(Kingdom),
    Phylum = first (Phylum),
    Class = first (Class),
    Order = first (Order),
    Family = first (Family),
    Genus = first (Genus),
    Species = first (Species),
    .groups = "drop"
  )

# Replace all NA values with blanks
df.summary[is.na(df.summary)] <- ""

# Sort by sample and read count
df.summary <- df.summary[order(df.summary$Sample, -df.summary$Reads),]

#####################################################################################
########################### IMPORT METADATA TABLE ###################################
#####################################################################################
# Add sample metadata to results file
metadata <- read.table(sprintf("../metadata_%s.txt", runid), header = T, sep = "\t", comment.char = "", fill = T, check.names = F)
names(metadata) <- c("Sample", "Collection Site", "Latitude", "Longitude", "Collection Start Date", "Collection End Date")

#####################################################################################
################## MERGE OTU COMPONENT READS BY RUN-WIDE OTU ########################
#####################################################################################
if(componentreads == 1){
  # Strip ";size=" from OTU name column to match component read file names
  df$Sample_OTU_Name <- sub(";size=.*$", "", df$Sample_OTU_Name)
  
  # Create import/export directory
  fasta_dir <- sprintf("./%s_%s_%s_OTU_Component_Reads", runid, marker, amp.size)
  if (!dir.exists(fasta_dir)) stop("FASTA directory does not exist!")
  
  # Group by Sample and OTU_ID
  df %>%
    group_by(Sample, OTU_ID) %>%
    group_walk(~ {
      sample <- .y$Sample
      otu <- .y$OTU_ID
      sum_reads <- sum(.x$Reads)
      out_file <- file.path(fasta_dir, sprintf("%s|%s|reads-%d_ComponentReads.fasta", sample, otu, sum_reads))
      
      # Construct input FASTA filenames
      fasta_files <- file.path(fasta_dir, paste0(.x$Sample_OTU_Name, "_ComponentReads.fasta"))
      fasta_files <- fasta_files[file.exists(fasta_files)]
      if (length(fasta_files) == 0) {
        warning(sprintf("No FASTA files found for %s|%s", sample, otu))
        return()
      }
      
      # Read all lines, flatten list into a single character vector
      all_lines <- unlist(lapply(fasta_files, readLines), use.names = FALSE)
      
      # Write merged FASTA
      writeLines(all_lines, out_file)
      
      # Delete original files
      file.remove(fasta_files)
      
      message(sprintf("Merged %d files → %s", length(fasta_files), out_file))
    })
}

# Remove Sample_OTU_Name column from df
df <- df[, -5]

#####################################################################################
########################### EXPORT AS EXCEL FILE ####################################
#####################################################################################
# Create workbook
wb <- createWorkbook()

# Add "By Sample" sheet
addWorksheet(wb, "By Sample", gridLines = TRUE)
writeData(wb, "By Sample", df.summary)
header.style <- createStyle(fontSize = 13, textDecoration = "bold", border = "TopBottomLeftRight")
addStyle(wb, "By Sample", cols = 1:ncol(df.summary), rows = 1, style = header.style, gridExpand = TRUE)
centre.style <- createStyle(halign = "center")
addStyle(wb, "By Sample", cols = c(4,5,8:11), rows = 1:nrow(df.summary)+1, style = centre.style, gridExpand = TRUE)
setColWidths(wb,"By Sample",cols = 1,widths = "13") #sample
setColWidths(wb,"By Sample",cols = 2,widths = "13") #otu 
setColWidths(wb,"By Sample",cols = 3,widths = "13") #reads
setColWidths(wb,"By Sample",cols = 4,widths = "25") #seq
setColWidths(wb,"By Sample",cols = 5,widths = "22") #length
setColWidths(wb,"By Sample",cols = 6,widths = "18") #N
setColWidths(wb,"By Sample",cols = 7,widths = "13") #replicate
setColWidths(wb,"By Sample",cols = 8,widths = "22") #k
setColWidths(wb,"By Sample",cols = 9,widths = "22") #p
setColWidths(wb,"By Sample",cols = 10,widths = "20") #c
setColWidths(wb,"By Sample",cols = 11,widths = "20") #o
setColWidths(wb,"By Sample",cols = 12,widths = "22") #f
setColWidths(wb,"By Sample",cols = 13,widths = "22") #g
setColWidths(wb,"By Sample",cols = 14,widths = "25") #s

# Add "By Replicate" sheet
addWorksheet(wb, "By Replicate", gridLines = TRUE)
writeData(wb, "By Replicate", df)
header.style <- createStyle(fontSize = 13, textDecoration = "bold")
addStyle(wb, "By Replicate", cols = 1:ncol(df), rows = 1, style = header.style, gridExpand = TRUE)
setColWidths(wb,"By Replicate",cols = 1,widths = "13") #sample
setColWidths(wb,"By Replicate",cols = 2,widths = "13") #replicate
setColWidths(wb,"By Replicate",cols = 3,widths = "13") #avs 
setColWidths(wb,"By Replicate",cols = 4,widths = "13") #reads
setColWidths(wb,"By Replicate",cols = 5,widths = "22") #length
setColWidths(wb,"By Replicate",cols = 6,widths = "18") #N
setColWidths(wb,"By Replicate",cols = 7,widths = "25") #seq
setColWidths(wb,"By Replicate",cols = 8,widths = "35") #runseq
setColWidths(wb,"By Replicate",cols = 9,widths = "22") #k
setColWidths(wb,"By Replicate",cols = 10,widths = "22") #p
setColWidths(wb,"By Replicate",cols = 11,widths = "20") #c
setColWidths(wb,"By Replicate",cols = 12,widths = "20") #o
setColWidths(wb,"By Replicate",cols = 13,widths = "22") #f
setColWidths(wb,"By Replicate",cols = 14,widths = "22") #g
setColWidths(wb,"By Replicate",cols = 15,widths = "25") #s

# Add "Sample Metadata" sheet
addWorksheet(wb, "Sample Metadata", gridLines = TRUE)
writeData(wb, "Sample Metadata", metadata)
header.style <- createStyle(fontSize = 13, textDecoration = "bold", border = "TopBottomLeftRight")
addStyle(wb, "Sample Metadata", cols = 1:ncol(metadata), rows = 1, style = header.style, gridExpand = TRUE)
setColWidths(wb,"Sample Metadata",cols = 1,widths = "13") #sample
setColWidths(wb,"Sample Metadata",cols = 2,widths = "20") #site
setColWidths(wb,"Sample Metadata",cols = 3,widths = "13") #lat
setColWidths(wb,"Sample Metadata",cols = 4,widths = "13") #long
setColWidths(wb,"Sample Metadata",cols = 5,widths = "21") #startdate
setColWidths(wb,"Sample Metadata",cols = 6,widths = "21") #enddate

# Save workbook
saveWorkbook(wb, sprintf("Metabarcoding Results - %s_%s_%s.xlsx", runid, marker, amp.size), overwrite = TRUE)

#####################################################################################
########################### PREPARE FOR BIN MATCH (COI-5P ONLY) #####################
#####################################################################################
if(marker == "COI-5P" & amp.size >= 300){
  output.fasta <- DNAStringSet(df.summary$OTU_Consensus_Sequence)
  names(output.fasta) <- df.summary$OTU_ID
  output.fasta <- output.fasta[!duplicated(output.fasta)]
  writeXStringSet(output.fasta, "temp_toBIN.fasta")
}

