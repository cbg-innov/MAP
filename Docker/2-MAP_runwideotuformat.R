##### MAP #####

# Get arguments from bash script
args <- commandArgs(trailingOnly = TRUE)
wkdir <- args[1]

setwd(wkdir)
library(Biostrings)

# Import run-wide OTU information
df <- read.table("all_otus_cluster_info.uc", header = FALSE, sep = "\t")
names(df) <- c("Record_type",
               "Run_Wide_OTU_Cluster",
               "AlignmentLength",
               "ID",
               "Orientation",
               "X1",
               "X2",
               "X3",
               "Sample_OTU_Name",
               "Run_Wide_OTU_Centroid")

# Remove unnecessary columns and rows
df <- df[,c(1,2,9,10)]
df <- df[df$Record_type != "C",]

# Reorder rows
df <- df[order(df$Run_Wide_OTU_Cluster),]

# Add centroid name to Centroid column
df$Run_Wide_OTU_Centroid[df$Run_Wide_OTU_Centroid == "*"] <- df$Sample_OTU_Name[df$Run_Wide_OTU_Centroid == "*"]
#drop ";size="
df$Run_Wide_OTU_Centroid <- sub(";.*", "", df$Run_Wide_OTU_Centroid)   

# Import run-wide OTU consensus sequence
fasta <- readDNAStringSet("all_otus_consensus.fasta")
df.fasta <- data.frame("Run_Wide_OTU_Name" = names(fasta),
                       "Run_Wide_OTU_Consensus_Sequence" = as.character(fasta, use.names = FALSE))


# Centroid join key = bare label to  centroid=<label>;seqs=N;size=SUM.
df.fasta$Run_Wide_OTU_Centroid <- sub(";.*", "", sub("^centroid=", "", df.fasta$Run_Wide_OTU_Name))

# Component count from the seqs= field by name, not position.
df.fasta$Run_Wide_OTU_ComponentOTU_Count <- sub(".*;seqs=([0-9]+).*", "\\1", df.fasta$Run_Wide_OTU_Name)
df.master <- merge(df, df.fasta, by = "Run_Wide_OTU_Centroid", all.x = TRUE)

# Add sample OTU sequences to master table and reformat
fasta2 <- readDNAStringSet("all_otus1_corrected_otus.fasta2")
df.fasta2 <- data.frame("Sample_OTU_Name" = names(fasta2),
                        "Sample_OTU_Consensus_Sequence" = as.character(fasta2, use.names = FALSE))
df.master2 <- merge(df.master, df.fasta2, by = "Sample_OTU_Name", all.x = TRUE)
df.master2$SampleTemp <- sapply(strsplit(df.master2$Sample_OTU_Name, "\\|"), "[", 1)
df.master2$Sample <- sub("^(.*)_Rep[0-9]+$", "\\1", df.master2$SampleTemp)
df.master2$Replicate <- ifelse(
  grepl("_(Rep[0-9]+)$", df.master2$SampleTemp),
  sub("^.*_(Rep[0-9]+)$", "\\1", df.master2$SampleTemp),
  "Rep1"
)
df.master2$Replicate[!grepl("^Rep[0-9]+$", df.master2$Replicate)] <- "Rep1"
df.master2$Sample_OTU_Reads <- sapply(strsplit(df.master2$Sample_OTU_Name, ";"), "[", 2)
df.master2$Sample_OTU_Reads <- gsub("size=", "", df.master2$Sample_OTU_Reads)
df.master2 <- df.master2[,c(10,11,1,12,8,5,7,6)]
df.master2 <- df.master2[order(-as.numeric(df.master2$Run_Wide_OTU_ComponentOTU_Count), df.master2$Run_Wide_OTU_Name, -as.numeric(df.master2$Sample_OTU_Reads), df.master2$Sample_OTU_Name),]

# Output master table for later use
write.table(df.master2, "all_otus_run_wide_OTU_info.txt", quote = F, row.names = F, sep = "\t")














