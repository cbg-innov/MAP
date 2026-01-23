##### MAPLE v 1.0 #####

# Get arguments from bash script
args <- commandArgs(trailingOnly = TRUE)
wkdir <- args[1]

setwd(wkdir)
library(Biostrings)

# Import run-wide ASV information
df <- read.table("all_asvs_cluster_info.uc", header = FALSE, sep = "\t")
names(df) <- c("Record_type",
               "Run_Wide_ASV_Cluster",
               "AlignmentLength",
               "ID",
               "Orientation",
               "X1",
               "X2",
               "X3",
               "Sample_ASV_Name",
               "Run_Wide_ASV_Centroid")

# Remove unnecessary columns and rows
df <- df[,c(1,2,9,10)]
df <- df[df$Record_type != "C",]

# Reorder rows
df <- df[order(df$Run_Wide_ASV_Cluster),]

# Add centroid name to Centroid column
df$Run_Wide_ASV_Centroid[df$Run_Wide_ASV_Centroid == "*"] <- df$Sample_ASV_Name[df$Run_Wide_ASV_Centroid == "*"]

# Import run-wide ASV consensus sequence
fasta <- readDNAStringSet("all_asvs_consensus.fasta")
df.fasta <- data.frame("Run_Wide_ASV_Name" = names(fasta),
                       "Run_Wide_ASV_Consensus_Sequence" = as.character(fasta, use.names = FALSE))

# Reformat run-wide ASV sequence headers and create master table
part1 <- sapply(strsplit(df.fasta$Run_Wide_ASV_Name, ";"), "[", 1)
part2 <- sapply(strsplit(df.fasta$Run_Wide_ASV_Name, ";"), "[", 2)
part3 <- sapply(strsplit(df.fasta$Run_Wide_ASV_Name, ";"), "[", 3)

df.fasta$Run_Wide_ASV_Centroid <- paste0(part1, ";", part2)
df.fasta$Run_Wide_ASV_Centroid <- gsub("centroid=", "", df.fasta$Run_Wide_ASV_Centroid)
df.fasta$Run_Wide_ASV_ComponentASV_Count <- gsub("seqs=", "", part3)

df.master <- merge(df, df.fasta, by = "Run_Wide_ASV_Centroid", all.x = TRUE)

# Add sample ASV sequences to master table and reformat
fasta2 <- readDNAStringSet("all_asvs.fasta2")
df.fasta2 <- data.frame("Sample_ASV_Name" = names(fasta2),
                        "Sample_ASV_Consensus_Sequence" = as.character(fasta2, use.names = FALSE))
df.master2 <- merge(df.master, df.fasta2, by = "Sample_ASV_Name", all.x = TRUE)
df.master2$SampleTemp <- sapply(strsplit(df.master2$Sample_ASV_Name, "\\|"), "[", 1)
df.master2$Sample <- sub("^(.*)_Rep[0-9]+$", "\\1", df.master2$SampleTemp)
df.master2$Replicate <- ifelse(
  grepl("_(Rep[0-9]+)$", df.master2$SampleTemp),
  sub("^.*_(Rep[0-9]+)$", "\\1", df.master2$SampleTemp),
  "Rep1"
)
df.master2$Replicate[!grepl("^Rep[0-9]+$", df.master2$Replicate)] <- "Rep1"
df.master2$Sample_ASV_Reads <- sapply(strsplit(df.master2$Sample_ASV_Name, ";"), "[", 2)
df.master2$Sample_ASV_Reads <- gsub("size=", "", df.master2$Sample_ASV_Reads)
df.master2 <- df.master2[,c(10,11,1,12,8,5,7,6)]
df.master2 <- df.master2[order(-as.numeric(df.master2$Run_Wide_ASV_ComponentASV_Count), df.master2$Run_Wide_ASV_Name, -as.numeric(df.master2$Sample_ASV_Reads), df.master2$Sample_ASV_Name),]

# Output master table for later use
write.table(df.master2, "all_asvs_run_wide_ASV_info.txt", quote = F, row.names = F, sep = "\t")














