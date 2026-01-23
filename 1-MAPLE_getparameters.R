##### MAPLE v 1.0 #####

# Get argument from bash script
args <- commandArgs(trailingOnly = TRUE)
wkdir <- args[1]

# Set working directory and load required libraries
setwd(wkdir)
library("readxl")

# Import then export run info from parameters file
runinfo <- data.frame(read_excel(list.files(pattern = "*.xlsx"), sheet = "Run Parameters", range = "B1:B6", col_names = F))
runinfo$...1 <- ifelse(
  grepl("^[0-9.]+$", runinfo$...1),
  format(round(as.numeric(runinfo$...1), 1), nsmall = 1),
  runinfo$...1
)
write.table(runinfo, "runinfo.txt", append = F, quote = F, row.names = F, sep = "\t", col.names = F)

# Import then export UMI map and related metadata
umimap <- data.frame(read_excel(list.files(pattern = "*.xlsx"), sheet = "Run Parameters", skip = 7))
write.table(umimap, sprintf("mapping_%s.txt", runinfo$...1[1]), append = F, quote = F, row.names = F, sep = "\t")

# Import then export sample metadata
metadata <- data.frame(read_excel(list.files(pattern = "*.xlsx"), sheet = "Sample Metadata"))
write.table(metadata, sprintf("metadata_%s.txt", runinfo$...1[1]), append = F, quote = F, row.names = F, sep = "\t")















