#get argument from bash script
args <- commandArgs(trailingOnly = TRUE)
runid <- args[1]
wkdir <- args[2]
marker <- args[3]
amp.size <- args[4]

# set working directory
setwd(wkdir)
library(readxl)
library(openxlsx)

# import Excel tabs
file <- paste("Metabarcoding Results - ", runid, "_", marker, "_", amp.size, ".xlsx", sep = "")
wb1 <- read_excel(file, sheet = "By Sample")
wb2 <- read_excel(file, sheet = "By Replicate")
wb3 <- read_excel(file, sheet = "Sample Metadata")

# import BIN match results
df <- read.delim("final_output_with_bin_match.txt", check.names = F, header = FALSE)
df <- df[,c(1:3,11,12,4:10,13)]
names(df) <- c("Query",
               "BIN BOLD ProcessID",
               "BIN Hit",
               "%ID match to BIN",
               "Overlap (bp) with BIN",
               "BIN Kingdom",
               "BIN Phylum",
               "BIN Class",
               "BIN Order",
               "BIN Family",
               "BIN Genus",
               "BIN Species",
               "BIN Match Status")

# add BIN match results to Excel output
output <- merge(wb1, df, by.x = "OTU_ID", by.y = "Query", all.x = TRUE)
# Order by NAME, not position, to account for extra confidence columns
bin_cols <- c("BIN Match Status", "BIN Hit", "BIN BOLD ProcessID", "%ID match to BIN",
              "Overlap (bp) with BIN", "BIN Kingdom", "BIN Phylum", "BIN Class", "BIN Order",
              "BIN Family", "BIN Genus", "BIN Species")
# Column order: sample/OTU info, then BIN results, then the SINTAX taxonomy and its confidences.
lead_cols  <- c("Sample", "Reads", "OTU_ID")
tax_cols   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
conf_cols  <- grep("_Confidence$", names(wb1), value = TRUE)
other_cols <- setdiff(names(wb1), c(lead_cols, tax_cols, conf_cols))
output <- output[, c(lead_cols, other_cols, bin_cols, tax_cols, conf_cols)]
output <- output[order(output$Sample, -output$Reads),]

wb <- createWorkbook()

addWorksheet(wb, "By Sample", gridLines = TRUE)
writeData(wb, "By Sample", output)
header.style <- createStyle(fontSize = 13, textDecoration = "bold", border = "TopBottomLeftRight")
addStyle(wb, "By Sample", cols = 1:ncol(output), rows = 1, style = header.style, gridExpand = TRUE)
centre.style <- createStyle(halign = "center")
conf.style <- createStyle(halign = "center", numFmt = "0.00")
col_of <- function(x, nm) match(nm, names(x))
addStyle(wb, "By Sample", cols = col_of(output, c("Reads", "Seq_Length", "Number_N", "Replicates", "BIN Match Status",
                                                   "%ID match to BIN", "Overlap (bp) with BIN")),
         rows = 1:nrow(output)+1, style = centre.style, gridExpand = TRUE)
conf_cols <- grep("_Confidence$", names(output), value = TRUE)
if (length(conf_cols)) addStyle(wb, "By Sample", cols = col_of(output, conf_cols), rows = 1:nrow(output)+1, style = conf.style, gridExpand = TRUE)
sample_widths <- c("Sample" = 13, "Reads" = 13, "OTU_ID" = 13, "OTU_Consensus_Sequence" = 25, "Seq_Length" = 22,
                   "Number_N" = 18, "Replicates" = 13, "Kingdom" = 22, "Phylum" = 22, "Class" = 20, "Order" = 20,
                   "Family" = 22, "Genus" = 22, "Species" = 25, "BIN Match Status" = 15, "BIN Hit" = 20,
                   "BIN BOLD ProcessID" = 20, "%ID match to BIN" = 13, "Overlap (bp) with BIN" = 13,
                   "BIN Kingdom" = 22, "BIN Phylum" = 22, "BIN Class" = 20, "BIN Order" = 20, "BIN Family" = 22,
                   "BIN Genus" = 22, "BIN Species" = 25)
sample_widths[conf_cols] <- 22
for (nm in names(sample_widths)) {
  i <- col_of(output, nm)
  if (!is.na(i)) setColWidths(wb, "By Sample", cols = i, widths = sample_widths[[nm]])
}

addWorksheet(wb, "By Replicate", gridLines = TRUE)
writeData(wb, "By Replicate", wb2)
header.style <- createStyle(fontSize = 13, textDecoration = "bold")
addStyle(wb, "By Replicate", cols = 1:ncol(wb2), rows = 1, style = header.style, gridExpand = TRUE)
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
rep_conf <- grep("_Confidence$", names(wb2), value = TRUE)
if (length(rep_conf)) {
  addStyle(wb, "By Replicate", cols = col_of(wb2, rep_conf), rows = 1:nrow(wb2)+1, style = conf.style, gridExpand = TRUE)
  setColWidths(wb, "By Replicate", cols = col_of(wb2, rep_conf), widths = 22)
}

addWorksheet(wb, "Sample Metadata", gridLines = TRUE)
writeData(wb, "Sample Metadata", wb3)
header.style <- createStyle(fontSize = 13, textDecoration = "bold", border = "TopBottomLeftRight")
addStyle(wb, "Sample Metadata", cols = 1:ncol(wb3), rows = 1, style = header.style, gridExpand = TRUE)
setColWidths(wb,"Sample Metadata",cols = 1,widths = "13") #sample
setColWidths(wb,"Sample Metadata",cols = 2,widths = "20") #site
setColWidths(wb,"Sample Metadata",cols = 3,widths = "13") #lat
setColWidths(wb,"Sample Metadata",cols = 4,widths = "13") #long
setColWidths(wb,"Sample Metadata",cols = 5,widths = "21") #startdate
setColWidths(wb,"Sample Metadata",cols = 6,widths = "21") #enddate

saveWorkbook(wb, sprintf("Metabarcoding Results - %s_%s_%s.xlsx", runid, marker, amp.size), overwrite = TRUE)









