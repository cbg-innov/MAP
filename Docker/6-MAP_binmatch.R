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
output <- output[,c(2,3,1,4:14,26,16,15,17:25)]
output <- output[order(output$Sample, -output$Reads),]

wb <- createWorkbook()

addWorksheet(wb, "By Sample", gridLines = TRUE)
writeData(wb, "By Sample", output)
header.style <- createStyle(fontSize = 13, textDecoration = "bold", border = "TopBottomLeftRight")
addStyle(wb, "By Sample", cols = 1:ncol(output), rows = 1, style = header.style, gridExpand = TRUE)
centre.style <- createStyle(halign = "center")
addStyle(wb, "By Sample", cols = c(2,5,6,7,15,18,19), rows = 1:nrow(output)+1, style = centre.style, gridExpand = TRUE)
setColWidths(wb,"By Sample",cols = 1,widths = "13") #sample
setColWidths(wb,"By Sample",cols = 2,widths = "13") #reads 
setColWidths(wb,"By Sample",cols = 3,widths = "13") #otu
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
setColWidths(wb,"By Sample",cols = 15,widths = "15") #bin status
setColWidths(wb,"By Sample",cols = 16,widths = "20") #bin hit
setColWidths(wb,"By Sample",cols = 17,widths = "20") #bin PID
setColWidths(wb,"By Sample",cols = 18,widths = "13") #bin %id
setColWidths(wb,"By Sample",cols = 19,widths = "13") #bin overlap
setColWidths(wb,"By Sample",cols = 20,widths = "22") #bin k
setColWidths(wb,"By Sample",cols = 21,widths = "22") #bin p
setColWidths(wb,"By Sample",cols = 22,widths = "20") #bin c
setColWidths(wb,"By Sample",cols = 23,widths = "20") #bin o
setColWidths(wb,"By Sample",cols = 24,widths = "22") #bin f
setColWidths(wb,"By Sample",cols = 25,widths = "22") #bin g
setColWidths(wb,"By Sample",cols = 26,widths = "25") #bin s

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









