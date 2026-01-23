#get argument from bash script (run name)
args <- commandArgs(trailingOnly = TRUE)
min.reads <- args[1]
min.reads <- as.numeric(min.reads)
num.reps <- args[2]
num.reps <- as.numeric(num.reps)
wkdir <- args[3]
runid <- args[4]
marker <- args[5]
ampsize <- args[6]

# set working directory
setwd(wkdir)
library("Biostrings")
library("plyr")
library("dplyr")
library("coil")
library("parallel")
library(openxlsx)
library(tidyr)


# import sintax results
sintax.input <- read.delim(
  "temp_sintax_output.txt",
  header       = FALSE,
  fill         = TRUE,
  comment.char = "",
  col.names    = paste0("V", 1:4),
  colClasses   = "character"
)
names(sintax.input) <- c("Query","Tax","strand","Tax_final")

# import master ASV table
df <- read.table("all_asvs_run_wide_ASV_info.txt", header = TRUE, sep = "\t")

# add sintax results to master ASV table
df2 <- merge(df, sintax.input[,c(1,4)], by.x = "Run_Wide_ASV_Name", by.y = "Query", all.x = TRUE)

# reformat master ASV table
df2 <- df2[,c(2:6,1,7:9)]
df2 <- df2[order(df2$Sample, df2$Replicate, -df2$Sample_ASV_Reads),]

###############################################################
############### DE-NOISE HITS #################################
###############################################################

# merge all identical hits within sample-rep into a single hit with reads summed
df2$SampleRep <- paste(df2$Sample, df2$Replicate, sep = "-")

df2 <- data.frame(df2 %>%
                    dplyr::group_by(SampleRep, Run_Wide_ASV_Name) %>%
                    dplyr::summarize(
                      Sample = dplyr::first(Sample), 
                      Replicate = dplyr::first(Replicate), 
                      Sample_ASV_Name = dplyr::first(Sample_ASV_Name),
                      Sample_ASV_Reads = sum(Sample_ASV_Reads),
                      Sample_ASV_Consensus_Sequence = dplyr::first(Sample_ASV_Consensus_Sequence),
                      Run_Wide_ASV_Name = dplyr::first(Run_Wide_ASV_Name),
                      Run_Wide_ASV_ComponentASV_Count = dplyr::first(Run_Wide_ASV_ComponentASV_Count),
                      Run_Wide_ASV_Consensus_Sequence = dplyr::first(Run_Wide_ASV_Consensus_Sequence),
                      Tax = dplyr::first(Tax_final),
                      .groups = 'drop'
                    ))

df2 <- df2[,c(3:7,2,8:10)]

#if negative controls used, calculate noise threshold and proportionally subtract reads in neg control wells
neg.rows <- grep("NEG|NEGATIVE|negative|CONTROL|BLANK|blank",df2$Sample)
if (length(neg.rows) != 0){
  #extract ASVs in negative control wells and output for user's files
  write.table(df2[neg.rows,], sprintf("%s_%s_%sbp_NegativeControlASVs.tsv", runid, marker, ampsize), quote = F, row.names = F, sep = "\t")
  
  #calculate mean reads per negative control well
  hits.in.neg <- df2[neg.rows,c("Sample_ASV_Reads","Run_Wide_ASV_Name")]
  hits.in.neg.sum <- aggregate(Sample_ASV_Reads~Run_Wide_ASV_Name, hits.in.neg, sum)
  hits.in.neg.mean <- hits.in.neg.sum
  hits.in.neg.mean$Mean.Reads.Per.Well <- hits.in.neg.mean$Sample_ASV_Reads/length(neg.rows)
  
  #calculate average noise threshold based on proportion of reads that are in negs and also samples
  sample.rows <- df2[-neg.rows,]
  Neg.hits.in.samples <- sample.rows[grep(paste(gsub("\\|", "\\\\|", hits.in.neg.mean$Run_Wide_ASV_Name),collapse = "|"), sample.rows$Run_Wide_ASV_Name),c(4,6)]
  combined <- merge(Neg.hits.in.samples, hits.in.neg.mean, by ="Run_Wide_ASV_Name")
  combined <- combined [,c(1,2,4)]
  names(combined) <- c("ASV","Reads","Mean.Reads.Per.NEG.Well")
  combined$Proportion <- combined$Mean.Reads.Per.NEG.Well/combined$Reads
  combined.filtered <- combined[combined$Proportion <= 0.0033,] #<<< this value was set based on observations of the data
  
  #if background noise was too low to detect, manually set it based on historic observations
  if (nrow(combined.filtered) == 0) {
    average.noise.threshold <- 0.0001
  }else{
    combined.filtered.mean <- aggregate(Proportion~ASV, combined.filtered, mean)
    average.noise.threshold <- mean(combined.filtered.mean$Proportion)
  }
  
  #proportionally subtract reads in control wells
  hits.in.neg.mean.2 <- aggregate(Sample_ASV_Reads~Run_Wide_ASV_Name, hits.in.neg, mean)
  hits.in.neg.mean.2$Reads <- round(hits.in.neg.mean.2$Sample_ASV_Reads, digits = 0)
  names(hits.in.neg.mean.2)[3] <- "Mean Reads"
  hits.in.neg.mean.2 <- hits.in.neg.mean.2[,-2]
  df2$Reads.to.subtract <- hits.in.neg.mean.2$`Mean Reads`[match(df2$Run_Wide_ASV_Name,hits.in.neg.mean.2$Run_Wide_ASV_Name)]
  df2$Reads.to.subtract[is.na(df2$Reads.to.subtract)] <- 0
  df2$Adjusted.Reads <- df2$Sample_ASV_Reads - df2$Reads.to.subtract
  df2 <- df2[df2$Adjusted.Reads > 0,]
} else { #this is done if there are no negative control wells
  
  #extract ASVs in negative control wells and output for user's files
  write.table(df2[0,], sprintf("%s_%s_%sbp_NegativeControlASVs.tsv", runid, marker, ampsize), quote = F, row.names = F, sep = "\t")
  
  #set average noise threshold to default
  average.noise.threshold <- 0.0001
  df2$Reads.to.subtract <- 0
  df2$Adjusted.Reads <- df2$Sample_ASV_Reads - df2$Reads.to.subtract
  df2 <- df2[df2$Adjusted.Reads > 0,]
  
}

#classify each detection as "keep" or "discard" based on average noise threshold and absolute minimum read count
df2 <- ddply(.data = df2, .variables = "Run_Wide_ASV_Name", .fun = transform, seq_mean  = mean(Adjusted.Reads), seq_max = max(Adjusted.Reads))
df2$position_about_mean <- with(df2, ifelse(Adjusted.Reads >= seq_mean, "Above", "Below"))
df2$replicate_count <- ave(seq(nrow(df2)), df2$Run_Wide_ASV_Name, df2$Sample, FUN = length)
df2$temp <- paste(df2$Run_Wide_ASV_Name, df2$Sample, sep = "|")
df2$replicate_count_below_min <- with(df2,ifelse(Adjusted.Reads >= min.reads,0,1))
temp.table <- aggregate(replicate_count_below_min ~ temp, df2, sum)
df2$replicate_count_below_min.total <- temp.table$replicate_count_below_min[match(df2$temp, temp.table$temp)]
df2$replicate_count_adjusted <- df2$replicate_count - df2$replicate_count_below_min.total
df2$replicate_count_about_mean <- ave(seq(nrow(df2)), df2$Run_Wide_ASV_Name, df2$Sample, df2$position_about_mean, FUN = length)
df2$replicate_ratio <- df2$replicate_count_about_mean/df2$replicate_count_adjusted

df2$Noise_Status <- with(df2, ifelse(Adjusted.Reads >= (seq_max*average.noise.threshold) & Adjusted.Reads >= min.reads, "Keep",
                                                       ifelse(Adjusted.Reads < (seq_max*average.noise.threshold) & Adjusted.Reads >= min.reads & replicate_count_adjusted > (num.reps/2), "Keep",
                                                              "Discard")))

# remove noise
df2 <- df2[df2$Noise_Status == "Keep",]

df2 <- df2[,c("Sample", "Replicate", "Sample_ASV_Name", "Adjusted.Reads", "Sample_ASV_Consensus_Sequence",
              "Run_Wide_ASV_Name", "Run_Wide_ASV_Consensus_Sequence", "Tax")]
names(df2)[4] <- "Reads"

write.table(df2, "df.txt", quote = F, row.names = F, sep = "\t")

###############################################################
###################### CHIMERA SCREEN PREP ####################
###############################################################
fasta <- DNAStringSet(df2$Run_Wide_ASV_Consensus_Sequence)
names(fasta) <- paste0(df2$Run_Wide_ASV_Name, ";ab=", df2$Reads)
writeXStringSet(fasta, "temp.fasta")








