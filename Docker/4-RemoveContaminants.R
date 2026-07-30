#get argument from bash script (run name)
args <- commandArgs(trailingOnly = TRUE)
wkdir <- args[1]
runid <- args[2]
marker <- args[3]
ampsize <-as.numeric(args[4])
min_rep_prop <-as.numeric(args[5])# proportion of replicates required
high_dens_prop <-as.numeric(args[6]) # proportion within high density area
alpha_default<-as.numeric(args[7])#Default alpha (if not using negative controls)
alpha_quantile <- as.numeric(args[8]) #quantile percentage used
k_multiplier   <- as.numeric(args[9])
mapping_file   <- args[10]   # parameters/mapping file containing the 'Negative Control' column

setwd(wkdir)
set.seed(123)

# ----------------------------
# LOAD REQUIRED PACKAGES
# ----------------------------
library(dplyr)
library(Biostrings)
library(purrr)
library(classInt)

# ----------------------------
# IMPORT DATA
# ----------------------------

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

# import master OTU table
df_raw <- read.table("all_otus_run_wide_OTU_info.txt", header = TRUE, sep = "\t")

# add sintax results to master OTU table
df <- merge(df_raw, sintax.input[,c(1,4)], by.x = "Run_Wide_OTU_Name", by.y = "Query", all.x = TRUE)
# ----------------------------
# Negative controls: read from the 'Negative Control' column of the parameters
# file (NOT from the sample name). A well is a negative control if that column
# is 'yes' (any case) or '1'; blank / 'no' / '0' -> sample.
# ----------------------------
neg_samples <- character(0)
if (length(mapping_file) == 1 && !is.na(mapping_file) && file.exists(mapping_file)) {
  map <- read.delim(mapping_file, header = TRUE, sep = "\t",
                    check.names = FALSE, colClasses = "character", fill = TRUE)
  if (all(c("Sample", "Negative Control") %in% names(map))) {
    is_neg <- toupper(trimws(as.character(map[["Negative Control"]]))) %in% c("YES", "1")
    # match the master table's Sample: strip _Rep<n> and replace '#' with '-'
    base_samp <- gsub("#", "-", sub("_Rep[0-9]+$", "", as.character(map[["Sample"]])))
    neg_samples <- unique(base_samp[is_neg])
  } else {
    warning("Mapping file lacks 'Sample' and/or 'Negative Control' column; treating all wells as samples.")
  }
} else {
  warning("Mapping file not found; treating all wells as samples.")
}

df$Well_Type <- ifelse(df$Sample %in% neg_samples, "neg", "sample")

# ----------------------------
# DETERMINE TOTAL REPLICATES
# ----------------------------
total_reps <- df %>%
  filter(Well_Type == "sample") %>%
  pull(Replicate) %>%
  unique() %>%
  length()

# ----------------------------
# STEP 1: total abundance per OTU across plate
# ----------------------------
otu_totals <- df %>%
  filter(Well_Type == "sample") %>%
  group_by(Run_Wide_OTU_Name) %>%
  dplyr::summarise(T_i = sum(Sample_OTU_Reads), .groups = "drop")

# ----------------------------
# STEP 2: estimate contamination rate (alpha)
# ----------------------------

neg_rates <- df %>%
  filter(Well_Type == "neg") %>%
  group_by(Run_Wide_OTU_Name) %>%
  dplyr::summarise(neg_reads = sum(Sample_OTU_Reads), .groups = "drop") %>%
  left_join(otu_totals, by = "Run_Wide_OTU_Name") %>%
  mutate(r_i = neg_reads / T_i) %>%
  filter(!is.na(r_i), T_i > 0)

# temporarily remove OTUs that are disproportionately highly contaminated in control relative to sample
neg_rates_filt<- neg_rates %>%
  filter(r_i < 1)

if(nrow(neg_rates_filt) == 0){
  alpha <- alpha_default
}else{
  alpha <- quantile(neg_rates_filt$r_i, probs = alpha_quantile, na.rm = TRUE)
}


# ----------------------------
# STEP 3: replicate support + IQR consistency
# ----------------------------

### requires library(classInt) and library(purrr)
#First calculate kernel density, using outlier-resistant"Sheather-Jones" alg.
#Second, find the number of  Jenks breaks in the kernel density dist using #replicates+1
#Third, calculate the points of kernel density x-axis where densest region borders are defined

get_hot_zone <- function(reads, n_breaks) {
  reads <- as.vector(reads)
  
  if (length(reads) < 3 || n_breaks < 3) {
    return(list(list(lower = NA, upper = NA, prop = NA)))
  }
  
  orig_reads <- reads
  
  if (length(unique(reads)) < 2 || sd(reads) < 0.01) {
    reads <- jitter(reads, amount = 0.005)
  }
  
  # Adaptive bandwidth: widen SJ when data is tight (low CV)
  # to prevent over-segmentation of naturally clustered values
  bw_sj <- tryCatch(bw.SJ(reads), error = function(e) bw.nrd0(reads))
  cv <- sd(reads) / abs(mean(reads))
  spread_factor <- dplyr::case_when(
    cv < 0.10 ~ 3.0,
    #increase value in the next line (cv < X) to be more forgiving of a spread-out distribution for middle value data.
    cv < 0.3 ~ 2.0,
    cv < 0.40 ~ 1.3,
    TRUE      ~ 1.0
  )
  # Floor scales with both range and mean — prevents over-segmentation
  # of low-magnitude vectors where the range is naturally small
  #Adjust abs(mean(reads)) (higher than 0.05) to be more aggressive at values consistently at lower end.
  bw_floor <- max(diff(range(reads)) * 0.05, abs(mean(reads)) * 0.05)
  bw_adaptive <- max(bw_sj * spread_factor, bw_floor)
  
  d      <- density(reads, bw = bw_adaptive)
  breaks <- classIntervals(d$y, n = n_breaks, style = "jenks")
  hot_x  <- d$x[d$y >= tail(breaks$brks, 2)[1]]
  
  lower <- max(round(hot_x[1], digits = 0), 0)
  upper <- round(tail(hot_x, 1), digits = 0)
  if (upper == lower) upper <- lower + 1
  prop  <- mean(orig_reads >= lower & orig_reads <= upper)
  
  list(list(lower = lower, upper = upper, prop = prop))
}


# then run the pipe
rep_support <- df %>%
  filter(Well_Type == "sample") %>%
  group_by(Run_Wide_OTU_Name, Sample) %>%
  summarise(
    reps_present    = n_distinct(Replicate),
    total_reads     = sum(Sample_OTU_Reads),
    hz              = get_hot_zone(Sample_OTU_Reads, reps_present),
    .groups         = "drop"
  ) %>%
  mutate(
    lower_bound     = map_dbl(hz, "lower"),
    upper_bound     = map_dbl(hz, "upper"),
    prop_within_IQR = map_dbl(hz, "prop")
  ) %>%
  select(-hz)


# ----------------------------
# STEP 4: classify checks + action
# ----------------------------

decision_table <- rep_support %>%
  left_join(otu_totals, by = "Run_Wide_OTU_Name") %>%
  left_join(neg_rates[-3], by = "Run_Wide_OTU_Name") %>%
  mutate(
    expected_contam = alpha * T_i,
    
    Contam_Check = ifelse(
      total_reads > (k_multiplier * expected_contam) &
        #remove Run_Wide OTUs from sample that are highly contaminated in controls  
        (is.na(r_i) | r_i <= 1),
      "pass", "fail"
    ),
    
    Replicate_Check = ifelse(
      (reps_present / total_reps >= min_rep_prop) &
        (prop_within_IQR >= high_dens_prop) , 
      "pass", "fail"
    ),
    
    #r_i over 1 is an automatic discard
    Action = case_when(
      !is.na(r_i) & r_i > 1                             ~ "discard",
      Contam_Check == "pass" | Replicate_Check == "pass" ~ "retain",
      TRUE                                               ~ "discard"
    )
  ) %>%
  select(
    Run_Wide_OTU_Name, Sample,
    total_reads, reps_present,
    total_reads, lower_bound, upper_bound, prop_within_IQR,
    expected_contam, 
    Contam_Check, Replicate_Check, Action,
    neg_reads, r_i, T_i
  )

# ----------------------------
# STEP 5: attach to full dataset
# ----------------------------
df_with_action <- df %>%
  left_join(decision_table, by = c("Run_Wide_OTU_Name", "Sample"))

# Label negatives explicitly
df_with_action$Action[is.na(df_with_action$Action)] <- "retain"
df_with_action$Contam_Check[is.na(df_with_action$Contam_Check)] <- "NA"
df_with_action$Replicate_Check[is.na(df_with_action$Replicate_Check)] <- "NA"

#Create table for retained sequences (negative-control wells are excluded here)
df_retain <- df_with_action[which(df_with_action$Action=="retain" & df_with_action$Well_Type!="neg"),c(2:6,1,8:9)]

#Create table for negative controls
df_neg <- df_with_action[which(df_with_action$Well_Type=="neg"),c(2:6,1,8:9)]

#Create retained sequences fasta for CHIMERA CHECK
fasta <- DNAStringSet(df_retain$Run_Wide_OTU_Consensus_Sequence)
names(fasta) <- paste0(df_retain$Run_Wide_OTU_Name, ";ab=", df_retain$Reads)

#create table for contaminants reads higher in control than in sample
df_ri_1plus <- df_with_action[which(df_with_action$r_i >1),]
write.table(
  df_ri_1plus,
  file = "df_contaminants_removed.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ----------------------------
# OUTPUTS
# ----------------------------

# Filtered fasta and table to be used further:
writeXStringSet(fasta, "temp.fasta")

write.table(
  df_retain,
  file = "df.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#Table showing sequence retention/rejection
write.table(
  df_with_action,
  file = "df_with_action.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



#Negative controls table
if (length(df_neg[,1]) != 0){
  #extract OTUs in negative control wells and output for user's files
  write.table(df_neg, sprintf("%s_%s_%sbp_NegativeControlOTUs.tsv", runid, marker, ampsize), quote = F, row.names = F, sep = "\t")
} else { #if no negative control wells
  #extract OTUs in negative control wells and output for user's files
  write.table(df_neg[0,], sprintf("%s_%s_%sbp_NegativeControlOTUs.tsv", runid, marker, ampsize), quote = F, row.names = F, sep = "\t")
}
