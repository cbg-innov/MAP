##### MAPLE v 1.0 #####

# Get arguments from bash
args <- commandArgs(trailingOnly = TRUE)
wkdir <- args[1]
UMImap_path <- args[2]
runid <- args[3]

setwd(wkdir)

# Load libraries
library("data.table")
library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(patchwork)

# Load mapping file
UMImap <- read.csv(UMImap_path, header = TRUE, sep = "\t")

# Create helper column to match FASTA filenames
UMImap$FASTA_File <- paste0(UMImap$Sample, "_", UMImap$Marker, "_", UMImap$Target.amplicon.length, "bp.fasta")

# Count sequences in each FASTA
fasta_files <- list.files(pattern = "*.fasta")
output <- data.frame()

for(f in fasta_files){
  fasta <- readDNAStringSet(f)
  seq_count <- length(fasta)
  
  # Ensure this FASTA exists in mapping helper column
  idx <- which(UMImap$FASTA_File == f)
  if(length(idx) == 0){
    warning(paste("FASTA file", f, "not found in mapping file helper column"))
    next
  }
  
  temp_df <- data.frame(
    "FASTA_File" = f,
    "Seq_Count"  = seq_count
  )
  output <- rbind(output, temp_df)
}

# Merge sequence counts with mapping info
df <- merge(
  UMImap[, c("Plate", "Well", "FASTA_File")],
  output,
  by = "FASTA_File",
  all.x = TRUE
)

# Extract original Sample name
df <- df %>%
  left_join(UMImap %>% select(FASTA_File, Sample),
            by = "FASTA_File")

# Prepare heatmap data
df_clean <- df %>%
  mutate(
    Row = str_sub(Well, 1, 1),
    Col = as.integer(str_sub(Well, 2))
  )

row_levels <- LETTERS[1:8]
col_levels <- 1:12

all_grid <- tidyr::expand_grid(
  Plate = unique(df_clean$Plate),
  Row   = row_levels,
  Col   = col_levels
)

df_plot <- all_grid %>%
  left_join(df_clean, by = c("Plate", "Row", "Col")) %>%
  mutate(Row = factor(Row, levels = rev(LETTERS[1:8])))

df_plot$Seq_Count[is.na(df_plot$Seq_Count) & !is.na(df_plot$FASTA_File)] <- 0

# Compute breaks for heatmap
max_val <- max(df_plot$Seq_Count, na.rm = TRUE)
breaks <- c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
breaks <- breaks[breaks <= max_val]

# Plot heatmap
plot1 <- ggplot(df_plot, aes(x = Col, y = Row, fill = Seq_Count)) +
  geom_tile(color = "grey70") +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = pseudo_log_trans(base = 10, sigma = 1),
    na.value = "white", 
    labels = comma,
    breaks = breaks
  ) +
  facet_wrap(~Plate) +
  coord_fixed() +
  labs(x = element_blank(), y = element_blank(), fill = "Reads", title = "Demultiplexed reads by well") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 15)))

# Read counts histogram
df_reads <- read.table(paste0(runid,"_readcounts.txt"), sep = "\t", header = FALSE)

total_rows <- nrow(df_reads)
for (i in 1:total_rows) {
  if (i == 1) {
    df_reads$V3[i] <- paste0("(", 100, "% remaining)")
    df_reads$V4[i] <- paste0("= ", 0,"% drop")
    df_reads$V5[i] <- 0
  } else {
    df_reads$V3[i] <- paste0("(", round(df_reads$V2[i] / df_reads$V2[1] * 100, digits = 0), "% remaining)")
    df_reads$V4[i] <- paste0("= ", 100 - round(df_reads$V2[i]/df_reads$V2[i-1]*100, digits = 0), "% drop")
    df_reads$V5[i] <- 100 - round(df_reads$V2[i] / df_reads$V2[i-1] * 100, digits = 0)
  }
}

plot2 <- ggplot(df_reads, aes(y = V2, x = factor(V1, levels = V1), fill = "#66BB6A")) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Demultiplexed read retention") +
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, color = "black", face = "bold"),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 15))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_identity() +
  geom_text(aes(label=paste(format(V2, big.mark = ",", trim = TRUE), V4, V3, "",sep = "\n"), fontface = "bold"),
            vjust=0,
            cex = 3) +
  coord_cartesian(ylim = c(0,max(df_reads$V2) *1.15))

# export to PDF
pdf(sprintf("Demultiplexing_Results_%s.pdf", runid), width = 11, height = 8.5)
print(plot1)
plot2_scaled <- plot2 + theme(plot.margin = margin(0, 0, 0, 0))
centered_plot2 <- plot_spacer() | plot2_scaled | plot_spacer()
centered_plot2 <- centered_plot2 + plot_layout(widths = c(0.17, 0.66, 0.17)) # middle = 66%

final_plot <- (plot_spacer() / centered_plot2 / plot_spacer()) +
  plot_layout(heights = c(0.2, 0.825, 0.2))

print(final_plot)

dev.off()
