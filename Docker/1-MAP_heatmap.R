##### MAP #####

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
UMImap$FASTA_File <- paste0(UMImap$Sample, "_", UMImap$Marker, "_", UMImap$Target.Amplicon.Length, "bp.fasta")

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

# Prepare heatmap data. Parse with regex rather than fixed positions so that
# A1 / A01 / a01 all work.
df_clean <- df %>%
  mutate(
    Row = toupper(str_extract(Well, "^[A-Za-z]+")),
    Col = as.integer(str_extract(Well, "[0-9]+$"))
  )

#### Detect plate layout
obs_rows <- suppressWarnings(max(match(df_clean$Row, LETTERS), na.rm = TRUE))
obs_cols <- suppressWarnings(max(df_clean$Col, na.rm = TRUE))
if (!is.finite(obs_rows)) obs_rows <- 1
if (!is.finite(obs_cols)) obs_cols <- 1

#### rows x cols for 6, 12, 24, 48, 96 and 384-well plates
plate_formats <- list(c(2, 3), c(3, 4), c(4, 6), c(6, 8), c(8, 12), c(16, 24))
fits <- Filter(function(f) obs_rows <= f[1] && obs_cols <= f[2], plate_formats)
fmt  <- if (length(fits) > 0) fits[[1]] else c(obs_rows, obs_cols)

row_levels <- LETTERS[1:fmt[1]]
col_levels <- 1:fmt[2]

cat(sprintf("Demultiplexing heatmap: %d-well layout (%s-%s, 1-%d)\n",
            fmt[1] * fmt[2], row_levels[1], row_levels[fmt[1]], fmt[2]))

#### Warn instead of silently dropping wells outside grid
outside <- df_clean[!(df_clean$Row %in% row_levels) | !(df_clean$Col %in% col_levels), ]
if (nrow(outside) > 0) {
  warning(sprintf("%d well(s) fall outside the detected %d-well layout and will not be plotted: %s",
                  nrow(outside), fmt[1] * fmt[2],
                  paste(unique(paste0(outside$Row, outside$Col)), collapse = ", ")))
}

all_grid <- tidyr::expand_grid(
  Plate = unique(df_clean$Plate),
  Row   = row_levels,
  Col   = col_levels
)

df_plot <- all_grid %>%
  left_join(df_clean, by = c("Plate", "Row", "Col")) %>%
  mutate(Row = factor(Row, levels = rev(row_levels)))

df_plot$Seq_Count[is.na(df_plot$Seq_Count) & !is.na(df_plot$FASTA_File)] <- 0

# Compute breaks for heatmap
max_val <- max(df_plot$Seq_Count, na.rm = TRUE)
breaks <- c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
breaks <- breaks[breaks <= max_val]

# Plot heatmap
#### Page and text sizing. facet_wrap lays panels out in roughly a sqrt(n) x sqrt(n)
#### grid, so the page grows with the grid to keep each panel physically usable.
#### Growth is capped so the PDF stays openable; text then tracks the panel size
#### the page actually delivers, rather than shrinking twice over.
n_plates   <- length(unique(df_plot$Plate))
facet_ncol <- ceiling(sqrt(max(n_plates, 1)))
facet_nrow <- ceiling(max(n_plates, 1) / facet_ncol)

pdf_width  <- min(30, max(11,  5.50 * facet_ncol))
pdf_height <- min(24, max(8.5, 4.25 * facet_nrow))

#### Panel width relative to a single-plate page, and page growth vs the 11x8.5 default
panel_scale <- min(1, (pdf_width / facet_ncol) / 11)
page_scale  <- pdf_width / 11

strip_size <- max(7, round(20 * panel_scale))
axis_size  <- max(4, round(11 * panel_scale * (12 / fmt[2])))
wrap_width <- max(10, round(26 * panel_scale))
tile_line  <- if (n_plates > 4 || fmt[2] > 12) 0.1 else 0.25

plot1 <- ggplot(df_plot, aes(x = Col, y = Row, fill = Seq_Count)) +
  geom_tile(color = "grey70", linewidth = tile_line) +
  scale_x_continuous(breaks = col_levels, expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = pseudo_log_trans(base = 10, sigma = 1),
    na.value = "white", 
    labels = comma,
    breaks = breaks
  ) +
  facet_wrap(~Plate, labeller = label_wrap_gen(width = wrap_width)) +
  coord_fixed() +
  labs(x = element_blank(), y = element_blank(), fill = "Reads", title = "Demultiplexed reads by well") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = axis_size),
        strip.text = element_text(size = strip_size, face = "bold", lineheight = 0.9),
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
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) +
  theme_bw() +
  labs(title = "Demultiplexed read retention") +
  theme(axis.text.x = element_text(size = 10 * page_scale, angle = 0, hjust = 0.5, color = "black", face = "bold"),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 22 * page_scale, face = "bold", hjust = 0.5, margin = margin(b = 15))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_identity() +
  geom_text(aes(label=paste(format(V2, big.mark = ",", trim = TRUE), V4, V3, "",sep = "\n"), fontface = "bold"),
            vjust=0,
            cex = 3 * page_scale) +
  coord_cartesian(ylim = c(0,max(df_reads$V2) *1.15))

# export to PDF
pdf(sprintf("Demultiplexing_Results_%s.pdf", runid), width = pdf_width, height = pdf_height)
print(plot1)
plot2_scaled <- plot2 + theme(plot.margin = margin(0, 0, 0, 0))
centered_plot2 <- plot_spacer() | plot2_scaled | plot_spacer()
centered_plot2 <- centered_plot2 + plot_layout(widths = c(0.17, 0.66, 0.17)) # middle = 66%

final_plot <- (plot_spacer() / centered_plot2 / plot_spacer()) +
  plot_layout(heights = c(0.2, 0.825, 0.2))

print(final_plot)

dev.off()
