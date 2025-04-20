# Step 5: CNV Visualization and Reporting

library(tidyverse)
library(ggplot2)

# Set directory containing CNVkit output
cnv_dir <- "../results/cnv_reports"
plot_dir <- "../results/cnv_plots"
dir.create(plot_dir, showWarnings = FALSE)

# Get all .cns files (segments)
cns_files <- list.files(cnv_dir, pattern = "*.cns$", full.names = TRUE)

for (cns_file in cns_files) {
  # Extract sample name
  sample <- tools::file_path_sans_ext(basename(cns_file))

  # Read .cns file
  cns <- read.delim(cns_file)

  # Basic plot: log2 copy number across chromosomes
  p <- ggplot(cns, aes(x = start, y = log2)) +
    geom_point(size = 0.5, color = "steelblue") +
    facet_wrap(~chrom, scales = "free_x", ncol = 6) +
    labs(title = paste("CNV Segments -", sample),
         x = "Genomic Position",
         y = "log2(Copy Ratio)") +
    theme_minimal()

  # Save plot
  ggsave(filename = file.path(plot_dir, paste0(sample, "_segments_plot.pdf")),
         plot = p, width = 12, height = 6)

  message("Saved plot for: ", sample)
}

cat("All CNV plots saved in results/cnv_plots\n")

