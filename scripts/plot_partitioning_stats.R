library(tidyverse)
library(ggplot2)
library(stringr)
library(patchwork)  # For combining plots
options(scipen=10000)

# Function to read and process a single BED file
process_bed_file <- function(file_path) {
  # Extract partition number from filename
  partition_num <- as.numeric(str_extract(basename(file_path), "\\d+"))
  
  # Read BED file
  bed_data <- read_tsv(file_path, 
                       col_names = c("name", "start", "end", "score", "strand", "orientation"),
                       col_types = "ciiccc")
  
  # Calculate interval lengths
  bed_data <- bed_data %>%
    mutate(length = end - start)
  
  # Extract sample and haplotype information
  bed_data <- bed_data %>%
    mutate(
      sample = str_extract(name, "^[^#]+"),
      haplotype = str_extract(name, "^[^#]+#([^#]+)#", group = 1)
    ) %>%
    mutate(
      haplotype = paste0(sample, "#", haplotype)
    )
  
  # Add partition information
  bed_data$partition <- partition_num
  
  return(bed_data)
}

# Modified function to create separate sample and haplotype count plots
create_sample_plot <- function(data) {
  sample_counts <- data %>%
    group_by(partition) %>%
    summarize(count = n_distinct(sample))
  
  # Calculate y-axis breaks with smaller intervals
  max_count <- max(sample_counts$count)
  y_breaks <- seq(0, ceiling(max_count), by = max(1, ceiling(max_count/20)))
  
  ggplot(sample_counts, aes(x = factor(partition), y = count)) +
    geom_bar(stat = "identity", fill = "#2E86C1", width = 0.8) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal() +
    labs(
      title = "Sample Counts by Partition",
      x = "Partition Number",
      y = "Number of Samples"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.minor = element_line(color = "gray90")
    )
}

create_haplotype_plot <- function(data) {
  haplotype_counts <- data %>%
    group_by(partition) %>%
    summarize(count = n_distinct(haplotype))
  
  # Calculate y-axis breaks with smaller intervals
  max_count <- max(haplotype_counts$count)
  y_breaks <- seq(0, ceiling(max_count), by = max(1, ceiling(max_count/20)))
  
  ggplot(haplotype_counts, aes(x = factor(partition), y = count)) +
    geom_bar(stat = "identity", fill = "#E74C3C", width = 0.8) +
    scale_y_continuous(
      breaks = y_breaks,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal() +
    labs(
      title = "Haplotype Counts by Partition",
      x = "Partition Number",
      y = "Number of Haplotypes"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.minor = element_line(color = "gray90")
    )
}

create_length_plot <- function(data) {
  lengths <- data %>%
    group_by(partition) %>%
    summarize(total_length = sum(length) / 1e6)  # Convert to megabases
  
  # Calculate y-axis breaks with smaller intervals
  max_length <- max(lengths$total_length)
  y_breaks <- seq(0, ceiling(max_length), by = max(1, ceiling(max_length/20)))
  
  ggplot(lengths, aes(x = factor(partition), y = total_length)) +
    geom_bar(stat = "identity", fill = "#27AE60", width = 0.8) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) format(x, big.mark = ",", scientific = FALSE),
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal() +
    labs(
      title = "Total Sequence Length by Partition",
      x = "Partition Number",
      y = "Total Length (Mb)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.grid.minor = element_line(color = "gray90")
    )
}

directory <- '/home/guarracino/Desktop/Garrison/impg/'

# List all partition BED files
bed_files <- list.files(directory, pattern = "partition\\d+\\.bed$", full.names = TRUE)

# Process each file and combine results
all_data <- map_df(bed_files, process_bed_file)

# Create individual plots
sample_plot <- create_sample_plot(all_data)
haplotype_plot <- create_haplotype_plot(all_data)
length_plot <- create_length_plot(all_data)

# Combine plots using patchwork
combined_plot <- sample_plot / haplotype_plot / length_plot +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Pangenome Partitions",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

# Save combined plot
ggsave("partition_analysis.pdf", combined_plot, width = 30, height = 12)
