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

# Modified function to create the sample/haplotype count plot
create_count_plot <- function(data) {
  counts <- data %>%
    group_by(partition) %>%
    summarize(
      sample_count = n_distinct(sample),
      haplotype_count = n_distinct(haplotype)
    ) %>%
    pivot_longer(
      cols = c(sample_count, haplotype_count),
      names_to = "count_type",
      values_to = "count"
    )
  
  ggplot(counts, aes(x = partition, y = count, fill = count_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(
      values = c("sample_count" = "#2E86C1", "haplotype_count" = "#E74C3C"),
      labels = c("Samples", "Haplotypes")
    ) +
    # Add more x-axis ticks
    scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
    # Add more y-axis ticks
    scale_y_continuous(breaks = function(x) seq(0, ceiling(max(x)), by = 2)) +
    theme_minimal() +
    labs(
      title = "Sample and Haplotype Counts by Partition",
      x = "Partition Number",
      y = "Count",
      fill = "Type"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_line(color = "gray90")
    )
}

# Modified function to create the sequence length plot
create_length_plot <- function(data) {
  lengths <- data %>%
    group_by(partition) %>%
    summarize(total_length = sum(length) / 1e6)  # Convert to megabases
  
  ggplot(lengths, aes(x = partition, y = total_length)) +
    geom_bar(stat = "identity", fill = "#27AE60") +
    # Add more x-axis ticks
    scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)) +
    # Format y-axis in megabases with comma separator
    scale_y_continuous(
      breaks = function(x) pretty(c(0, x), n = 10),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    theme_minimal() +
    labs(
      title = "Total Sequence Length by Partition",
      x = "Partition Number",
      y = "Total Length (Mb)",  # Updated to show Mb
      caption = "Values shown in megabases (Mb)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_line(color = "gray90")
    )
}



directory <- '/home/guarracino/Desktop/Garrison/impg/aaa'

# List all partition BED files
bed_files <- list.files(directory, pattern = "partition\\d+\\.bed$", full.names = TRUE)

# Process each file and combine results
all_data <- map_df(bed_files, process_bed_file)


# Create plots
count_plot <- create_count_plot(all_data)
length_plot <- create_length_plot(all_data)

# Save plots
#ggsave("partition_counts.pdf", count_plot, width = 10, height = 6)
#ggsave("partition_lengths.pdf", length_plot, width = 10, height = 6)

# Combine plots using patchwork
combined_plot <- count_plot / length_plot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Partition Analysis",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

# Save combined plot
ggsave("partition_analysis.pdf", combined_plot, width = 12, height = 10)





