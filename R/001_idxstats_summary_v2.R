# compare with the previous version, here process two separate sets of inputs

library(dplyr)

# Function to read and process a single idxstats file
process_idxstats <- function(file) {
  read_counts_data <- read.table(file, header = FALSE) %>% filter(V1 != "*")
  
  print(head(read_counts_data)) # Debug: Print first few rows of the read data
  
  # Extracting index name (first column) and read counts (last column)
  index_name <- read_counts_data$V1
  read_counts <- read_counts_data$V3

  # Extract the barcode from the file path
  path_elements <- strsplit(file, "/")[[1]]
  barcode <- path_elements[length(path_elements) - 1]
  innerbc <- stringr::str_extract(basename(file), "[A-Z][0-9]{2}")

  data.frame(Barcode = barcode, IndexName = index_name, inner_bc = innerbc, Counts = read_counts)
}

# Apply the function to all fil1.fil2.fil3.pri.idxstats files
files_pri <- snakemake@input$pri_idxstats
counts_pri <- lapply(files_pri, process_idxstats)
# Combine the results into a data frame
combined_counts_pri <- do.call(rbind, counts_pri)

# Apply the function to all fil1.fil2.indel.idxstats files
files_indel <- snakemake@input$indel_idxstats
counts_indel <- lapply(files_indel, process_idxstats)
# Combine the results into a data frame
combined_counts_indel <- do.call(rbind, counts_indel)



# Write to CSV for both sets of inputs
write.csv(combined_counts_pri, snakemake@output[[1]], row.names = FALSE)
write.csv(combined_counts_indel, snakemake@output[[2]], row.names = FALSE)
