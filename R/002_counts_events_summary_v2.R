# Author: Dan Lou 
# Heidelberg University, ZMBH, Knop Lab
# 2025-09-02


library(tidyverse)
library(dplyr)

# Transform the "fil3_pri_idxstats_summary.csv" table into wide format
# Reading the data from the previous rule's output
combined_counts <- read.csv(snakemake@input[[1]])

# Transform to wide format
df_gene <- combined_counts %>% 
  filter (! (IndexName %in% c("impLbCas12a", "PCRcassette")))  %>%
  separate(IndexName, into = c("Gene", "Component"), sep = "_", fill = "right") %>%
  spread(key = Component, value = Counts)

df_cas <- combined_counts %>% 
  filter (IndexName %in% c("impLbCas12a", "PCRcassette"))  %>%
  spread(key = IndexName, value = Counts)



# Transform the "fil2_indel_idxstats_summary.csv" table into wide format
# Reading the data from the previous rule's output
indel_counts <- read.csv(snakemake@input[[2]])

# Transform to wide format
df_indel <- indel_counts %>% 
  filter (! (IndexName %in% c("impLbCas12a", "PCRcassette")))  %>%
  separate(IndexName, into = c("Gene", "Component"), sep = "_", fill = "right") %>%
  spread(key = Component, value = Counts) %>%
  select(-M1M1, -m2M1) %>%
  rename(Indel = HDR)


# Summarize the read counts from "off_target_summary.tsv" in "readcount_offtarget.tsv" 
# Reading the data from the previous rule's output
offtar_counts <- read.table(snakemake@input[[3]], header = TRUE)

df_offtar <- offtar_counts %>% 
  select(barcode, inner_bc, gene, read_id) %>%
  group_by(barcode, inner_bc, gene) %>%
  summarize(count = n(), .groups = 'drop') %>%
  rename(offtar = count) %>%
  rename(Barcode = barcode) %>%
  rename(Gene = gene)


# Summarize the read counts from "misAssignedM1M1_summary.tsv" in "readcount_misAssignedM1M1.csv"
# Reading the data from the previous rule's output
misM1M1_counts <- read.table(snakemake@input[[4]], header = TRUE)

df_misM1M1 <- misM1M1_counts %>% 
  select(barcode, inner_bc, gene, read_id) %>%
  group_by(barcode, inner_bc, gene) %>%
  summarize(count = n(), .groups = 'drop') %>%
  rename(misM1M1 = count) %>%
  rename(Barcode = barcode) %>%
  rename(Gene = gene)


# combine df_gene and df_indel into one data frame
df <- left_join(df_gene, df_indel) %>%
  left_join(., df_misM1M1) %>%
  mutate(misM1M1 = replace_na(misM1M1, 0)) %>%
  mutate(INDEL = M1M1 + Indel + misM1M1)


### merge with sample metadata
metadata <- read.csv(snakemake@input$METADATA, header = TRUE) %>%
  rename(Gene = gene)

df_combine <- left_join(df_cas, metadata) %>%
  left_join(., df, by = join_by(Barcode, inner_bc, Gene)) %>%
  left_join(., df_offtar, by = join_by(Barcode, inner_bc, Gene)) %>%
  mutate(offtar = replace_na(offtar, 0)) %>%
  mutate(concatemer = PCRcassette + m2M1) %>%
  select(-PCRcassette, -m2M1, -M1M1, -Indel, -misM1M1) %>%
  mutate(total = HDR + INDEL + concatemer + offtar) %>%
  mutate(genomic_junction = HDR + INDEL + offtar)


# Write the transformed data frames to new CSV files
write.csv(df, snakemake@output[[1]], row.names = FALSE)
write.csv(df_cas, snakemake@output[[2]], row.names = FALSE)
write.csv(df_offtar, snakemake@output[[3]], row.names = FALSE)
write.csv(df_misM1M1, snakemake@output[[4]], row.names = FALSE)
write.csv(df_combine, snakemake@output[[5]], row.names = FALSE)
