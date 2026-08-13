# run_data_extraction.R
# Orchestrates the genome-level summary (sequence, GC%, taxonomy) and full
# protocluster / proto_core annotation extraction.
#
# Expects, in the working directory:
#   ./input/*.fna
#   merged*.gbk
#
# This script contains NO function definitions -- all logic lives in R/.

library(rentrez)
library(stringr)
library(dplyr)
library(tidyr)
library(Biostrings)

source("R/io_read.R")
source("R/ncbi.R")
source("R/genome_stats.R")
source("R/extract_features.R")

## 1. Input FASTA headers + sequence stats + taxonomy ---------------------------
fna_files <- list.files(path = "./input", pattern = "*.fna", full.names = TRUE)
ref_data <- parse_fasta_headers(fna_files)
ref_data <- add_genome_stats(ref_data, fna_files)

write.csv(ref_data, file = "./1_input_data_sequence.csv", row.names = FALSE)

ref_data_ws <- subset(ref_data, select = -sequence)
write.csv(ref_data_ws, file = "./1_input_data.csv", row.names = FALSE)

## 2. Merged GenBank features ----------------------------------------------------
GBFFfeat <- load_merged_gbff()

## 3. Protocluster data -----------------------------------------------------------
result_list <- mapply(extract_protocluster_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined_result_df <- do.call("rbind", result_list)
write.csv(combined_result_df, file = "./2_results_protocluster.csv")

## 4. Proto-core data ---------------------------------------------------------------
result_list <- mapply(extract_proto_core_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined_result_df <- do.call("rbind", result_list)
write.csv(combined_result_df, file = "./3_results_protocore.csv")
