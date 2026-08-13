# run_genome_mining.R
# Orchestrates the exploratory genome mining pipeline: protoclusters,
# candidate clusters, PFAM/aSDomain annotations, PKS shortlisting, and
# KnownClusterBlast compound hits.
#
# Expects, in the working directory:
#   ./input/*.fna
#   merged*.gbk
#   ./knownclusterblast/*.txt
#
# This script contains NO function definitions -- all logic lives in R/.

library(dplyr)

source("R/io_read.R")
source("R/ncbi.R")
source("R/extract_features.R")
source("R/filters.R")

cat("Start of run_genome_mining.R\n")

## 1. Input FASTA headers ----------------------------------------------------
fna_files <- list.files(path = "./input", pattern = "*.fna", full.names = TRUE)
ref_data <- parse_fasta_headers(fna_files)
write.csv(ref_data, file = "./1_input_data.csv")

## 2. Merged GenBank features -------------------------------------------------
GBFFfeat <- load_merged_gbff()

## 3. Protocluster data --------------------------------------------------------
result_list <- mapply(extract_protocluster_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined_result_df <- do.call("rbind", result_list)
write.csv(combined_result_df, file = "./2_results_protocluster.csv")

## 4. Candidate cluster data ---------------------------------------------------
result_list_cand_cluster <- mapply(extract_cand_cluster_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined_result_cand_cluster <- do.call("rbind", result_list_cand_cluster)
write.csv(combined_result_cand_cluster, file = "./3_results_cand_cluster.csv", row.names = FALSE)

## 5. PKS shortlist -------------------------------------------------------------
selected_rows <- select_pks_candidates(combined_result_cand_cluster)
write.csv(selected_rows, file = "./3_results_cand_cluster_Selected.csv", row.names = FALSE)

## 6. PFAM domain data -----------------------------------------------------------
result_list_PFAM_data <- mapply(extract_PFAM_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined_result_PFAM_data <- do.call("rbind", result_list_PFAM_data)
write.csv(combined_result_PFAM_data, file = "./4_results_PFAM_data.csv", row.names = FALSE)

## 7. aSDomain data ----------------------------------------------------------------
result_list_aSDomain_data <- mapply(extract_aSDomain_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined_result_aSDomain_data <- do.call("rbind", result_list_aSDomain_data)
write.csv(combined_result_aSDomain_data, file = "./5_results_aSDomain_data.csv", row.names = FALSE)

## 8. KnownClusterBlast compound hits ----------------------------------------------
Compounds_data <- read_knownclusterblast("./knownclusterblast/")

## 9. Merge protocluster + compound data -------------------------------------------
merged_data <- merge(combined_result_df, Compounds_data,
  by = c("refseq", "Protocluster_number"), sort = FALSE
)
write.csv(merged_data, file = "./2_mergedResults.csv")

cat("End of run_genome_mining.R\n")
