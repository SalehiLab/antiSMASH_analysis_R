# genome_stats.R
# Genome-level sequence statistics: raw sequence extraction, length, GC
# content, and taxonomy lookup. Depends on get_taxonomy() from ncbi.R.

library(stringr)

#' Extract each genome's raw sequence from its FASTA file, and append
#' sequence, length, GC content, and taxonomy columns to a refseq/Deffinition
#' table.
#'
#' @param ref_data Data frame with columns `refseq`, `Deffinition`, as
#'   returned by parse_fasta_headers().
#' @param fna_files Character vector of paths to the .fna files ref_data was
#'   built from (re-read here to pull sequence lines).
#' @return ref_data with added `sequence`, `Length`, `GC_content`, `taxonomy`
#'   columns (all coerced to character, matching the original pipeline).
add_genome_stats <- function(ref_data, fna_files) {
  Data0 <- readLines(fna_files)
  ref_data$refseq <- as.character(ref_data$refseq)
  ref_data$sequence <- NA

  Fasta_header <- grep(">", Data0)

  for (i in seq_along(Fasta_header)) {
    start <- Fasta_header[i] + 1
    end <- ifelse(i < length(Fasta_header), Fasta_header[i + 1] - 1, length(Data0))
    header <- Data0[Fasta_header[i]]

    for (j in seq_along(ref_data$refseq)) {
      if (!is.na(ref_data$refseq[j]) && str_detect(header, ref_data$refseq[j])) {
        ref_data$sequence[j] <- paste(Data0[start:end], collapse = "")
      }
    }
  }

  G_count <- str_count(ref_data$sequence, "G")
  C_count <- str_count(ref_data$sequence, "C")
  ref_data$Length <- str_length(ref_data$sequence)
  ref_data$GC_content <- (G_count + C_count) / ref_data$Length * 100

  taxonomy <- lapply(ref_data$refseq, get_taxonomy)
  ref_data$taxonomy <- lapply(taxonomy, function(x) x[2])

  ref_data <- apply(ref_data, 2, as.character)
  ref_data <- as.data.frame(ref_data)

  ref_data
}
