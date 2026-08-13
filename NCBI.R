# ncbi.R
# Helpers for querying NCBI via rentrez.

library(rentrez)

#' Retrieve the full taxonomic lineage for one or more NCBI RefSeq
#' accessions, by fetching the GenBank record and parsing its ORGANISM
#' field.
#'
#' @param refseq_ids A RefSeq accession (or vector of accessions).
#' @return A character vector: c(refseq_ids, taxonomy).
get_taxonomy <- function(refseq_ids) {
  result <- strsplit(
    entrez_fetch(db = "nuccore", id = refseq_ids, rettype = "gb", retmode = "text"),
    "\n"
  )[[1]]
  taxonomy_lines <- grep("^  ORGANISM", result)
  REFERENCE_lines <- grep("REFERENCE", result)
  taxonomy <- gsub(
    "\\s", "",
    paste0(result[(taxonomy_lines + 1):(REFERENCE_lines[1] - 1)], collapse = "")
  )
  c(refseq_ids, taxonomy)
}
