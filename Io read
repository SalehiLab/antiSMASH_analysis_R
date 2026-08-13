# io_read.R
# Low-level file reading utilities: antiSMASH merged GenBank files, FASTA
# headers, and KnownClusterBlast text reports. No antiSMASH-domain parsing
# lives here -- see extract_features.R for that.

#' Read a (possibly multi-record) GenBank flat file and split it into one
#' character vector per record, using each "LOCUS" line as a record boundary.
#'
#' @param file Path to a .gbk/.gbff file.
#' @param text Optionally, pre-loaded lines of the file.
#' @param verbose Unused; reserved for future logging.
#' @return A list of character vectors, one per GenBank record.
readGBFF <- function(file, text = readLines(file), verbose = FALSE) {
  if (is(text, "character")) {
    split_indices <- cumsum(grepl("^LOCUS", text))
    text <- base::split(text, split_indices)
  }
  return(text)
}

#' Parse the ">accession, description" header lines of one or more FASTA
#' (.fna) files into a refseq / Deffinition lookup table.
#'
#' @param fna_files Character vector of paths to .fna files.
#' @return A data frame with columns `refseq` and `Deffinition`.
parse_fasta_headers <- function(fna_files) {
  Data0 <- readLines(fna_files)
  headerLines <- strsplit(grep("^>", Data0, value = TRUE), ",")
  headerData <- lapply(headerLines, strsplit, " ")

  ref_data <- data.frame(matrix(ncol = 2, nrow = length(headerData)), stringsAsFactors = FALSE)
  names(ref_data) <- c("refseq", "Deffinition")

  for (i in seq_along(headerData)) {
    ref_data[i, "refseq"] <- gsub(">", "", headerData[[i]][[1]][[1]])
    ref_data[i, "Deffinition"] <- paste(headerData[[i]][[1]][-1], collapse = " ")
  }

  ref_data
}

#' Load all antiSMASH "merged*.gbk" files in the working directory, name each
#' record by its VERSION accession, and extract the feature block between
#' LOCUS and ORIGIN for each record.
#'
#' @param pattern Filename pattern identifying merged GenBank files.
#' @return A named list of character vectors (the GBK feature block for each
#'   genome), named by accession.
load_merged_gbff <- function(pattern = "^merged.*\\.gbk$") {
  file_list <- list.files(pattern = pattern)
  Data1 <- readGBFF(file_list)

  names(Data1) <- sapply(Data1, function(x) {
    version_line <- x[grep("^VERSION", x)]
    gsub("VERSION\\s+", "", version_line)
  })

  GBFFfeat <- lapply(Data1, function(x) {
    i_start <- grep("^LOCUS ", x)
    i_end <- grep("^ORIGIN", x)
    paste0(x[(i_start[1]):(i_end[1] - 1)])
  })

  GBFFfeat
}

#' Read antiSMASH KnownClusterBlast .txt reports from a directory and extract
#' the "Source:" compound name reported under each cluster's "Details:"
#' section.
#'
#' @param dir Path to the knownclusterblast directory.
#' @return A data frame with columns `refseq`, `Protocluster_number`, `Compounds`.
#'   (Note: column renamed from the original `protocluster_number` to
#'   `Protocluster_number` so it matches the canonical
#'   `extract_protocluster_data()` output for merging.)
read_knownclusterblast <- function(dir = "./knownclusterblast/") {
  file_list <- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)
  file_names <- list.files(path = dir, pattern = "\\.txt$")
  refseq <- strsplit(file_names, "_")
  file_contents <- lapply(file_list, readLines)

  Compounds <- lapply(file_contents, function(x) {
    details_line_index <- grep("Details:", x, fixed = TRUE)
    source_lines_index <- details_line_index + 1:4
    source_lines <- x[source_lines_index]
    result <- grep("Source:", source_lines, fixed = TRUE, value = TRUE)
    ifelse(length(result) > 0, sub("\\s", "", sub("Source:", "", result)), NA)
  })

  Compounds_data <- data.frame(matrix(ncol = 3, nrow = length(file_list)), stringsAsFactors = FALSE)
  names(Compounds_data) <- c("refseq", "Protocluster_number", "Compounds")

  for (i in seq_along(file_list)) {
    Compounds_data[i, "refseq"] <- paste(refseq[[i]][-3], collapse = "_")
    Compounds_data[i, "Protocluster_number"] <- gsub("^c|\\.txt$", "", refseq[[i]][3])
    Compounds_data[i, "Compounds"] <- Compounds[[i]]
  }

  Compounds_data
}
