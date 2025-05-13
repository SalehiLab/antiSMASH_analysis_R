library(rentrez)
library(stringr)
library(dplyr)
library(tidyr)
library(Biostrings)

#####fna files##########
##################input file########################
# Get a list of file paths to all .fna files in the "input" directory
file_list <- list.files(path = "./input", pattern = "*.fna", full.names = TRUE)
# Read the contents of all .fna files into a single character vector
Data0 <- readLines(file_list)
# Split header lines at commas to separate refseq and explanation data
headerLines <- strsplit(grep("^>", Data0, value = TRUE), ",")
headerData<- lapply(headerLines,strsplit, " ")
# Create a data frame to store refseq and explanation data
ref_data <- data.frame(matrix(ncol = 2, nrow = length(headerData)), stringsAsFactors = FALSE)
names(ref_data) <- c("refseq", "explanation")
# Fill in the data frame with extracted information
for (i in 1:length(headerData)) {
  ref_data[i, "refseq"] <- gsub(">", "", headerData[[i]][[1]][[1]])
  ref_data[i, "explanation"] <- paste(headerData[[i]][[1]][-1], collapse = " ")
}

# add a new column for microorganism's sequence
ref_data$refseq <- as.character(ref_data$refseq)
ref_data$sequence <- NA
# find the beginning of a fasta 
Fasta_header <- grep(">",Data0)

# create a loop for sequence extraction
for (i in seq_along(Fasta_header)) {
  start <- Fasta_header[i] + 1
  end <- ifelse(i < length(Fasta_header), Fasta_header[i + 1] - 1, length(Data0))
  
  
  header <- Data0[Fasta_header[i]]
  
  for (j in seq_along(ref_data$refseq)) {
    if (!is.na(ref_data$refseq[j]) && str_detect(header, ref_data$refseq[j])) {
      sequence <- paste(Data0[start:end], collapse = "")
      ref_data$sequence[j] <- sequence
    }
  }
}
# find genome size and GC content of a microorganism
ref_data$G_count <- str_count(ref_data$sequence, "G")
ref_data$C_count <- str_count(ref_data$sequence, "C")
ref_data$Length <- str_length(ref_data$sequence)
ref_data$GC_content <- (ref_data$G_count + ref_data$C_count) / ref_data$Length * 100

# find taxonomy of each microorganism
get_taxonomy <- function(refseq_ids) {
  result <- strsplit(entrez_fetch(db = "nuccore", id = refseq_ids, rettype = "gb", retmode = "text"), "\n")[[1]]
  taxonomy_lines <- grep("^  ORGANISM", result)
  REFERENCE_lines <- grep("REFERENCE", result)
  taxonomy <- gsub("\\s", "", paste0(result[(taxonomy_lines+1 ):(REFERENCE_lines[1]-1)],collapse = ""))
  return(c(refseq_ids,taxonomy))
}
taxonomy<- lapply(ref_data$refseq, get_taxonomy)
#ref_data$taxonomy <- c(taxonomy1, taxonomy2, taxonomy3, taxonomy4, taxonomy5, taxonomy6, taxonomy7)
second_elements <- lapply(taxonomy, function(x) x[2])
ref_data$taxonomy <- second_elements

ref_data <- apply(ref_data,2,as.character)
ref_data <- as.data.frame(ref_data)

# create a csv file containing sequence
write.csv(ref_data,file = "./1_input_data_sequence.csv", row.names = FALSE)

ref_data_ws <- subset(ref_data, select = -sequence)
# create a csv file without sequence(due to the huge size of sequence)          
write.csv(ref_data_ws,file = "./1_input_data.csv" , row.names = FALSE)


##################GBK files########################
# Get a list of file paths for GenBank files with names starting with "merged"
file_list2 <- list.files(pattern = "^merged.*\\.gbk$")
# Read the contents of GenBank files
Data1 <- readGBFF(file_list2)

# Assign accession numbers to the names of each genbank record 
names(Data1) <- sapply(Data1, function(x) {
  version_line <- x[grep("^VERSION", x)]
  gsub("VERSION\\s+", "", version_line)
})
names(Data1)

# Capture the features for each GenBank record
GBFFfeat <- lapply(Data1, function(x) {
  i_start <- grep("^LOCUS ", x)
  i_end <- grep("^ORIGIN", x)
  paste0(x[(i_start[1]):(i_end[1]-1)])
}
)
###################extract protocluster data################
# Apply the extract_protocluster_data function to each element in GBFFfeat
result_list <- mapply(extract_protocluster_data, GBFFfeat, ref_data$refseq, ref_data$explanation, SIMPLIFY = FALSE)
# Combine the resulting data frames into one
combined_result_df <- do.call("rbind", lapply(result_list, function(x) cbind(x[[1]], x[, -1])))
# Rename the column names
names(combined_result_df) <- c("refseq","explanation", "Length","Protocluster_number","FromTo","Product"
                               ,"Category","Core_location","aStool","Contig_edge","cutoff",
                               "Tool","Neighbourhood","Detection_rule")
write.csv(combined_result_df, file= "./2_results_protocluster.csv")

###################extract protocore data################
# Apply the extract_proto_core_data function to each element in GBFFfeat
result_list <- mapply(extract_proto_core_data, GBFFfeat, ref_data$refseq, ref_data$explanation, SIMPLIFY = FALSE)
# Combine the resulting data frames into one
combined_result_df <- do.call("rbind", lapply(result_list, function(x) cbind(x[[1]], x[, -1])))
# Rename the column names
names(combined_result_df) <- c("refseq","explanation", "Length","Protocluster_number","FromTo","Product",
                               "aStool","cutoff","Tool","Neighbourhood","Detection_rule")
write.csv(combined_result_df, file= "./3_results_protocore.csv")
