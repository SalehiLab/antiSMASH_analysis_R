
cat("Start of GenomeMining.R script\n")
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
# Display the ref_data data frame
write.csv(ref_data, file= "./1_input_data.csv")
##################GBK file########################
# Define a function to read a GenBank file
readGBFF <- function (file, text = readLines(file), verbose = FALSE){
  if (is(text, "character")) {
    split_indices <- cumsum(grepl("^LOCUS", text))
    text <- base::split(text, split_indices)
  }
  return(text)
}
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

###################Function to extraxt protocluster data################
extract_protocluster_data <- function(GBFFfeat, refseq, explanation) {
  # Check if GBFFfeat is a list and extract the first element
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  # Convert GBFFfeat to a character vector
  GBFFfeat <- as.character(GBFFfeat)
  # Find the positions of "protocluster" and "proto_core" in the GBK feature
  protocluster_start <- grep("protocluster ", GBFFfeat)
  proto_core_start <- grep("proto_core ", GBFFfeat)
  
  if (length(protocluster_start) > 0) {
    # Initialize a data frame to store protocluster data
    protocluster_data <- data.frame(matrix(ncol = 7, nrow = length(protocluster_start)), stringsAsFactors = FALSE)
    # Set column names
    names(protocluster_data) <- c("refseq", "explanation", "protocluster_number", "FromTo", "Product","Category","core_location") 
    for (i in 1:length(protocluster_start)) {
      # Extract the text between "protocluster" and "proto_core"
      protocluster_text <- GBFFfeat[protocluster_start[i]:proto_core_start[i]]
      # Set values for refseq and explanation columns
      protocluster_data[i, "refseq"] <- refseq
      protocluster_data[i, "explanation"] <- explanation
      # Extract and clean protocluster number
      protocluster_number_line <- grep("/protocluster_number=", protocluster_text)
      protocluster_data[i, "protocluster_number"] <- as.numeric(gsub("[^0-9]", "", protocluster_text[protocluster_number_line]))
      # Extract FromTo value
      protocluster_data[i, "FromTo"] <- strsplit(protocluster_text, "\\s+")[[1]][3]
      # Extract and clean Product value
      protocluster_Type_line <- grep("/product=", protocluster_text)
      protocluster_data[i, "Product"] <- gsub("^\\s*/product=\"|\"$", "", protocluster_text[protocluster_Type_line])
      # Extract and clean Category value
      protocluster_Category_line<- grep("/category=", protocluster_text)
      protocluster_data[i, "Category"] <- gsub("^\\s*/category=\"|\"$", "", protocluster_text[protocluster_Category_line])
      # Extract and clean core_location value
      protocluster_core_line<- grep("/core_location=", protocluster_text)
      protocluster_data[i, "core_location"] <- gsub("[^0-9:]", "", protocluster_text[ protocluster_core_line])
    }
    return(protocluster_data)
  } else {
    # If no protocluster data found, return a data frame with NA values
    return(data.frame(refseq = NA, explanation = NA, protocluster_number = NA, FromTo = NA, Product = NA, Category = NA, core_location= NA))
  }
}
# Apply the extract_protocluster_data function to each element in GBFFfeat
result_list <- mapply(extract_protocluster_data, GBFFfeat, ref_data$refseq, ref_data$explanation, SIMPLIFY = FALSE)
# Combine the resulting data frames into one
combined_result_df <- do.call("rbind", lapply(result_list, function(x) cbind(x[[1]], x[, -1])))
# Rename the column names
names(combined_result_df) <- c("refseq", "explanation", "protocluster_number", "FromTo", "Product","Category","core_location")


write.csv(combined_result_df, file= "./2_results_protocluster.csv")

###################Function to extraxt cand_cluster data################
extract_cand_cluster_data <- function(GBFFfeat, refseq, explanation) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)
  cand_cluster_start<- grep("cand_cluster ", GBFFfeat)
  
  if (length(cand_cluster_start) > 0) {
    cand_cluster_data <- data.frame(matrix(ncol = 7, nrow = length(cand_cluster_start)), stringsAsFactors = FALSE)
    names(cand_cluster_data) <- c("refseq", "explanation","FromTo", "SMILES", "candidate_cluster_number", "product", "protocluster_number")
    
    for (i in 1:length(cand_cluster_start)) {
      cand_cluster_data[i, "refseq"] <- refseq
      cand_cluster_data[i, "explanation"] <- explanation
      
      # Determine the end of the current cand_cluster block
      end_line <- grep("^\\s+/tool=", GBFFfeat[cand_cluster_start[i]:(cand_cluster_start[i]+25)])
      if (length(end_line) == 0) {
        end_line <- cand_cluster_start[i] +25
      } else {
        end_line <- cand_cluster_start[i] + end_line[1] - 1
      }
      cand_cluster_text <- GBFFfeat[cand_cluster_start[i]:end_line]
      # Extract FromTo value
      cand_cluster_data[i, "FromTo"] <- strsplit(cand_cluster_text, "\\s+")[[1]][3]
      
      # Extract SMILES values that span multiple lines
      SMILES_values <- c()
      protocluster_SMILES_lines <- grep("/SMILES=", cand_cluster_text)
      if (length(protocluster_SMILES_lines) == 0) {
        cand_cluster_data[i, "SMILES"] <- NA
      } else {
        while (protocluster_SMILES_lines <= length(cand_cluster_text)) {
          SMILES_line <- gsub("\"", "", trimws(gsub("/SMILES=", "", cand_cluster_text[protocluster_SMILES_lines])))
          SMILES_values <- c(SMILES_values, SMILES_line)
          if (grepl("^\\s+/", cand_cluster_text[protocluster_SMILES_lines + 1])){
            break  # Stop if the next line starts with "/" and ends with a quote
          }
          protocluster_SMILES_lines <- protocluster_SMILES_lines + 1
        }
        cand_cluster_data[i, "SMILES"] <- paste(SMILES_values, collapse = "")
      }
      # Extract candidate_cluster_number values
      candidate_cluster_number_line <- grep("/candidate_cluster_number=", cand_cluster_text)
      cand_cluster_data_values <- as.numeric(gsub("[^0-9]", "", cand_cluster_text[candidate_cluster_number_line])) 
      cand_cluster_data[i, "candidate_cluster_number"]<-paste(cand_cluster_data_values, collapse = ", ")
      # Extract product values
      product_line <- grep("/product=", cand_cluster_text)
      cand_cluster_product_values<- gsub("^\\s*/product=\"|\"$", "", cand_cluster_text[product_line])
      cand_cluster_data[i, "product"] <- paste(cand_cluster_product_values, collapse = ", ")
      # Extract protocluster_number values
      cand_cluster_number_line <- grep("/protoclusters=", cand_cluster_text)
      cand_cluster_number_values<- as.numeric(gsub("[^0-9]", "", cand_cluster_text[cand_cluster_number_line]))
      cand_cluster_data[i, "protocluster_number"] <- paste(cand_cluster_number_values, collapse = ", ")
    }
    return(cand_cluster_data)
  } else {
    return(data.frame(refseq = NA, explanation = NA, FromTo = NA, SMILES= NA, candidate_cluster_number= NA, product = NA, protocluster_number = NA))
  }
}

# Apply the extract_cand_cluster_data function to each element in GBFFfeat
result_list_cand_cluster <- mapply(extract_cand_cluster_data, GBFFfeat, ref_data$refseq, ref_data$explanation, SIMPLIFY = FALSE)
# Combine the resulting data frames into one
combined_result_cand_cluster <- do.call("rbind", lapply(result_list_cand_cluster, function(x) cbind(x[[1]], x[, -1])))
# Rename the column names
names(combined_result_cand_cluster) <- c("refseq", "explanation","FromTo", "SMILES", "candidate_cluster_number", "product", "protocluster_number")


write.csv(combined_result_cand_cluster, file= "./3_results_cand_cluster.csv",row.names = FALSE)
####Selection
counts_C <- sapply(gregexpr("C", combined_result_cand_cluster$SMILES), function(x) sum(x >= 0))
counts_equals <- sapply(gregexpr("=", combined_result_cand_cluster$SMILES), function(x) sum(x >= 0))
selected_rows <- combined_result_cand_cluster[!grepl("N", combined_result_cand_cluster$SMILES) &
                                                counts_C >= 12 & counts_C <= 20 &
                                                counts_equals <= 4 &
                                                grepl("T1PKS|T2PKS|T3PKS|PKS-like", combined_result_cand_cluster$product), ]

write.csv(selected_rows, file= "./3_results_cand_cluster_Selected.csv",row.names = FALSE)

###################Function to extraxt PFam data################
extract_PFAM_data <- function(GBFFfeat, refseq, explanation) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)
  PFAM_start<- grep("PFAM_domain ", GBFFfeat)
  
  if (length(PFAM_start) > 0) {
    PFAM_data <- data.frame(matrix(ncol = 11, nrow = length(PFAM_start)), stringsAsFactors = FALSE)
    names(PFAM_data) <- c("refseq", "explanation","FromTo", "description", "evalue", "gene_ontologies", "label", "protein_start", "protein_end", "score", "translation")
    
    for (i in 1:length(PFAM_start)) {
      PFAM_data[i, "refseq"] <- refseq
      PFAM_data[i, "explanation"] <- explanation
      # Determine the end of the current cand_cluster block
      translation_line <- grep("^\\s+/translation=", GBFFfeat[PFAM_start[i]:(PFAM_start[i]+35)])
      if ((PFAM_start[i]+ translation_line[1]+20)>length(GBFFfeat)) {
        quote_search_line <- length(GBFFfeat)
      } else {
        quote_search_line <-(PFAM_start[i]+ translation_line[1]+30)
      }
      quote_counts <- sapply(gregexpr('"', GBFFfeat[(PFAM_start[i]+ translation_line[1]-1):quote_search_line]), function(x) sum(x >= 0))
      quote_counts[is.na(quote_counts)] <- 0
      end_line <- unlist(lapply(seq_along(quote_counts), function(i) rep(i, quote_counts[i])))
      PFAM_text <- GBFFfeat[PFAM_start[i]:(PFAM_start[i]+translation_line[1]+end_line[2]-2)]
      
      PFAM_data[i, "FromTo"] <- strsplit(PFAM_text, "\\s+")[[1]][3]
      # Extract and clean description value
      PFAM_description_line <- grep("/description=", PFAM_text)
      PFAM_data[i, "description"] <- paste(gsub("^\\s*/description=\"|\"$", "", PFAM_text[PFAM_description_line]), collapse = ", ")
      # Extract and clean evalue 
      PFAM_evalue_line <- grep("/evalue=", PFAM_text)
      PFAM_data[i, "evalue"] <- gsub("^\\s*/evalue=\"|\"$", "", PFAM_text[PFAM_evalue_line])[1]
      # Extract gene_ontologies values
      PFAM_ontologies_line <- grep("/gene_ontologies=", PFAM_text)
      PFAM_data[i, "gene_ontologies"]<- paste(gsub("^\\s*/gene_ontologies=\"|\"$", "", PFAM_text[PFAM_ontologies_line]) , collapse = ", ")
      # Extract and clean label value 
      PFAM_label_line <- grep("/label=", PFAM_text)
      PFAM_data[i, "label"] <- gsub("^\\s*/label=\"|\"$", "", PFAM_text[PFAM_label_line])[1]
      # Extract and clean protein_end value
      PFAM_protein_end_line <- grep("/protein_end=", PFAM_text)
      PFAM_data[i, "protein_end"] <- gsub("^\\s*/protein_end=\"|\"$", "", PFAM_text[PFAM_protein_end_line])[1]
      # Extract and clean protein_start value 
      PFAM_protein_start_line <- grep("/protein_start=", PFAM_text)
      PFAM_data[i, "protein_start"] <- as.numeric(gsub("^\\s*/protein_start=\"|\"$", "", PFAM_text[PFAM_protein_start_line]))[1]
      # Extract and clean score value 
      PFAM_score_line <- grep("/score=", PFAM_text)
      PFAM_data[i, "score"] <- gsub("^\\s*/score=\"|\"$", "", PFAM_text[PFAM_score_line])[1]
      # Extract and clean translation value 
      translation_text <- PFAM_text[translation_line[1]:length(PFAM_text)]
      PFAM_data[i, "translation"] <- paste(gsub("\\s+", "",gsub("^\\s*/translation=\"|\"$", "", translation_text)), collapse = "")
    }
    return(PFAM_data)
  } else {
    return(data.frame(refseq = NA, explanation = NA, FromTo = NA, description= NA, evalue = NA, gene_ontologies = NA, label= NA, protein_start= NA, protein_end= NA, score= NA, translation= NA))
  }
}
# Apply the extract_PFAM_data function to each element in GBFFfeat
result_list_PFAM_data <- mapply(extract_PFAM_data, GBFFfeat, ref_data$refseq, ref_data$explanation, SIMPLIFY = FALSE)

# Combine the resulting data frames into one
combined_result_PFAM_data <- do.call("rbind", lapply(result_list_PFAM_data, function(x) cbind(x[[1]], x[, -1])))
# Rename the column names
names(combined_result_PFAM_data) <- c("refseq", "explanation","FromTo",  "description", "evalue", "gene_ontologies", "label", "protein_start", "protein_end", "score", "translation")

write.csv(combined_result_PFAM_data, file= "./4_results_PFAM_data.csv",row.names = FALSE)

###################Function to extraxt aSDomain data################
extract_aSDomain_data <- function(GBFFfeat, refseq, explanation) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)
  aSDomain_start<- grep("aSDomain ", GBFFfeat)
  
  if (length(aSDomain_start) > 0) {
    aSDomain_data <- data.frame(matrix(ncol = 12, nrow = length(aSDomain_start)), stringsAsFactors = FALSE)
    names(aSDomain_data) <- c("refseq", "explanation", "FromTo", "description", "evalue","identifier", "label", "protein_start", "protein_end", "score","specificity", "translation")
    
    for (i in 1:length(aSDomain_start)) {
      aSDomain_data[i, "refseq"] <- refseq
      aSDomain_data[i, "explanation"] <- explanation
      # Determine the end of the current cand_cluster block
      translation_line <- grep("^\\s+/translation=", GBFFfeat[aSDomain_start[i]:(aSDomain_start[i]+35)])
      if ((aSDomain_start[i]+ translation_line[1]+20)>length(GBFFfeat)){
        quote_search_line <- length(GBFFfeat)
      } else {
        quote_search_line <-(aSDomain_start[i]+ translation_line[1]+30)
      }
      quote_counts <- sapply(gregexpr('"', GBFFfeat[(aSDomain_start[i]+ translation_line[1]-1):quote_search_line]), function(x) sum(x >= 0))
      quote_counts[is.na(quote_counts)] <- 0
      end_line <- unlist(lapply(seq_along(quote_counts), function(i) rep(i, quote_counts[i])))
      aSDomain_text <- GBFFfeat[aSDomain_start[i]:(aSDomain_start[i]+translation_line[1]+end_line[2]-2)]
      
      aSDomain_data[i, "FromTo"] <- strsplit(aSDomain_text, "\\s+")[[1]][3]
      # Extract and clean description value
      aSDomain_description_line <- grep("/description=", aSDomain_text)
      aSDomain_data[i, "description"] <- paste(gsub("^\\s*/description=\"|\"$", "", aSDomain_text[aSDomain_description_line]), collapse = ", ")
      # Extract and clean evalue
      aSDomain_evalue_line <- grep("/evalue=", aSDomain_text)
      aSDomain_data[i, "evalue"] <- gsub("^\\s*/evalue=\"|\"$", "", aSDomain_text[aSDomain_evalue_line])[1]
      # Extract and clean identifier
      aSDomain_identifier_line <- grep("/identifier=", aSDomain_text)
      aSDomain_data[i, "identifier"] <- gsub("^\\s*/identifier=\"|\"$", "", aSDomain_text[aSDomain_identifier_line])[1]
      # Extract and clean label value
      aSDomain_label_line <- grep("/label=", aSDomain_text)
      aSDomain_data[i, "label"] <- gsub("^\\s*/label=\"|\"$", "", aSDomain_text[aSDomain_label_line])[1]
      # Extract and clean protein_end value
      aSDomain_protein_end_line <- grep("/protein_end=", aSDomain_text)
      aSDomain_data[i, "protein_end"] <- gsub("^\\s*/protein_end=\"|\"$", "", aSDomain_text[aSDomain_protein_end_line])[1]
      # Extract and clean protein_start value
      aSDomain_protein_start_line <- grep("/protein_start=", aSDomain_text)
      aSDomain_data[i, "protein_start"] <- as.numeric(gsub("^\\s*/protein_start=\"|\"$", "", aSDomain_text[aSDomain_protein_start_line]))[1]
      # Extract and clean score value
      aSDomain_score_line <- grep("/score=", aSDomain_text)
      aSDomain_data[i, "score"] <- gsub("^\\s*/score=\"|\"$", "", aSDomain_text[aSDomain_score_line])[1]
      # Extract specificity values
      aSDomain_specificity_line <- grep("/specificity=", aSDomain_text)
      aSDomain_data[i, "specificity"] <- paste(gsub("^\\s*/specificity=\"|\"$", "", aSDomain_text[aSDomain_specificity_line]) , collapse = ", ")
      # Extract and clean translation value
      translation_text <- aSDomain_text[translation_line[1]:length(aSDomain_text)]
      aSDomain_data[i, "translation"] <- paste(gsub("\\s+", "", gsub("^\\s*/translation=\"|\"$", "", translation_text)), collapse = "")
    }
    return(aSDomain_data)
  } else {
    return(data.frame(refseq = NA, explanation = NA, FromTo = NA, description = NA, evalue = NA,identifier = NA, label = NA, protein_start = NA, protein_end = NA, score = NA, specificity = NA, translation = NA))
  }
}

# Apply the extract_aSDomain_data function to each element in GBFFfeat
result_list_aSDomain_data <- mapply(extract_aSDomain_data, GBFFfeat, ref_data$refseq, ref_data$explanation, SIMPLIFY = FALSE)
# Combine the resulting data frames into one
combined_result_aSDomain_data <- do.call("rbind", lapply(result_list_aSDomain_data, function(x) cbind(x[[1]], x[, -1])))
# Rename the column names
names(combined_result_aSDomain_data) <- c("refseq", "explanation", "FromTo", "description", "evalue", "identifier" , "label", "protein_start", "protein_end", "score","specificity", "translation")


write.csv(combined_result_aSDomain_data, file = "./5_results_aSDomain_data.csv", row.names = FALSE)

#######################text file########################
# Get a list of file paths to all .txt files in the "knownclusterblast" directory
file_list <- list.files(path = "./knownclusterblast/",pattern = "\\.txt$",full.names = TRUE)
# Split the file names by underscores to extract refseq information
refseq <- strsplit(list.files(path = "./knownclusterblast/",pattern = "\\.txt$"), "_")
# Read the contents of each .txt file into a lis
file_contents <- lapply(file_list, readLines)
# Extract "Compounds" information from each file
Compounds <- lapply(file_contents, function(x) {
  # Find the line containing "Details:"
  details_line_index <- grep("Details:", x, fixed = TRUE)
  # Determine the index of lines containing "Source:"
  source_lines_index <- details_line_index + 1:4
  source_lines <- x[source_lines_index]
  # Extract and clean the "Source" information
  result <- grep("Source:", source_lines, fixed = TRUE, value = TRUE)
  ifelse(length(result) > 0, sub("\\s", "", sub("Source:", "", result)), NA)
})
# Create a data frame to store the extracted Compounds data
Compounds_data <- data.frame(matrix(ncol = 3, nrow = length(file_list)), stringsAsFactors = FALSE)
names(Compounds_data) <- c("refseq", "protocluster_number", "Compounds")
# Fill in the data frame with the extracted information
for (i in 1:length(file_list)) {
  Compounds_data[i, "refseq"] <- paste(refseq[[i]][-3], collapse = "_")
  Compounds_data[i, "protocluster_number"] <- gsub("^c|\\.txt$", "", refseq[[i]][3])
  Compounds_data[i, "Compounds"] <- Compounds[[i]]}

##########################Create a result file in each merged folder
# Merge the combined_result_df data with Compounds_data based on refseq and protocluster_number
merged_data <- merge(combined_result_df, Compounds_data, by = c("refseq", "protocluster_number"), sort = FALSE) 
# Write the merged data to a CSV file
write.csv(merged_data, file= "./2_mergedResults.csv")
cat("End of GenomeMining.R script\n")
69
82
