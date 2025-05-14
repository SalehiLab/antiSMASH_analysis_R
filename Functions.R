library(dplyr)
library(rentrez)
###################function to find taxonomy################ 
get_taxonomy <- function(refseq_ids) {
  result <- strsplit(entrez_fetch(db = "nuccore", id = refseq_ids, rettype = "gb", retmode = "text"), "\n")[[1]]
  taxonomy_lines <- grep("^  ORGANISM", result)
  REFERENCE_lines <- grep("REFERENCE", result)
  taxonomy <- gsub("\\s", "", paste0(result[(taxonomy_lines+1 ):(REFERENCE_lines[1]-1)],collapse = ""))
  return(c(refseq_ids,taxonomy))
}
###################function to read a GenBank file################
readGBFF <- function (file, text = readLines(file), verbose = FALSE){
  if (is(text, "character")) {
    split_indices <- cumsum(grepl("^LOCUS", text))
    text <- base::split(text, split_indices)
  }
  return(text)
}
###################Function to extract protocluster data################
extract_protocluster_data <- function(GBFFfeat, refseq, Deffinition) {
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
    protocluster_data <- data.frame(matrix(ncol = 14, nrow = length(protocluster_start)), stringsAsFactors = FALSE)
    # Set column names
    names(protocluster_data) <- c("refseq","Deffinition", "Length","Protocluster_number","FromTo","Product",
                                  "Category","Core_location","aStool","Contig_edge","Cutoff",
                                  "Tool","Neighbourhood","Detection_rule") 
    for (i in 1:length(protocluster_start)) {
      # Extract the text between "protocluster" and "proto_core"
      protocluster_text <- GBFFfeat[protocluster_start[i]:proto_core_start[i]]
      # Set values for refseq and Deffinition columns
      protocluster_data[i, "refseq"] <- refseq
      protocluster_data[i, "Deffinition"] <- Deffinition
      # Extract and clean length
      protocluster_data[i,"Length"] <- as.numeric(strsplit(grep("LOCUS", GBFFfeat,value = TRUE), "\\s+")[[1]][3])
      # Extract and clean protocluster number
      protocluster_number_line <- grep("/protocluster_number=", protocluster_text)
      protocluster_data[i, "Protocluster_number"] <- as.numeric(gsub("[^0-9]", "", protocluster_text[protocluster_number_line]))
      # Extract FromTo value
      protocluster_data[i, "FromTo"] <- strsplit(protocluster_text, "\\s+")[[1]][3]
      # Extract and clean Product value
      protocluster_Type_line <- grep("/product=", protocluster_text)
      protocluster_data[i, "Product"] <- gsub("^\\s*/product=\"|\"$", "", protocluster_text[protocluster_Type_line])
      # Extract and clean Category value
      protocluster_Category_line<- grep("/category=", protocluster_text)
      protocluster_data[i, "Category"] <- gsub("^\\s*/category=\"|\"$", "", protocluster_text[protocluster_Category_line])
      # Extract and clean core_location value
      protocluster_core_line <- grep("/core_location=", protocluster_text)
      protocluster_data[i, "Core_location"] <- gsub("[^0-9:]", "", protocluster_text[ protocluster_core_line])
      # Extract and clean contig_edge value
      protocluster_contig_edge_line<- grep("/contig_edge=", protocluster_text)
      protocluster_data[i, "Contig_edge"] <- gsub("^\\s*/contig_edge=\"|\"$", "", protocluster_text[protocluster_contig_edge_line])
      # Extract and clean aStool value
      protocluster_aStool_line <- grep("/aStool=", protocluster_text)
      protocluster_data[i, "aStool"] <- gsub("^\\s*/aStool=\"|\"$", "", protocluster_text[protocluster_aStool_line])
      # Extract and clean cutoff value
      protocluster_cutoff_line <- grep("/cutoff=", protocluster_text)
      protocluster_data[i, "Cutoff"] <- as.numeric(gsub("[^0-9:]", "", protocluster_text[protocluster_cutoff_line]))
      # Extract and clean tool value
      protocluster_tool_line<- grep("/tool=", protocluster_text)
      protocluster_data[i, "Tool"] <- gsub("^\\s*/tool=\"|\"$", "", protocluster_text[protocluster_tool_line])
      # Extract and clean neighbourhood value
      protocluster_neighbourhood_line<- grep("/neighbourhood=", protocluster_text)
      protocluster_data[i, "Neighbourhood"] <- as.numeric(gsub("[^0-9:]", "", protocluster_text[protocluster_neighbourhood_line]))
      # Extract and clean detection_rule value
      protocluster_detection_rule_line<- grep("/detection_rule=", protocluster_text)
      # Find the start of the detection_rule
      start_line <- grep("/detection_rule=", protocluster_text)
      end_line <- start_line
      while (!grepl("\"$", protocluster_text[end_line])) {
        end_line <- end_line + 1
      }
      detection_rule_raw <- paste(protocluster_text[start_line:end_line], collapse = " ")
      protocluster_data[i, "Detection_rule"] <- gsub("^\\s*/detection_rule=\"|\"$", "", detection_rule_raw)
    }
    return(protocluster_data)
  } else {
    # If no protocluster data found, return a data frame with NA values
    return(data.frame(refseq = NA, Deffinition = NA, Length= NA, Protocluster_number = NA, FromTo = NA, 
                      Product = NA, Category = NA, Core_location= NA, aStool = NA, Contig_edge = NA, Cutoff = NA,
                      Tool=NA ,Neighbourhood = NA, Detection_rule = NA, stringsAsFactors = FALSE))
  }
}

###################Function to extract proto_core data################
extract_proto_core_data <- function(GBFFfeat, refseq, Deffinition) {
  # Check if GBFFfeat is a list and extract the first element
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  # Convert GBFFfeat to a character vector
  GBFFfeat <- as.character(GBFFfeat)
  # Find the positions of "proto_core" 
  proto_core_start <- grep("proto_core ", GBFFfeat)
  
  
  if (length(proto_core_start) > 0) {
    # Initialize a data frame to store proto_core data
    proto_core_data <- data.frame(matrix(ncol = 11, nrow = length(proto_core_start)), stringsAsFactors = FALSE)
    # Set column names
    names(proto_core_data) <- c("refseq","Deffinition", "Length","Protocluster_number","FromTo","Product",
                                "aStool","Cutoff","Tool","Neighbourhood","Detection_rule")
    
    for (i in 1:length(proto_core_start)) {
      # Set values for refseq and Deffinition columns
      proto_core_data[i, "refseq"] <- refseq
      proto_core_data[i, "Deffinition"] <- Deffinition
      
      # Determine the end of the current cand_cluster block
      end_line <- grep("^\\s+/protocluster_number=", GBFFfeat[proto_core_start[i]:(proto_core_start[i]+25)])
      if (length(end_line) == 0) {
        end_line <- proto_core_start[i] +25
      } else {
        end_line <- proto_core_start[i] + end_line[1] - 1
      }
      proto_core_text <- GBFFfeat[proto_core_start[i]:end_line]
      
      # Extract and clean length
      proto_core_data[i,"Length"] <- as.numeric(strsplit(grep("LOCUS", GBFFfeat,value = TRUE), "\\s+")[[1]][3])
      # Extract and clean protocluster number
      proto_core_number_line <- grep("/protocluster_number=", proto_core_text)
      proto_core_data[i, "Protocluster_number"] <- as.numeric(gsub("[^0-9]", "", proto_core_text[proto_core_number_line]))
      # Extract FromTo value
      proto_core_data[i, "FromTo"] <- strsplit(proto_core_text, "\\s+")[[1]][3]
      # Extract and clean Product value
      proto_core_Type_line <- grep("/product=", proto_core_text)
      proto_core_data[i, "Product"] <- paste(gsub("^\\s*/product=\"|\"$", "", proto_core_text[proto_core_Type_line]), collapse = "; ")
      # Extract and clean aStool value
      proto_core_aStool_line <- grep("/aStool=", proto_core_text)
      proto_core_data[i, "aStool"] <- gsub("^\\s*/aStool=\"|\"$", "", proto_core_text[proto_core_aStool_line])
      # Extract and clean cutoff value
      proto_core_cutoff_line <- grep("/cutoff=", proto_core_text)
      proto_core_data[i, "Cutoff"] <- as.numeric(gsub("[^0-9:]", "", proto_core_text[proto_core_cutoff_line]))
      # Extract and clean tool value
      proto_core_tool_line<- grep("/tool=", proto_core_text)
      proto_core_data[i, "Tool"] <- gsub("^\\s*/tool=\"|\"$", "", proto_core_text[proto_core_tool_line])
      # Extract and clean neighbourhood value
      proto_core_neighbourhood_line<- grep("/neighbourhood=", proto_core_text)
      proto_core_data[i, "Neighbourhood"] <- as.numeric(gsub("[^0-9:]", "", proto_core_text[proto_core_neighbourhood_line]))
      # Extract and clean detection_rule value
      proto_core_detection_rule_line<- grep("/detection_rule=", proto_core_text)
      # Find the start of the detection_rule
      start_line <- grep("/detection_rule=", proto_core_text)
      end_line <- start_line
      while (!grepl("\"$", proto_core_text[end_line])) {
        end_line <- end_line + 1
      }
      protocore_detection_rule_raw <- paste(proto_core_text[start_line:end_line], collapse = " ")
      proto_core_data[i, "Detection_rule"] <- gsub("^\\s*/detection_rule=\"|\"$", "",protocore_detection_rule_raw)
    }
    return(proto_core_data)
  } else {
    # If no protocluster data found, return a data frame with NA values
    return(data.frame(refseq = NA, Deffinition = NA, Length= NA, Protocluster_number = NA, FromTo = NA, 
                      Product = NA, aStool = NA, Cutoff = NA,
                      Tool=NA ,Neighbourhood = NA, Detection_rule = NA, stringsAsFactors = FALSE))
  }
}


