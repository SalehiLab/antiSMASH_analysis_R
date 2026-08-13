# extract_features.R
# Parsers for antiSMASH GenBank feature blocks: protocluster, proto_core,
# cand_cluster, PFAM_domain, and aSDomain. Each function takes a single
# genome's GBFFfeat (as returned by load_merged_gbff()) plus its refseq /
# Deffinition, and returns one row per feature found.
#
# These are the ONE canonical definition of each function -- previously
# extract_protocluster_data() existed in two incompatible versions (14
# columns here vs. 7 columns in GenomeMining.R). This file keeps the fuller
# 14-column version; run_genome_mining.R has been updated accordingly.

#' Extract "protocluster" feature data from a single genome's GBK feature
#' block.
#'
#' @param GBFFfeat A GenBank record (or list containing one).
#' @param refseq Accession/identifier to tag output rows with.
#' @param Deffinition Definition/description string to tag output rows with.
#' @return A data frame with one row per protocluster.
extract_protocluster_data <- function(GBFFfeat, refseq, Deffinition) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)

  protocluster_start <- grep("protocluster ", GBFFfeat)
  proto_core_start <- grep("proto_core ", GBFFfeat)

  if (length(protocluster_start) == 0) {
    return(data.frame(
      refseq = NA, Deffinition = NA, Length = NA, Protocluster_number = NA, FromTo = NA,
      Product = NA, Category = NA, Core_location = NA, aStool = NA, Contig_edge = NA, Cutoff = NA,
      Tool = NA, Neighbourhood = NA, Detection_rule = NA, stringsAsFactors = FALSE
    ))
  }

  protocluster_data <- data.frame(matrix(ncol = 14, nrow = length(protocluster_start)), stringsAsFactors = FALSE)
  names(protocluster_data) <- c(
    "refseq", "Deffinition", "Length", "Protocluster_number", "FromTo", "Product",
    "Category", "Core_location", "aStool", "Contig_edge", "Cutoff",
    "Tool", "Neighbourhood", "Detection_rule"
  )

  for (i in seq_along(protocluster_start)) {
    protocluster_text <- GBFFfeat[protocluster_start[i]:proto_core_start[i]]

    protocluster_data[i, "refseq"] <- refseq
    protocluster_data[i, "Deffinition"] <- Deffinition
    protocluster_data[i, "Length"] <- as.numeric(strsplit(grep("LOCUS", GBFFfeat, value = TRUE), "\\s+")[[1]][3])

    protocluster_number_line <- grep("/protocluster_number=", protocluster_text)
    protocluster_data[i, "Protocluster_number"] <- as.numeric(gsub("[^0-9]", "", protocluster_text[protocluster_number_line]))

    protocluster_data[i, "FromTo"] <- strsplit(protocluster_text, "\\s+")[[1]][3]

    protocluster_Type_line <- grep("/product=", protocluster_text)
    protocluster_data[i, "Product"] <- gsub("^\\s*/product=\"|\"$", "", protocluster_text[protocluster_Type_line])

    protocluster_Category_line <- grep("/category=", protocluster_text)
    protocluster_data[i, "Category"] <- gsub("^\\s*/category=\"|\"$", "", protocluster_text[protocluster_Category_line])

    protocluster_core_line <- grep("/core_location=", protocluster_text)
    protocluster_data[i, "Core_location"] <- gsub("[^0-9:]", "", protocluster_text[protocluster_core_line])

    protocluster_contig_edge_line <- grep("/contig_edge=", protocluster_text)
    protocluster_data[i, "Contig_edge"] <- gsub("^\\s*/contig_edge=\"|\"$", "", protocluster_text[protocluster_contig_edge_line])

    protocluster_aStool_line <- grep("/aStool=", protocluster_text)
    protocluster_data[i, "aStool"] <- gsub("^\\s*/aStool=\"|\"$", "", protocluster_text[protocluster_aStool_line])

    protocluster_cutoff_line <- grep("/cutoff=", protocluster_text)
    protocluster_data[i, "Cutoff"] <- as.numeric(gsub("[^0-9:]", "", protocluster_text[protocluster_cutoff_line]))

    protocluster_tool_line <- grep("/tool=", protocluster_text)
    protocluster_data[i, "Tool"] <- gsub("^\\s*/tool=\"|\"$", "", protocluster_text[protocluster_tool_line])

    protocluster_neighbourhood_line <- grep("/neighbourhood=", protocluster_text)
    protocluster_data[i, "Neighbourhood"] <- as.numeric(gsub("[^0-9:]", "", protocluster_text[protocluster_neighbourhood_line]))

    start_line <- grep("/detection_rule=", protocluster_text)
    end_line <- start_line
    while (!grepl("\"$", protocluster_text[end_line])) {
      end_line <- end_line + 1
    }
    detection_rule_raw <- paste(protocluster_text[start_line:end_line], collapse = " ")
    protocluster_data[i, "Detection_rule"] <- gsub("^\\s*/detection_rule=\"|\"$", "", detection_rule_raw)
  }

  protocluster_data
}

#' Extract "proto_core" feature data (the core biosynthetic region within
#' each protocluster) from a single genome's GBK feature block.
#'
#' @inheritParams extract_protocluster_data
#' @return A data frame with one row per proto_core.
extract_proto_core_data <- function(GBFFfeat, refseq, Deffinition) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)

  proto_core_start <- grep("proto_core ", GBFFfeat)

  if (length(proto_core_start) == 0) {
    return(data.frame(
      refseq = NA, Deffinition = NA, Length = NA, Protocluster_number = NA, FromTo = NA,
      Product = NA, aStool = NA, Cutoff = NA,
      Tool = NA, Neighbourhood = NA, Detection_rule = NA, stringsAsFactors = FALSE
    ))
  }

  proto_core_data <- data.frame(matrix(ncol = 11, nrow = length(proto_core_start)), stringsAsFactors = FALSE)
  names(proto_core_data) <- c(
    "refseq", "Deffinition", "Length", "Protocluster_number", "FromTo", "Product",
    "aStool", "Cutoff", "Tool", "Neighbourhood", "Detection_rule"
  )

  for (i in seq_along(proto_core_start)) {
    proto_core_data[i, "refseq"] <- refseq
    proto_core_data[i, "Deffinition"] <- Deffinition

    end_line <- grep("^\\s+/protocluster_number=", GBFFfeat[proto_core_start[i]:(proto_core_start[i] + 25)])
    if (length(end_line) == 0) {
      end_line <- proto_core_start[i] + 25
    } else {
      end_line <- proto_core_start[i] + end_line[1] - 1
    }
    proto_core_text <- GBFFfeat[proto_core_start[i]:end_line]

    proto_core_data[i, "Length"] <- as.numeric(strsplit(grep("LOCUS", GBFFfeat, value = TRUE), "\\s+")[[1]][3])

    proto_core_number_line <- grep("/protocluster_number=", proto_core_text)
    proto_core_data[i, "Protocluster_number"] <- as.numeric(gsub("[^0-9]", "", proto_core_text[proto_core_number_line]))

    proto_core_data[i, "FromTo"] <- strsplit(proto_core_text, "\\s+")[[1]][3]

    proto_core_Type_line <- grep("/product=", proto_core_text)
    proto_core_data[i, "Product"] <- paste(gsub("^\\s*/product=\"|\"$", "", proto_core_text[proto_core_Type_line]), collapse = "; ")

    proto_core_aStool_line <- grep("/aStool=", proto_core_text)
    proto_core_data[i, "aStool"] <- gsub("^\\s*/aStool=\"|\"$", "", proto_core_text[proto_core_aStool_line])

    proto_core_cutoff_line <- grep("/cutoff=", proto_core_text)
    proto_core_data[i, "Cutoff"] <- as.numeric(gsub("[^0-9:]", "", proto_core_text[proto_core_cutoff_line]))

    proto_core_tool_line <- grep("/tool=", proto_core_text)
    proto_core_data[i, "Tool"] <- gsub("^\\s*/tool=\"|\"$", "", proto_core_text[proto_core_tool_line])

    proto_core_neighbourhood_line <- grep("/neighbourhood=", proto_core_text)
    proto_core_data[i, "Neighbourhood"] <- as.numeric(gsub("[^0-9:]", "", proto_core_text[proto_core_neighbourhood_line]))

    start_line <- grep("/detection_rule=", proto_core_text)
    end_line2 <- start_line
    while (!grepl("\"$", proto_core_text[end_line2])) {
      end_line2 <- end_line2 + 1
    }
    protocore_detection_rule_raw <- paste(proto_core_text[start_line:end_line2], collapse = " ")
    proto_core_data[i, "Detection_rule"] <- gsub("^\\s*/detection_rule=\"|\"$", "", protocore_detection_rule_raw)
  }

  proto_core_data
}

#' Extract "cand_cluster" feature data (candidate cluster boundaries,
#' predicted SMILES, and linked protocluster numbers) from a single genome's
#' GBK feature block.
#'
#' @inheritParams extract_protocluster_data
#' @return A data frame with one row per candidate cluster.
extract_cand_cluster_data <- function(GBFFfeat, refseq, Deffinition) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)

  cand_cluster_start <- grep("cand_cluster ", GBFFfeat)

  if (length(cand_cluster_start) == 0) {
    return(data.frame(
      refseq = NA, Deffinition = NA, FromTo = NA, SMILES = NA,
      candidate_cluster_number = NA, product = NA, protocluster_number = NA
    ))
  }

  cand_cluster_data <- data.frame(matrix(ncol = 7, nrow = length(cand_cluster_start)), stringsAsFactors = FALSE)
  names(cand_cluster_data) <- c(
    "refseq", "Deffinition", "FromTo", "SMILES", "candidate_cluster_number", "product", "protocluster_number"
  )

  for (i in seq_along(cand_cluster_start)) {
    cand_cluster_data[i, "refseq"] <- refseq
    cand_cluster_data[i, "Deffinition"] <- Deffinition

    end_line <- grep("^\\s+/tool=", GBFFfeat[cand_cluster_start[i]:(cand_cluster_start[i] + 25)])
    if (length(end_line) == 0) {
      end_line <- cand_cluster_start[i] + 25
    } else {
      end_line <- cand_cluster_start[i] + end_line[1] - 1
    }
    cand_cluster_text <- GBFFfeat[cand_cluster_start[i]:end_line]

    cand_cluster_data[i, "FromTo"] <- strsplit(cand_cluster_text, "\\s+")[[1]][3]

    SMILES_values <- c()
    protocluster_SMILES_lines <- grep("/SMILES=", cand_cluster_text)
    if (length(protocluster_SMILES_lines) == 0) {
      cand_cluster_data[i, "SMILES"] <- NA
    } else {
      while (protocluster_SMILES_lines <= length(cand_cluster_text)) {
        SMILES_line <- gsub("\"", "", trimws(gsub("/SMILES=", "", cand_cluster_text[protocluster_SMILES_lines])))
        SMILES_values <- c(SMILES_values, SMILES_line)
        if (grepl("^\\s+/", cand_cluster_text[protocluster_SMILES_lines + 1])) {
          break
        }
        protocluster_SMILES_lines <- protocluster_SMILES_lines + 1
      }
      cand_cluster_data[i, "SMILES"] <- paste(SMILES_values, collapse = "")
    }

    candidate_cluster_number_line <- grep("/candidate_cluster_number=", cand_cluster_text)
    cand_cluster_data_values <- as.numeric(gsub("[^0-9]", "", cand_cluster_text[candidate_cluster_number_line]))
    cand_cluster_data[i, "candidate_cluster_number"] <- paste(cand_cluster_data_values, collapse = ", ")

    product_line <- grep("/product=", cand_cluster_text)
    cand_cluster_product_values <- gsub("^\\s*/product=\"|\"$", "", cand_cluster_text[product_line])
    cand_cluster_data[i, "product"] <- paste(cand_cluster_product_values, collapse = ", ")

    cand_cluster_number_line <- grep("/protoclusters=", cand_cluster_text)
    cand_cluster_number_values <- as.numeric(gsub("[^0-9]", "", cand_cluster_text[cand_cluster_number_line]))
    cand_cluster_data[i, "protocluster_number"] <- paste(cand_cluster_number_values, collapse = ", ")
  }

  cand_cluster_data
}

#' Extract "PFAM_domain" feature data from a single genome's GBK feature
#' block.
#'
#' @inheritParams extract_protocluster_data
#' @return A data frame with one row per PFAM domain hit.
extract_PFAM_data <- function(GBFFfeat, refseq, Deffinition) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)

  PFAM_start <- grep("PFAM_domain ", GBFFfeat)

  if (length(PFAM_start) == 0) {
    return(data.frame(
      refseq = NA, Deffinition = NA, FromTo = NA, description = NA, evalue = NA,
      gene_ontologies = NA, label = NA, protein_start = NA, protein_end = NA, score = NA, translation = NA
    ))
  }

  PFAM_data <- data.frame(matrix(ncol = 11, nrow = length(PFAM_start)), stringsAsFactors = FALSE)
  names(PFAM_data) <- c(
    "refseq", "Deffinition", "FromTo", "description", "evalue", "gene_ontologies",
    "label", "protein_start", "protein_end", "score", "translation"
  )

  for (i in seq_along(PFAM_start)) {
    PFAM_data[i, "refseq"] <- refseq
    PFAM_data[i, "Deffinition"] <- Deffinition

    translation_line <- grep("^\\s+/translation=", GBFFfeat[PFAM_start[i]:(PFAM_start[i] + 35)])
    if ((PFAM_start[i] + translation_line[1] + 20) > length(GBFFfeat)) {
      quote_search_line <- length(GBFFfeat)
    } else {
      quote_search_line <- (PFAM_start[i] + translation_line[1] + 30)
    }
    quote_counts <- sapply(gregexpr('"', GBFFfeat[(PFAM_start[i] + translation_line[1] - 1):quote_search_line]), function(x) sum(x >= 0))
    quote_counts[is.na(quote_counts)] <- 0
    end_line <- unlist(lapply(seq_along(quote_counts), function(i) rep(i, quote_counts[i])))
    PFAM_text <- GBFFfeat[PFAM_start[i]:(PFAM_start[i] + translation_line[1] + end_line[2] - 2)]

    PFAM_data[i, "FromTo"] <- strsplit(PFAM_text, "\\s+")[[1]][3]

    PFAM_description_line <- grep("/description=", PFAM_text)
    PFAM_data[i, "description"] <- paste(gsub("^\\s*/description=\"|\"$", "", PFAM_text[PFAM_description_line]), collapse = ", ")

    PFAM_evalue_line <- grep("/evalue=", PFAM_text)
    PFAM_data[i, "evalue"] <- gsub("^\\s*/evalue=\"|\"$", "", PFAM_text[PFAM_evalue_line])[1]

    PFAM_ontologies_line <- grep("/gene_ontologies=", PFAM_text)
    PFAM_data[i, "gene_ontologies"] <- paste(gsub("^\\s*/gene_ontologies=\"|\"$", "", PFAM_text[PFAM_ontologies_line]), collapse = ", ")

    PFAM_label_line <- grep("/label=", PFAM_text)
    PFAM_data[i, "label"] <- gsub("^\\s*/label=\"|\"$", "", PFAM_text[PFAM_label_line])[1]

    PFAM_protein_end_line <- grep("/protein_end=", PFAM_text)
    PFAM_data[i, "protein_end"] <- gsub("^\\s*/protein_end=\"|\"$", "", PFAM_text[PFAM_protein_end_line])[1]

    PFAM_protein_start_line <- grep("/protein_start=", PFAM_text)
    PFAM_data[i, "protein_start"] <- as.numeric(gsub("^\\s*/protein_start=\"|\"$", "", PFAM_text[PFAM_protein_start_line]))[1]

    PFAM_score_line <- grep("/score=", PFAM_text)
    PFAM_data[i, "score"] <- gsub("^\\s*/score=\"|\"$", "", PFAM_text[PFAM_score_line])[1]

    translation_text <- PFAM_text[translation_line[1]:length(PFAM_text)]
    PFAM_data[i, "translation"] <- paste(gsub("\\s+", "", gsub("^\\s*/translation=\"|\"$", "", translation_text)), collapse = "")
  }

  PFAM_data
}

#' Extract "aSDomain" feature data from a single genome's GBK feature block.
#'
#' @inheritParams extract_protocluster_data
#' @return A data frame with one row per antiSMASH domain hit.
extract_aSDomain_data <- function(GBFFfeat, refseq, Deffinition) {
  if (is.list(GBFFfeat)) {
    GBFFfeat <- GBFFfeat[[1]]
  }
  GBFFfeat <- as.character(GBFFfeat)

  aSDomain_start <- grep("aSDomain ", GBFFfeat)

  if (length(aSDomain_start) == 0) {
    return(data.frame(
      refseq = NA, Deffinition = NA, FromTo = NA, description = NA, evalue = NA, identifier = NA,
      label = NA, protein_start = NA, protein_end = NA, score = NA, specificity = NA, translation = NA
    ))
  }

  aSDomain_data <- data.frame(matrix(ncol = 12, nrow = length(aSDomain_start)), stringsAsFactors = FALSE)
  names(aSDomain_data) <- c(
    "refseq", "Deffinition", "FromTo", "description", "evalue", "identifier",
    "label", "protein_start", "protein_end", "score", "specificity", "translation"
  )

  for (i in seq_along(aSDomain_start)) {
    aSDomain_data[i, "refseq"] <- refseq
    aSDomain_data[i, "Deffinition"] <- Deffinition

    translation_line <- grep("^\\s+/translation=", GBFFfeat[aSDomain_start[i]:(aSDomain_start[i] + 35)])
    if ((aSDomain_start[i] + translation_line[1] + 20) > length(GBFFfeat)) {
      quote_search_line <- length(GBFFfeat)
    } else {
      quote_search_line <- (aSDomain_start[i] + translation_line[1] + 30)
    }
    quote_counts <- sapply(gregexpr('"', GBFFfeat[(aSDomain_start[i] + translation_line[1] - 1):quote_search_line]), function(x) sum(x >= 0))
    quote_counts[is.na(quote_counts)] <- 0
    end_line <- unlist(lapply(seq_along(quote_counts), function(i) rep(i, quote_counts[i])))
    aSDomain_text <- GBFFfeat[aSDomain_start[i]:(aSDomain_start[i] + translation_line[1] + end_line[2] - 2)]

    aSDomain_data[i, "FromTo"] <- strsplit(aSDomain_text, "\\s+")[[1]][3]

    aSDomain_description_line <- grep("/description=", aSDomain_text)
    aSDomain_data[i, "description"] <- paste(gsub("^\\s*/description=\"|\"$", "", aSDomain_text[aSDomain_description_line]), collapse = ", ")

    aSDomain_evalue_line <- grep("/evalue=", aSDomain_text)
    aSDomain_data[i, "evalue"] <- gsub("^\\s*/evalue=\"|\"$", "", aSDomain_text[aSDomain_evalue_line])[1]

    aSDomain_identifier_line <- grep("/identifier=", aSDomain_text)
    aSDomain_data[i, "identifier"] <- gsub("^\\s*/identifier=\"|\"$", "", aSDomain_text[aSDomain_identifier_line])[1]

    aSDomain_label_line <- grep("/label=", aSDomain_text)
    aSDomain_data[i, "label"] <- gsub("^\\s*/label=\"|\"$", "", aSDomain_text[aSDomain_label_line])[1]

    aSDomain_protein_end_line <- grep("/protein_end=", aSDomain_text)
    aSDomain_data[i, "protein_end"] <- gsub("^\\s*/protein_end=\"|\"$", "", aSDomain_text[aSDomain_protein_end_line])[1]

    aSDomain_protein_start_line <- grep("/protein_start=", aSDomain_text)
    aSDomain_data[i, "protein_start"] <- as.numeric(gsub("^\\s*/protein_start=\"|\"$", "", aSDomain_text[aSDomain_protein_start_line]))[1]

    aSDomain_score_line <- grep("/score=", aSDomain_text)
    aSDomain_data[i, "score"] <- gsub("^\\s*/score=\"|\"$", "", aSDomain_text[aSDomain_score_line])[1]

    aSDomain_specificity_line <- grep("/specificity=", aSDomain_text)
    aSDomain_data[i, "specificity"] <- paste(gsub("^\\s*/specificity=\"|\"$", "", aSDomain_text[aSDomain_specificity_line]), collapse = ", ")

    translation_text <- aSDomain_text[translation_line[1]:length(aSDomain_text)]
    aSDomain_data[i, "translation"] <- paste(gsub("\\s+", "", gsub("^\\s*/translation=\"|\"$", "", translation_text)), collapse = "")
  }

  aSDomain_data
}
