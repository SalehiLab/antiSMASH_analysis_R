# antiSMASH Data Extraction and Postprocessing (R)

R utilities for extracting, parsing, and postprocessing [antiSMASH](https://antismash.secondarymetabolites.org/) genome mining output, plus helper functions for retrieving NCBI taxonomy metadata.

## Repository Structure

```
.
├── R/
│   ├── io_read.R           # File reading: GenBank files, FASTA headers, KnownClusterBlast reports
│   ├── ncbi.R                # NCBI taxonomy lookup (rentrez)
│   ├── genome_stats.R         # Sequence extraction, length, GC content
│   ├── extract_features.R     # antiSMASH feature parsers (protocluster, proto_core, cand_cluster, PFAM, aSDomain)
├── run_genome_mining.R      # Orchestration: exploratory genome mining pipeline
└── run_data_extraction.R    # Orchestration: genome summary + full BGC annotation extraction
```

Everything in `R/` is a pure function library — no file I/O side effects beyond what each function is explicitly meant to do, and no top-level script logic. The two `run_*.R` scripts at the repo root are the only files that actually execute a pipeline: they `source()` the relevant files in `R/`, call functions in order, and write the numbered output CSVs.


## Dependencies

```r
install.packages(c("dplyr", "rentrez", "stringr", "tidyr"))

# Biostrings is distributed via Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```

---

## `R/io_read.R`

Low-level file reading utilities. No antiSMASH-domain parsing logic lives here — just getting raw text and records off disk.

### `readGBFF(file, text = readLines(file), verbose = FALSE)`

Reads a GenBank flat file (GBFF) and splits it into a list of individual records, using each `LOCUS` line as a record boundary. This allows multi-record GBFF files to be processed one record at a time.

**Returns:** a list of character vectors, one per GenBank record.

```r
records <- readGBFF("cluster_001.region001.gbk")
```

### `parse_fasta_headers(fna_files)`

Parses the `>accession, description` header lines of one or more FASTA (`.fna`) files into a lookup table.

**Returns:** a data frame with columns `refseq`, `Deffinition`.

```r
fna_files <- list.files("./input", pattern = "*.fna", full.names = TRUE)
ref_data <- parse_fasta_headers(fna_files)
```

### `load_merged_gbff(pattern = "^merged.*\\.gbk$")`

Loads all antiSMASH `merged*.gbk` files in the working directory, names each record by its `VERSION` accession, and extracts the feature block between `LOCUS` and `ORIGIN` for each record.

**Returns:** a named list of character vectors (one per genome), ready to pass into the `extract_*_data()` functions.

```r
GBFFfeat <- load_merged_gbff()
```

### `read_knownclusterblast(dir = "./knownclusterblast/")`

Reads antiSMASH KnownClusterBlast `.txt` reports from a directory and extracts the `Source:` compound name reported under each cluster's `Details:` section.

**Returns:** a data frame with columns `refseq`, `Protocluster_number`, `Compounds`.

```r
compounds <- read_knownclusterblast("./knownclusterblast/")
```

---

## `R/ncbi.R`

### `get_taxonomy(refseq_ids)`

Retrieves the full taxonomic lineage for one or more NCBI RefSeq accessions by querying the `nuccore` database via `rentrez::entrez_fetch()` and parsing the `ORGANISM` field from the returned GenBank record.

**Returns:** a character vector `c(refseq_ids, taxonomy)`.

```r
get_taxonomy("NC_003888.3")
```

---

## `R/genome_stats.R`

### `add_genome_stats(ref_data, fna_files)`

Given a `refseq`/`Deffinition` table (from `parse_fasta_headers()`) and the FASTA files it came from, extracts each genome's raw sequence, computes its length and GC content, and looks up its taxonomy via `get_taxonomy()`.

**Depends on:** `get_taxonomy()` (`R/ncbi.R`) — source that file first.

**Returns:** `ref_data` with `sequence`, `Length`, `GC_content`, and `taxonomy` columns added (all coerced to character).

```r
ref_data <- parse_fasta_headers(fna_files)
ref_data <- add_genome_stats(ref_data, fna_files)
```

---

## `R/extract_features.R`

The core antiSMASH domain-knowledge parsers. Each function takes one genome's GBK feature block (as produced by `load_merged_gbff()`), plus a `refseq` and `Deffinition` to tag rows with, and returns a tidy data frame — one row per feature instance found. All five share the same underlying pattern: locate a feature block by its antiSMASH feature key, then pull out `/key="value"` qualifiers by regex.

| Function | Feature parsed | Columns returned |
|---|---|---|
| `extract_protocluster_data()` | `protocluster` | `refseq`, `Deffinition`, `Length`, `Protocluster_number`, `FromTo`, `Product`, `Category`, `Core_location`, `aStool`, `Contig_edge`, `Cutoff`, `Tool`, `Neighbourhood`, `Detection_rule` |
| `extract_proto_core_data()` | `proto_core` | `refseq`, `Deffinition`, `Length`, `Protocluster_number`, `FromTo`, `Product`, `aStool`, `Cutoff`, `Tool`, `Neighbourhood`, `Detection_rule` |
| `extract_cand_cluster_data()` | `cand_cluster` | `refseq`, `Deffinition`, `FromTo`, `SMILES`, `candidate_cluster_number`, `product`, `protocluster_number` |
| `extract_PFAM_data()` | `PFAM_domain` | `refseq`, `Deffinition`, `FromTo`, `description`, `evalue`, `gene_ontologies`, `label`, `protein_start`, `protein_end`, `score`, `translation` |
| `extract_aSDomain_data()` | `aSDomain` | `refseq`, `Deffinition`, `FromTo`, `description`, `evalue`, `identifier`, `label`, `protein_start`, `protein_end`, `score`, `specificity`, `translation` |

Every function returns a single `NA` row if the corresponding feature isn't present in that genome's record, so batch calls via `mapply()` never fail on an empty genome.

```r
GBFFfeat <- load_merged_gbff()
protoclusters <- mapply(extract_protocluster_data, GBFFfeat, ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
combined <- do.call("rbind", protoclusters)
```


---

## `run_genome_mining.R`

End-to-end pipeline producing tidy CSV tables covering biosynthetic gene clusters, candidate clusters, protein domains, and known-cluster hits — including a filtered shortlist of PKS candidates.

### Expected input structure

```
working directory/
├── input/
│   └── *.fna                     # FASTA files; headers used to map refseq -> Deffinition
├── merged*.gbk                   # antiSMASH merged GenBank output, one per genome
└── knownclusterblast/
    └── *.txt                     # antiSMASH KnownClusterBlast hit reports
```

### Workflow

1. **Parse input FASTA headers** (`parse_fasta_headers()`). Saved as `1_input_data.csv`.
2. **Load merged GenBank files** (`load_merged_gbff()`).
3. **Extract protocluster data** (`extract_protocluster_data()`). Saved as `2_results_protocluster.csv`.
4. **Extract candidate cluster data** (`extract_cand_cluster_data()`), including multi-line SMILES strings. Saved as `3_results_cand_cluster.csv`.
5. **Extract PFAM domain data** (`extract_PFAM_data()`). Saved as `4_results_PFAM_data.csv`.
6. **Extract antiSMASH domain data** (`extract_aSDomain_data()`). Saved as `5_results_aSDomain_data.csv`.
7. **Parse KnownClusterBlast hits** (`read_knownclusterblast()`).
8. **Merge results** — joins the protocluster table (step 3) with the compound hits (step 7) on `refseq` and `Protocluster_number`. Saved as `2_mergedResults.csv`.

### Output files

| File | Contents |
|---|---|
| `1_input_data.csv` | RefSeq accession ↔ definition lookup from input FASTA headers |
| `2_results_protocluster.csv` | One row per protocluster: coordinates, product, category, core location, detection tool/rule |
| `2_mergedResults.csv` | Protocluster data merged with KnownClusterBlast compound hits |
| `3_results_cand_cluster.csv` | One row per candidate cluster: coordinates, SMILES, product, linked protoclusters |
| `4_results_PFAM_data.csv` | PFAM domain annotations (description, e-value, GO terms, coordinates, translation) |
| `5_results_aSDomain_data.csv` | antiSMASH domain annotations (description, e-value, specificity, coordinates, translation) |

### Usage

```r
# Run from a directory containing ./input/*.fna, merged*.gbk, and ./knownclusterblast/*.txt
source("run_genome_mining.R")
```

---

## `run_data_extraction.R`

Builds a genome-level summary table (sequence, length, GC content, taxonomy) from input FASTA files, then extracts the full protocluster and proto_core annotations from merged antiSMASH GenBank output.

### Expected input structure

```
working directory/
├── input/
│   └── *.fna              # FASTA files; headers used to map refseq -> Deffinition
└── merged*.gbk             # antiSMASH merged GenBank output, one per genome
```

### Workflow

1. **Parse input FASTA headers** (`parse_fasta_headers()`).
2. **Compute genome sequence stats** — sequence, length, GC%, taxonomy (`add_genome_stats()`).
3. **Save genome summary tables** — full table (including sequence) to `1_input_data_sequence.csv`; a lighter version without the sequence column to `1_input_data.csv`.
4. **Load merged GenBank files** (`load_merged_gbff()`).
5. **Extract protocluster data** (`extract_protocluster_data()`). Saved as `2_results_protocluster.csv`.
6. **Extract proto_core data** (`extract_proto_core_data()`). Saved as `3_results_protocore.csv`.

### Output files

| File | Contents |
|---|---|
| `1_input_data_sequence.csv` | Genome summary including raw sequence: refseq, definition, length, GC%, taxonomy, sequence |
| `1_input_data.csv` | Same genome summary, without the sequence column |
| `2_results_protocluster.csv` | Full protocluster annotations (product, category, core location, detection tool/rule, cutoff, etc.) |
| `3_results_protocore.csv` | Full proto_core annotations (product, detection tool/rule, cutoff, neighbourhood, etc.) |

### Usage

```r
# Run from a directory containing ./input/*.fna and merged*.gbk
source("run_data_extraction.R")
```

---

## End-to-End Workflow

```r
source("run_genome_mining.R")    # exploratory pipeline: protoclusters, candidate clusters,
                                  # PFAM/aSDomain annotations, PKS shortlist, KnownClusterBlast hits

source("run_data_extraction.R")  # genome-level summary (sequence, GC%, taxonomy) +
                                  # full protocluster / proto_core annotation tables
```

Each script sources everything it needs from `R/` internally — you don't need to `source()` the `R/` files yourself before running either pipeline. Both scripts expect the same general input layout: FASTA files in `./input` and antiSMASH `merged*.gbk` files in the working directory; `run_genome_mining.R` additionally expects KnownClusterBlast reports in `./knownclusterblast/`.

## Running Step by Step (Interactively)

Because every piece of logic is now a standalone function in `R/`, you don't have to run a full `run_*.R` script end to end — you can source just the pieces you need and call them one at a time in the console, inspecting each result before moving on. This is the recommended way to explore a new dataset or debug a specific step.

### `run_genome_mining.R`, step by step

```r
library(dplyr)
source("R/io_read.R")
source("R/ncbi.R")
source("R/extract_features.R")
source("R/filters.R")

# Step 1: parse FASTA headers
fna_files <- list.files(path = "./input", pattern = "*.fna", full.names = TRUE)
ref_data <- parse_fasta_headers(fna_files)
head(ref_data)                       # <- inspect before continuing

# Step 2: load merged GenBank feature blocks
GBFFfeat <- load_merged_gbff()
names(GBFFfeat)                      # <- check which genomes loaded

# Step 3: extract protoclusters for all genomes
protoclusters <- mapply(extract_protocluster_data, GBFFfeat,
                         ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
protoclusters_df <- do.call("rbind", protoclusters)
View(protoclusters_df)               # <- inspect before writing out

# Step 4: extract candidate clusters, then filter to PKS candidates
cand_clusters <- mapply(extract_cand_cluster_data, GBFFfeat,
                         ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE)
cand_clusters_df <- do.call("rbind", cand_clusters)
pks_candidates <- select_pks_candidates(cand_clusters_df)
nrow(pks_candidates)                 # <- sanity-check the filter before saving

# ...continue with extract_PFAM_data(), extract_aSDomain_data(),
# read_knownclusterblast(), and merge() as needed, writing CSVs only
# for the steps you actually want to keep.
```

### `run_data_extraction.R`, step by step

```r
library(rentrez); library(stringr); library(dplyr); library(tidyr); library(Biostrings)
source("R/io_read.R")
source("R/ncbi.R")
source("R/genome_stats.R")
source("R/extract_features.R")

fna_files <- list.files(path = "./input", pattern = "*.fna", full.names = TRUE)
ref_data <- parse_fasta_headers(fna_files)

# add_genome_stats() calls get_taxonomy() per genome, which hits NCBI --
# run this on a small subset first if you're testing
ref_data <- add_genome_stats(ref_data, fna_files)
head(ref_data[, c("refseq", "Length", "GC_content", "taxonomy")])

GBFFfeat <- load_merged_gbff()
protoclusters_df <- do.call("rbind", mapply(extract_protocluster_data, GBFFfeat,
                             ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE))
protocore_df <- do.call("rbind", mapply(extract_proto_core_data, GBFFfeat,
                         ref_data$refseq, ref_data$Deffinition, SIMPLIFY = FALSE))
```

**Why this matters in practice:** `add_genome_stats()` calls NCBI once per genome via `get_taxonomy()`, which can be slow or hit rate limits on a large batch. Running it interactively — one genome or a small subset at a time — makes it much easier to catch issues (missing accessions, network errors) before committing to a full run.



