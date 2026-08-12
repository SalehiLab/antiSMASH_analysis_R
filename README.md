# antiSMASH Data Extraction and Postprocessing (R)

R utilities for extracting, parsing, and postprocessing [antiSMASH](https://antismash.secondarymetabolites.org/) genome mining output, plus helper functions for retrieving NCBI taxonomy metadata.

## Repository Structure

| File | Description |
|---|---|
| `Functions.R` | Core helper functions used across the pipeline (taxonomy lookup, GenBank file parsing, biosynthetic gene cluster data extraction). Sourced by the other two scripts. |
| `GenomeMining.R` | Exploratory genome mining pipeline: protoclusters, candidate clusters, PFAM/aSDomain annotations, PKS candidate shortlisting, KnownClusterBlast hits |
| `aS_DataExtraction.R` | Genome-level summary (sequence, GC%, taxonomy) plus full protocluster/proto_core annotation extraction |

## Dependencies

```r
install.packages(c("dplyr", "rentrez", "stringr", "tidyr"))

# Biostrings is distributed via Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```

---

## `Functions.R`

Four core functions supporting the antiSMASH postprocessing pipeline: one for NCBI taxonomy retrieval, one for parsing multi-record GenBank flat files, and two for extracting biosynthetic gene cluster (BGC) region annotations produced by antiSMASH.

### `get_taxonomy(refseq_ids)`

Retrieves the full taxonomic lineage for one or more NCBI RefSeq accessions by querying the `nuccore` database via `rentrez::entrez_fetch()` and parsing the `ORGANISM` field from the returned GenBank record.

**Arguments**
- `refseq_ids` — a RefSeq accession (or vector of accessions) to query.

**Returns**
A character vector: `c(refseq_ids, taxonomy)`, where `taxonomy` is the concatenated, whitespace-stripped lineage string.

**Example**
```r
get_taxonomy("NC_003888.3")
```

---

### `readGBFF(file, text = readLines(file), verbose = FALSE)`

Reads a GenBank flat file (GBFF) and splits it into a list of individual records, using each `LOCUS` line as a record boundary. This allows multi-record GBFF files (e.g. antiSMASH region outputs) to be processed one record at a time by the extraction functions below.

**Arguments**
- `file` — path to a `.gbk`/`.gbff` file.
- `text` — optionally, pre-loaded lines of the file (defaults to `readLines(file)`).
- `verbose` — currently unused; reserved for future logging.

**Returns**
A list of character vectors, each containing the lines of a single GenBank record.

**Example**
```r
records <- readGBFF("cluster_001.region001.gbk")
```

---

### `extract_protocluster_data(GBFFfeat, refseq, Deffinition)`

Parses the `protocluster` features of an antiSMASH-annotated GenBank record and returns one row of metadata per protocluster (e.g. NRPS, PKS, terpene, etc.), including product type, category, genomic coordinates, detection tool/rule, and cutoff values.

**Arguments**
- `GBFFfeat` — a GenBank record (or list containing one, as returned by `readGBFF()`).
- `refseq` — accession/identifier to tag output rows with.
- `Deffinition` — definition/description string to tag output rows with.

**Returns**
A data frame with one row per protocluster and the following columns:

| Column | Description |
|---|---|
| `refseq` | Accession passed in |
| `Deffinition` | Definition passed in |
| `Length` | Sequence length (from `LOCUS` line) |
| `Protocluster_number` | antiSMASH protocluster index |
| `FromTo` | Genomic coordinates of the protocluster |
| `Product` | Predicted BGC product type |
| `Category` | antiSMASH product category |
| `Core_location` | Coordinates of the core region |
| `aStool` | antiSMASH tool version tag |
| `Contig_edge` | Whether the cluster touches a contig edge (`True`/`False`) |
| `Cutoff` | Detection cutoff distance (bp) |
| `Tool` | Detection tool used |
| `Neighbourhood` | Neighbourhood extension distance (bp) |
| `Detection_rule` | Rule string used for detection |

Returns a single-row `NA` data frame if no protocluster features are found.

**Example**
```r
records <- readGBFF("cluster_001.region001.gbk")
protoclusters <- extract_protocluster_data(records, refseq = "NC_003888.3", Deffinition = "Streptomyces coelicolor")
```

---

### `extract_proto_core_data(GBFFfeat, refseq, Deffinition)`

Parses the `proto_core` features (the core biosynthetic region within each protocluster) of an antiSMASH-annotated GenBank record and returns one row of metadata per proto_core.

**Arguments**
Same as `extract_protocluster_data()`.

**Returns**
A data frame with one row per proto_core and the following columns:

| Column | Description |
|---|---|
| `refseq` | Accession passed in |
| `Deffinition` | Definition passed in |
| `Length` | Sequence length (from `LOCUS` line) |
| `Protocluster_number` | Associated protocluster index |
| `FromTo` | Genomic coordinates of the proto_core |
| `Product` | Predicted BGC product type(s), semicolon-separated |
| `aStool` | antiSMASH tool version tag |
| `Cutoff` | Detection cutoff distance (bp) |
| `Tool` | Detection tool used |
| `Neighbourhood` | Neighbourhood extension distance (bp) |
| `Detection_rule` | Rule string used for detection |

Returns a single-row `NA` data frame if no proto_core features are found.

**Example**
```r
records <- readGBFF("cluster_001.region001.gbk")
proto_cores <- extract_proto_core_data(records, refseq = "NC_003888.3", Deffinition = "Streptomyces coelicolor")
```

---

## `GenomeMining.R`

End-to-end pipeline that parses merged antiSMASH GenBank output for a batch of genomes and produces a set of tidy CSV tables covering biosynthetic gene clusters, candidate clusters, protein domains, and known-cluster hits — including a filtered shortlist of polyketide (PKS) candidate clusters based on SMILES composition.

### Expected input structure

```
working directory/
├── input/
│   └── *.fna                     # FASTA files; headers used to map refseq -> explanation
├── merged*.gbk                   # antiSMASH merged GenBank output, one per genome
└── knownclusterblast/
    └── *.txt                     # antiSMASH KnownClusterBlast hit reports
```

### Workflow

1. **Parse input FASTA headers** — reads every `.fna` file in `./input`, splits `>` header lines to build a `refseq` / `explanation` lookup table (`ref_data`). Saved as `1_input_data.csv`.
2. **Load merged GenBank files** — reads all `merged*.gbk` files with `readGBFF()`, names each record by its `VERSION` accession, and extracts the feature block between `LOCUS` and `ORIGIN` (`GBFFfeat`).
3. **Extract protocluster data** — `extract_protocluster_data()` (a lighter 7-column variant of the one in `Functions.R`) pulls protocluster number, coordinates, product, category, and core location for every genome. Saved as `2_results_protocluster.csv`.
4. **Extract candidate cluster data** — `extract_cand_cluster_data()` parses `cand_cluster` features, including multi-line SMILES strings, candidate cluster number, product, and linked protocluster number. Saved as `3_results_cand_cluster.csv`.
5. **Filter PKS candidates** — candidate clusters are screened by SMILES composition (12–20 carbon atoms, ≤4 double bonds `=`, no nitrogen `N`) and restricted to `T1PKS`/`T2PKS`/`T3PKS`/`PKS-like` products. Saved as `3_results_cand_cluster_Selected.csv`.
6. **Extract PFAM domain data** — `extract_PFAM_data()` parses `PFAM_domain` features (description, e-value, GO terms, protein coordinates, score, translation). Saved as `4_results_PFAM_data.csv`.
7. **Extract antiSMASH domain data** — `extract_aSDomain_data()` parses `aSDomain` features (description, e-value, identifier, specificity, protein coordinates, score, translation). Saved as `5_results_aSDomain_data.csv`.
8. **Parse KnownClusterBlast hits** — reads each `.txt` report in `./knownclusterblast/`, extracts the `Source:` compound name(s) following the `Details:` line, and links them back to `refseq` + `protocluster_number`.
9. **Merge results** — joins the protocluster table (step 3) with the compound hits (step 8) on `refseq` and `protocluster_number`. Saved as `2_mergedResults.csv`.

### Output files

| File | Contents |
|---|---|
| `1_input_data.csv` | RefSeq accession ↔ explanation lookup from input FASTA headers |
| `2_results_protocluster.csv` | One row per protocluster: coordinates, product, category, core location |
| `2_mergedResults.csv` | Protocluster data merged with KnownClusterBlast compound hits |
| `3_results_cand_cluster.csv` | One row per candidate cluster: coordinates, SMILES, product, linked protoclusters |
| `3_results_cand_cluster_Selected.csv` | Candidate clusters filtered to likely PKS compounds by SMILES composition |
| `4_results_PFAM_data.csv` | PFAM domain annotations (description, e-value, GO terms, coordinates, translation) |
| `5_results_aSDomain_data.csv` | antiSMASH domain annotations (description, e-value, specificity, coordinates, translation) |

### Usage

```r
# Run from a directory containing ./input/*.fna, merged*.gbk, and ./knownclusterblast/*.txt
source("Functions.R")
source("GenomeMining.R")
```

---

## `aS_DataExtraction.R`

Builds a genome-level summary table (sequence, length, GC content, taxonomy) from input FASTA files, then extracts full protocluster and proto_core annotations from merged antiSMASH GenBank output using the functions defined in `Functions.R`.

> **Depends on `Functions.R`** — this script calls `get_taxonomy()`, `readGBFF()`, `extract_protocluster_data()`, and `extract_proto_core_data()`, so `Functions.R` must be sourced first.

### Expected input structure

```
working directory/
├── input/
│   └── *.fna              # FASTA files; headers used to map refseq -> Deffinition
└── merged*.gbk             # antiSMASH merged GenBank output, one per genome
```

### Workflow

1. **Parse input FASTA headers** — reads every `.fna` file in `./input` and builds a `refseq` / `Deffinition` lookup table (`ref_data`), identical in approach to `GenomeMining.R`.
2. **Extract genome sequences** — for each FASTA record, concatenates the sequence lines following its header and stores the full sequence in `ref_data$sequence`.
3. **Compute genome statistics** — calculates sequence `Length` and `GC_content` (%) for each genome from its base composition.
4. **Look up taxonomy** — calls `get_taxonomy()` (from `Functions.R`) for every `refseq` accession and appends the parsed lineage as `ref_data$taxonomy`.
5. **Save genome summary tables** — writes the full table (including raw sequence) to `1_input_data_sequence.csv`, and a lighter version without the sequence column to `1_input_data.csv`.
6. **Load merged GenBank files** — reads all `merged*.gbk` files with `readGBFF()`, names each record by its `VERSION` accession, and extracts the feature block between `LOCUS` and `ORIGIN` (`GBFFfeat`).
7. **Extract protocluster data** — applies `extract_protocluster_data()` (the full 14-column version from `Functions.R`) across all genomes. Saved as `2_results_protocluster.csv`.
8. **Extract proto_core data** — applies `extract_proto_core_data()` across all genomes. Saved as `3_results_protocore.csv`.

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
source("Functions.R")
source("aS_DataExtraction.R")
```

---

## End-to-End Workflow

The three scripts build on each other and are intended to be run against antiSMASH output as follows:

```r
source("Functions.R")          # load shared helper functions

source("GenomeMining.R")       # exploratory pipeline: protoclusters, candidate clusters,
                                # PFAM/aSDomain annotations, PKS shortlist, KnownClusterBlast hits

source("aS_DataExtraction.R")  # genome-level summary (sequence, GC%, taxonomy) +
                                # full protocluster / proto_core annotation tables
```

Both `GenomeMining.R` and `aS_DataExtraction.R` expect the same general input layout: FASTA files in `./input`, antiSMASH `merged*.gbk` files in the working directory, and (for `GenomeMining.R`) KnownClusterBlast reports in `./knownclusterblast/`. Each script can be run independently once `Functions.R` has been sourced, and each writes its outputs as numbered CSV files in the working directory.

