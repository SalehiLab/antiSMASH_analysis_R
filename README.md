# antiSMASH Data Extraction and Postprocessing (R)

R utilities for extracting, parsing, and postprocessing [antiSMASH](https://antismash.secondarymetabolites.org/) genome mining output, plus helper functions for retrieving NCBI taxonomy metadata.

## Repository Structure

| File | Description |
|---|---|
| `Functions.R` | Core helper functions used across the pipeline (taxonomy lookup, GenBank file parsing, biosynthetic gene cluster data extraction) |
| `GenomeMining.R` | Main genome mining workflow script |
| `aS_DataExtraction.R` | antiSMASH output data extraction script |

## Dependencies

```r
install.packages("dplyr")
install.packages("rentrez")
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

*More sections (GenomeMining.R, aS_DataExtraction.R, usage workflow) coming next.*

