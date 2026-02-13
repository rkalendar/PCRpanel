# PCRpanel — Custom Amplicon Panel Designer

[![Java](https://img.shields.io/badge/Java-25%2B-orange)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Independent-blue)]()
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

PCRpanel designs **custom amplicon panels** for **NGS** and **Oxford Nanopore (ONT)** sequencing, including **ultra-high multiplex tiling PCR** and standard singleplex/multiplex PCR workflows. It accepts GenBank or FASTA input, automatically detects target regions (exons, introns, promoters, or arbitrary coordinates), and outputs optimised primer sets ready for wet-lab validation.

🌐 **Online version:** <https://primerdigital.com/tools/panel.html>  
👤 **Author:** Ruslan Kalendar  
📧 **Contact:** ruslan.kalendar@helsinki.fi

---

## Table of Contents

- [Key Capabilities](#key-capabilities)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Memory Configuration](#memory-configuration)
- [Configuration Reference](#configuration-reference)
- [Input Formats](#input-formats)
- [Target Region Detection](#target-region-detection)
- [Use Cases](#use-cases)
- [Related Tools](#related-tools)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Key Capabilities

- **Multiplex tiling PCR** and conventional PCR primer design
- **Single-plex panel** mode for non-multiplexed workflows
- Works with **any sequence length** and **any number of targets**
- Target regions can be **exons, introns, promoters, or arbitrary coordinates**
- Design **common primers** based on shared sequences across multiple input files (homology mode)
- Design primers starting from an **existing primer/probe list**
- **Repeated sequences control** to avoid off-target amplification
- Optional alignment against a **reference genome** to check for gene duplications and unknown repeats
- Support for **5′/3′ tails** (adapters, UMIs, barcodes, etc.)

---

## Requirements

| Requirement          | Details                                        |
|----------------------|------------------------------------------------|
| **Operating System** | Platform-independent (Windows, macOS, Linux)   |
| **Java**             | Version 25 or higher                           |

**Java downloads:**

- Oracle JDK: <https://www.oracle.com/java/technologies/downloads/>
- Setting `JAVA_HOME` / `PATH`: <https://www.java.com/en/download/help/path.html>

---

## Installation

### Standard Installation

1. Download or clone this repository.
2. Ensure Java 25+ is installed and available on your `PATH`.
3. The executable JAR is located at `dist/PCRpanel.jar`.

### Installing Java 25 with Conda (Optional)

```bash
# 1. Add conda-forge and set strict channel priority
conda config --add channels conda-forge
conda config --set channel_priority strict

# 2. Create an environment with OpenJDK 25
conda create -n java25 openjdk=25
conda activate java25

# 3. Verify the installation
java -version
```

---

## Quick Start

Run PCRpanel from the command line:

```bash
java -jar <path>/dist/PCRpanel.jar <config.file>
```

**Windows:**

```bash
java -jar C:\PCRpanel\dist\PCRpanel.jar C:\PCRpanel\test\config.file
```

**Linux / macOS:**

```bash
java -jar /data/soft/PCRpanel.jar /data/soft/config.file
```

**Linux / macOS (with memory allocation for large genomes):**

```bash
java -Xms32g -Xmx128g -jar /data/soft/PCRpanel.jar /data/soft/config.file
```

No additional dependencies are required.

---

## Memory Configuration

When using a reference genome via `genome_path`, you must allocate additional heap memory. Depending on genome size, up to 256 GB of RAM may be needed.

```bash
java -Xms32g -Xmx128g -jar PCRpanel.jar config.file
```

| Genome Size    | Recommended RAM | JVM Flags               |
|----------------|-----------------|-------------------------|
| < 100 MB       | Default         | None needed             |
| 100–500 MB     | 16–32 GB        | `-Xms8g -Xmx16g`       |
| 500 MB – 2 GB  | 64 GB           | `-Xms32g -Xmx64g`      |
| > 2 GB         | 128+ GB         | `-Xms64g -Xmx128g`     |

**JVM memory parameters:**

- **`-Xms`** (initial heap size) — memory allocated at startup; avoids allocation delays for large genomes.
- **`-Xmx`** (maximum heap size) — upper limit of heap memory; prevents `OutOfMemoryError`.

---

## Configuration Reference

All parameters are specified in a plain-text configuration file.

### Basic Parameters

| Parameter      | Description                                                          | Default / Example |
|----------------|----------------------------------------------------------------------|-------------------|
| `minPCR`       | Minimum amplicon size (bp)                                           | `250`             |
| `maxPCR`       | Maximum amplicon size (bp)                                           | `500`             |
| `minLen`       | Minimum primer length (nt)                                           | `18`              |
| `maxLen`       | Maximum primer length (nt)                                           | `24`              |
| `minTm`        | Minimum melting temperature (°C)                                     | `60`              |
| `maxTm`        | Maximum melting temperature (°C)                                     | `62`              |
| `minLC`        | Minimum linguistic complexity (%)                                    | `80`              |
| `3end`         | 3′ end constraint                                                    | `w`               |
| `5end`         | 5′ end constraint                                                    | *(empty)*         |
| `forwardtail`  | 5′ adapter for forward primers                                       | Illumina P5 adapter |
| `reversetail`  | 5′ adapter for reverse primers                                       | Illumina P7 adapter |
| `multiplex`    | Generate two overlapping multiplex-compatible panels                 | `true`            |
| `homology`     | Design common primers from shared sequences across input files       | `false`           |

### Input / Output Paths

| Parameter        | Description |
|------------------|-------------|
| `target_path`    | Path to an individual target file (can be specified multiple times) |
| `target_primers` | *(Optional)* Path to an existing primer/probe list |
| `folder_path`    | *(Optional)* Path to a folder of target files (subdirectories included) |
| `folder_out`     | *(Optional)* Output directory for results |
| `genome_path`    | *(Optional)* Path to a folder of reference genome FASTA files (subdirectories included) |

### Example Configuration File (Windows)

```ini
# Target sequences (multiple files supported)
target_path=C:\PCRpanel\test\NG_013019.gb
target_path=C:\PCRpanel\test\NG_011731.gb
target_path=C:\PCRpanel\test\NG_013019.gb
target_path=C:\PCRpanel\test\NG_008847.gb
target_path=C:\PCRpanel\test\NC_000002.gb

# Optional: existing primers to incorporate
target_primers=C:\PCRpanel\test\primers.txt

# Panel mode
homology=false
multiplex=true

# Amplicon size constraints
minPCR=250
maxPCR=500

# Primer parameters
minLen=18
maxLen=24
minTm=60
maxTm=62

# Primer end constraints
3end=w
5end=

# Adapter tails (Illumina example)
forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

### Example Configuration File (Linux)

```ini
# Batch processing
folder_path=/data/genes/
folder_out=/data/report/
genome_path=/data/t2t/

# Optional: existing primers to incorporate
target_primers=/data/primers/primers.txt

# Panel mode
homology=false
multiplex=true

# Amplicon size constraints
minPCR=250
maxPCR=500

# Primer parameters
minLen=18
maxLen=24
minTm=60
maxTm=62

# Primer end constraints
3end=w
5end=

# Adapter tails (Illumina example)
forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

### Batch Processing

```ini
# Reference genome directory (optional)
genome_path=C:\PCRpanel\HumanChromosomes\

# Process all files in a folder
folder_path=C:\PCRpanel\test\

# Specify an output directory
folder_out=C:\PCRpanel\report\
```

**Output behaviour:**

| Configuration                              | Output Location                      |
|--------------------------------------------|--------------------------------------|
| `folder_out` specified                     | All results go to `folder_out`       |
| `folder_out` omitted        | The output directory will be created    |

> **Note:** When both `target_path` and `folder_path` are provided, PCRpanel processes the **union** of all explicitly listed files plus all files discovered in `folder_path`.

---

## Input Formats

### GenBank Files

Recommended (but not strictly required) format. Supported extensions: `*.gb`, `*.gbff`.

PCRpanel parses RefSeqGene / RefSeq records (typically downloaded with `rettype=gbwithparts`). Each file must contain:

- A **FEATURES** table with gene annotations
- An **ORIGIN** section with the nucleotide sequence

> **Fallback:** If a GenBank file lacks an ORIGIN section, PCRpanel falls back to FASTA parsing.

### FASTA Files

Standard FASTA format with one or more records per file. Each record is treated as a separate full-length target region.

---

## Target Region Detection

PCRpanel designs primers **only within target coordinates**. These can be derived automatically from GenBank annotations or specified as arbitrary regions.

### Exon Recognition Priority

Many RefSeq/GenBank records lack explicit `exon` features. PCRpanel applies the following recognition hierarchy:

| Priority | Feature Type                       | Status      | Notes |
|----------|------------------------------------|-------------|-------|
| 1        | `exon`                             | **Active**  | Used directly if present (primary feature) |
| 2        | `mRNA`, `ncRNA`, `rRNA`, `tRNA`    | *Suspended* | `join(…)` blocks extracted as exons |
| 3        | `CDS`                              | *Suspended* | `join(…)` blocks used when transcript features are absent |
| 4        | *Fallback*                         | **Active**  | Full sequence treated as one target region |

> **Note:** Priorities 2 and 3 are currently suspended. Only explicit `exon` features (priority 1) and the full-sequence fallback (priority 4) are active.

This strategy works across RefSeqGene (`NG_*`), genomic (`NC_*`), and transcript records.

### GenBank Coordinate Rules

Coordinates are **1-based, inclusive** (e.g., `5049..5095` includes both endpoints).

**Supported location qualifiers:**

| Qualifier           | Behaviour |
|---------------------|-----------|
| `join(…)`           | Multi-exon features — individual intervals are extracted as targets |
| `order(…)`          | Treated identically to `join` for target extraction |
| `complement(…)`     | Indicates reverse strand; coordinates are extracted normally |
| `<` / `>` (partial) | Partial-bound symbols are ignored; numeric positions are used as-is |

**Examples from FEATURES tables:**

```
# mRNA with multiple exons
mRNA            join(5049..5095,26596..26643,30691..30834,46554..46688)

# CDS on reverse strand
CDS             complement(join(330..403,1050..1120))
```

### Coordinate Conversion (Implementation Note)

To convert GenBank 1-based inclusive intervals (`a..b`) to internal 0-based half-open coordinates (used by Java arrays):

```
start = a − 1    (0-based)
end   = b        (exclusive)
```

**Example:** `5049..5095` → `[5048, 5095)` (length = 47 bp)

---

## Use Cases

### Gene Panel Design

Design primers targeting specific exons across multiple genes for applications such as hereditary disease screening panels, cancer hotspot panels, and pharmacogenomics panels.

### Whole-Genome Tiling Panels

For complete genome coverage (e.g., viral sequencing), the entire genome becomes the target region. PCRpanel designs tiled amplicons across the full sequence, independent of exon annotations.

**Example — SARS-CoV-2 ([NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)):**

```ini
target_path=/path/to/NC_045512.2.gb
minPCR=400
maxPCR=500
```

---

## Related Tools

- **[NCBI RefSeq GenBank Downloader](https://github.com/rkalendar/genbanktools)** — Batch-download GenBank records by gene abbreviation or accession number. Useful for assembling large target lists for PCRpanel.

---

## Troubleshooting

### No exons or targets detected

1. Verify the GenBank record contains both a **FEATURES** table and an **ORIGIN** sequence section.
2. Ensure transcript features use `join(…)` syntax (required for multi-exon genes).
3. Check that the file encoding is UTF-8.

### FASTA fallback behaviour

When FASTA input is supplied, each record is treated as a separate full-length target region. Exon-level targeting is not available without GenBank annotations.

### Java version errors

Ensure Java 25+ is installed:

```bash
java -version
```

If multiple Java versions are installed, verify that `JAVA_HOME` points to the correct version.

### Memory issues with large genomes

Increase the Java heap size as described in [Memory Configuration](#memory-configuration):

```bash
java -Xmx64g -jar PCRpanel.jar config.file
```

---

## Citation

If you use PCRpanel in your research, please cite:

> Kalendar, R. (2025). PCRpanel: Custom Amplicon Panel Designer. Available at <https://primerdigital.com/tools/panel.html>

<!-- TODO: Replace the citation above with the published paper reference once available. -->

---

## License

This project is distributed under the terms of the [MIT License](LICENSE).

---

<p align="center">
  <em>PCRpanel — Designed for researchers, by researchers.</em>
</p>
