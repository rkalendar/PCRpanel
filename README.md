# PCRpanel — Custom Amplicon Panel Designer

[![Java](https://img.shields.io/badge/Java-25%2B-orange)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Independent-blue)]()

PCRpanel designs **custom amplicon panels** for **NGS** and **Oxford Nanopore (ONT)**, including **ultra-high multiplex tiling PCR** and standard singleplex/multiplex PCR tasks.

🌐 **Online version:** https://primerdigital.com/tools/panel.html  
👤 **Author:** Ruslan Kalendar  
📧 **Contact:** ruslan.kalendar@helsinki.fi

---

## Table of Contents

- [Key Capabilities](#key-capabilities)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Standard Installation](#standard-installation)
  - [Installing Java 25 with Conda](#installing-java-25-with-conda-optional)
- [Quick Start](#quick-start)
- [Configuration Reference](#configuration-reference)
  - [Basic Parameters](#basic-parameters)
  - [Input/Output Paths](#inputoutput-paths)
- [Input Formats](#input-formats)
  - [GenBank Files](#genbank-files)
  - [FASTA Files](#fasta-files)
- [Target Region Detection](#target-region-detection)
  - [Exon Recognition Priority](#exon-recognition-priority)
  - [GenBank Coordinate Rules](#genbank-coordinate-rules)
  - [Coordinate Conversion](#coordinate-conversion-implementation-note)
- [Use Cases](#use-cases)
  - [Gene Panel Design](#gene-panel-design)
  - [Whole-Genome Tiling](#whole-genome-tiling-panels)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Key Capabilities

- **Multiplex tiling PCR** and conventional PCR primer design
- Works with **any sequence length** and **any number of targets**
- Target regions can be **exons, introns, promoters, or arbitrary coordinates**
- Design primers starting from **existing primer/probe lists**
- Supports optional **5′/3′ tails** (adapters, UMIs, barcodes, etc.)

---

## Requirements

| Requirement | Details |
|-------------|---------|
| **Operating System** | Platform independent (Windows, macOS, Linux) |
| **Java** | Version 25 or higher |

**Java Downloads:**
- Oracle JDK: https://www.oracle.com/java/technologies/downloads/
- Setting JAVA_HOME/PATH: https://www.java.com/en/download/help/path.html

---

## Installation

### Standard Installation

1. Download or clone this repository
2. Ensure Java 25+ is installed and available in your PATH
3. The executable JAR is located at `dist/PCRpanel.jar`

### Installing Java 25 with Conda (Optional)

1. Add conda-forge and set strict priority:

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

2. Create an environment and install OpenJDK 25:

```bash
conda create -n java25 openjdk=25
conda activate java25
```

3. Verify the installation:

```bash
java -version
```

---

## Quick Start

Run PCRpanel from the command line:

```bash
java -jar <path>/dist/PCRpanel.jar <config.file>
```

**Example (Windows):**

```bash
java -jar C:\PCRpanel\dist\PCRpanel.jar C:\PCRpanel\test\config.file
```

**Example (Linux/macOS):**

```bash
java -jar ~/PCRpanel/dist/PCRpanel.jar ~/PCRpanel/test/config.file
```

---

## Configuration Reference

All parameters are specified in a plain-text configuration file.

### Basic Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `minPCR` | Minimum amplicon size (bp) | `250` |
| `maxPCR` | Maximum amplicon size (bp) | `500` |
| `minLen` | Minimum primer length (nt) | `18` |
| `maxLen` | Maximum primer length (nt) | `24` |
| `minTm` | Minimum melting temperature (°C) | `60` |
| `maxTm` | Maximum melting temperature (°C) | `62` |
| `minLC` | Minimum linguistic complexity (%) | `80` |
| `3end` | 3′ end constraint | `w` |
| `5end` | 5′ end constraint | *(empty)* |
| `forwardtail` | 5′ adapter for forward primers | Illumina P5 adapter |
| `reversetail` | 5′ adapter for reverse primers | Illumina P7 adapter |

### Input/Output Paths

| Parameter | Description |
|-----------|-------------|
| `target_path` | Path to individual target file (can be repeated) |
| `target_primers` | Path to existing primer/probe list |
| `folder_path` | Path to folder containing multiple target files |
| `folder_out` | Output directory for results |

**Example Configuration File:**

```ini
# Target sequences (multiple files supported)
target_path=C:\PCRpanel\test\NG_008690.txt
target_path=C:\PCRpanel\test\NG_011731.txt
target_path=C:\PCRpanel\test\NG_013019.txt
target_path=C:\PCRpanel\test\NG_008847.txt
target_path=C:\PCRpanel\test\NC_000002.txt

# Optional: existing primers to incorporate
target_primers=C:\PCRpanel\test\primers.txt

# Amplicon size constraints
minPCR=250
maxPCR=500

# Primer parameters
minLen=18
maxLen=24
minTm=60
maxTm=62
minLC=78

# Primer end constraints
3end=w
5end=

# Adapter tails (Illumina example)
forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

### Batch Processing

**Using a folder for input:**

```ini
folder_path=C:\PCRpanel\test\
```

**Specifying output directory:**

```ini
folder_out=C:\PCRpanel\report\
```

**Output behavior:**

| Configuration | Output Location |
|---------------|-----------------|
| `folder_out` specified | All results go to `folder_out` |
| `folder_out` not specified + `target_path` | Same directory as each input file |
| `folder_out` not specified + `folder_path` | Input folder |

> **Note:** When both `target_path` and `folder_path` are provided, PCRpanel processes the **union** of all listed files plus all files discovered in `folder_path`.

---

## Input Formats

### GenBank Files

**Recommended format.** Supported extensions: `*.gb`, `*.gbff`

PCRpanel parses RefSeqGene/RefSeq records (typically downloaded with `rettype=gbwithparts`). The file must contain:

- A **FEATURES** table with gene annotations
- An **ORIGIN** section with the nucleotide sequence

### FASTA Files

Standard FASTA format with one or multiple records per file. Each FASTA record is treated as a separate target (full-length region).

> **Fallback behavior:** If a GenBank file lacks an ORIGIN sequence section, PCRpanel attempts FASTA parsing.

---

## Target Region Detection

PCRpanel designs primers **only within target coordinates**. These can be derived automatically from GenBank annotations or specified as arbitrary regions.

### Exon Recognition Priority

Many RefSeq/GenBank records lack explicit `exon` features. PCRpanel uses this recognition hierarchy:

| Priority | Feature Type | Notes |
|----------|--------------|-------|
| 1 | `exon` | Used directly if present |
| 2 | `mRNA`, `ncRNA`, `rRNA`, `tRNA` | `join(...)` blocks extracted as exons |
| 3 | `CDS` | `join(...)` blocks used when transcript features absent |
| 4 | *Fallback* | Full sequence treated as one target region |

This strategy works across RefSeqGene (NG_\*), genomic (NC_\*), and transcript records.

### GenBank Coordinate Rules

- Locations are **1-based, inclusive** (e.g., `5049..5095` includes both endpoints)
- Supported location wrappers:
  - `join(...)` — multi-exon features
  - `order(...)` — treated like join for target extraction
  - `complement(...)` — indicates reverse strand; coordinates extracted normally
- Partial bounds (`<123..456` or `123..>456`) — the `<`/`>` symbols are ignored; numeric position is used

**Examples from FEATURES tables:**

```
# mRNA with multiple exons
mRNA            join(5049..5095,26596..26643,30691..30834,46554..46688)

# CDS on reverse strand
CDS             complement(join(330..403,1050..1120))
```

The extracted target blocks are the individual intervals inside `join(...)`.

### Coordinate Conversion (Implementation Note)

For internal 0-based half-open coordinates (common in Java arrays), convert GenBank intervals `a..b` as follows:

```
start = a - 1    (0-based)
end   = b        (exclusive)
```

**Example:** `5049..5095` → `[5048, 5095)` (length: 47 bp)

---

## Use Cases

### Gene Panel Design

Design primers targeting specific exons across multiple genes for applications such as:

- Hereditary disease screening panels
- Cancer hotspot panels
- Pharmacogenomics panels

### Whole-Genome Tiling Panels

For complete genome coverage (e.g., viral sequencing), the entire genome becomes the target region.

**Example: SARS-CoV-2 (NC_045512.2)**

```ini
target_path=/path/to/NC_045512.2.gb
minPCR=400
maxPCR=500
```

Reference: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2

PCRpanel designs tiled amplicons across the full sequence, independent of exon annotations.

---

## Troubleshooting

### No exons/targets detected

**Check the following:**

1. The GenBank record contains both a **FEATURES** table and an **ORIGIN** sequence section
2. Transcript features use `join(...)` syntax (required for multi-exon genes)
3. File encoding is correct (UTF-8 recommended)

### FASTA fallback behavior

If you supply FASTA input, each record is treated as a separate full-length target region. No exon-level targeting is available without GenBank annotations.

### Java version errors

Ensure Java 25+ is installed:

```bash
java -version
```

If multiple Java versions are installed, verify `JAVA_HOME` points to the correct version.

### Memory issues with large genomes

For large genomes, increase Java heap size:

```bash
java -Xmx4g -jar PCRpanel.jar config.file
```

---

---

## Citation

If you use PCRpanel in your research, please cite:

**

---

## Contributing

Contributions are welcome! Please submit issues and pull requests to the project repository.

---

<p align="center">
  <i>PCRpanel — Designed for researchers, by researchers</i>
</p>
