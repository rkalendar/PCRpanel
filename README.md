# PCRpanel — Custom Amplicon Panel Designer

[![Java](https://img.shields.io/badge/Java-25%2B-orange?logo=openjdk)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-blue)]()
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Docs](https://img.shields.io/badge/Docs-Online-informational)](https://primerdigital.com/tools/panel.html)

**PCRpanel** designs custom amplicon panels for **next-generation sequencing (NGS)** and **Oxford Nanopore (ONT)** platforms. It supports **ultra-high multiplex tiling PCR**, standard singleplex workflows, and everything in between. Give it a GenBank or FASTA file, and it returns optimised primer sets ready for wet-lab validation.

🌐 **Online version:** <https://primerdigital.com/tools/panel.html>
👤 **Author:** Ruslan Kalendar
📧 **Contact:** ruslan.kalendar@helsinki.fi

---

## Table of Contents

- [Key Capabilities](#key-capabilities)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Memory Configuration](#memory-configuration)
- [Configuration Reference](#configuration-reference)
- [Input Formats](#input-formats)
- [Target Region Detection](#target-region-detection)
- [Output Files](#output-files)
- [Use Cases](#use-cases)
- [Related Tools](#related-tools)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Citation](#citation)
- [License](#license)

---

## Key Capabilities

| Feature | Description |
|---------|-------------|
| **Multiplex tiling PCR** | Designs two overlapping pools of amplicons for complete, gap-free target coverage |
| **Single-plex mode** | Non-multiplexed primer design for simpler workflows |
| **Any-scale input** | Works with sequences of any length and any number of targets |
| **Flexible targeting** | Target exons, introns, promoters, UTRs, or arbitrary coordinate ranges |
| **Homology mode** | Designs common (consensus) primers from shared regions across multiple input files |
| **Existing primer import** | Starts from your current primer/probe list and fills gaps |
| **Repeat filtering** | Avoids off-target amplification by screening for repeated sequences |
| **Genome alignment** | Optional alignment against a reference genome to flag gene duplications and unknown repeats |
| **Adapter support** | Attach 5′/3′ tails — adapters, UMIs, barcodes, indexing sequences — to any primer |

---

## How It Works

```
┌─────────────────────────────────────────────────────────────────────┐
│                         PCRpanel Workflow                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│   1. INPUT                 2. TARGET DETECTION                      │
│   ┌──────────────┐         ┌──────────────────────┐                 │
│   │ GenBank (.gb)│────────▶│ Parse exons, introns,│                │
│   │ FASTA  (.fa) │         │ CDS, or full sequence │                │
│   └──────────────┘         └──────────┬───────────┘                 │
│                                       │                             │
│   3. PRIMER DESIGN                    ▼                             │
│   ┌──────────────────────────────────────────────────┐              │
│   │ • Generate candidate primers per target region   │              │
│   │ • Filter by Tm, length, GC%, complexity, ends    │              │
│   │ • Screen against repeats & reference genome      │              │
│   │ • Optimise for multiplex compatibility           │              │
│   └──────────────────────────┬───────────────────────┘              │
│                              │                                      │
│   4. OUTPUT                  ▼                                      │
│   ┌──────────────────────────────────────────────────┐              │
│   │ Primer report  •  Pool assignments  •  Sequences │              │
│   └──────────────────────────────────────────────────┘              │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

1. **Parse** — Reads GenBank or FASTA input; identifies target coordinates from feature annotations.
2. **Detect** — Extracts target regions (exons by default) using a priority-based recognition hierarchy.
3. **Design** — Generates candidate primers, applies thermodynamic and specificity filters, and arranges amplicons into multiplex-compatible pools.
4. **Report** — Outputs primer sequences, amplicon coordinates, pool assignments, and validation metrics.

---

## Requirements

| Requirement | Details |
|---|---|
| **Java** | Version **25** or higher ([download](https://www.oracle.com/java/technologies/downloads/)) |
| **OS** | Windows, macOS, or Linux — platform-independent |
| **RAM** | Default heap is sufficient for most genes; large genomes may need 16–128 GB (see [Memory Configuration](#memory-configuration)) |

> **Verify your Java version:**
> ```bash
> java -version
> ```
> If you see a version below 25, update Java or use Conda (see below).

---

## Installation

### Option A — Direct Download

```bash
# 1. Clone the repository
git clone https://github.com/<username>/PCRpanel.git
cd PCRpanel

# 2. Confirm Java 25+ is on your PATH
java -version

# 3. Run (no build step required — the JAR is pre-built)
java -jar dist/PCRpanel.jar test/config.file
```

### Option B — Install Java 25 via Conda

If you don't have Java 25 installed system-wide:

```bash
# Create a dedicated environment
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n java25 openjdk=25
conda activate java25

# Verify
java -version
```

---

## Quick Start

### 1. Prepare a configuration file

Create a plain-text file (e.g., `my_panel.conf`) specifying your targets and primer parameters:

```ini
# --- Input ---
target_path=/data/genes/NG_013019.gb
target_path=/data/genes/NG_011731.gb

# --- Panel mode ---
multiplex=true
homology=false

# --- Amplicon size ---
minPCR=250
maxPCR=500

# --- Primer constraints ---
minLen=18
maxLen=24
minTm=60
maxTm=62
minLC=80

# --- End constraints ---
3end=w
5end=

# --- Adapters (Illumina) ---
forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

### 2. Run PCRpanel

```bash
java -jar dist/PCRpanel.jar my_panel.conf
```

### 3. Inspect results

Primer reports and pool assignments are written alongside your input files (or to `folder_out` if specified). See [Output Files](#output-files) for details.

---

## Memory Configuration

When using a **reference genome** via `genome_path`, allocate additional heap memory with JVM flags. The table below provides recommended settings:

| Genome Size | Recommended RAM | JVM Flags |
|---|---|---|
| < 100 MB | Default | None needed |
| 100–500 MB | 16–32 GB | `-Xms8g -Xmx16g` |
| 500 MB – 2 GB | 64 GB | `-Xms32g -Xmx64g` |
| > 2 GB | 128+ GB | `-Xms64g -Xmx128g` |

**Example — human genome:**

```bash
java -Xms32g -Xmx128g -jar PCRpanel.jar config.file
```

**What the flags mean:**

- **`-Xms`** — Initial heap size. Pre-allocating avoids costly resizing during execution.
- **`-Xmx`** — Maximum heap size. Prevents `OutOfMemoryError` for large genomes.

> **Tip:** Set `-Xms` to roughly one-quarter of `-Xmx` for a good balance between startup speed and memory efficiency.

---

## Configuration Reference

All parameters are specified in a plain-text configuration file. Lines beginning with `#` are comments.

### Primer & Amplicon Parameters

| Parameter | Description | Default |
|---|---|---|
| `minPCR` | Minimum amplicon size (bp) | `250` |
| `maxPCR` | Maximum amplicon size (bp) | `500` |
| `minLen` | Minimum primer length (nt) | `18` |
| `maxLen` | Maximum primer length (nt) | `24` |
| `minTm` | Minimum melting temperature (°C) | `60` |
| `maxTm` | Maximum melting temperature (°C) | `62` |
| `minLC` | Minimum linguistic complexity (%) | `80` |
| `3end` | 3′ end constraint (e.g., `w` for weak bases A/T) | `w` |
| `5end` | 5′ end constraint | *(none)* |
| `forwardtail` | 5′ adapter tail for forward primers | Illumina P5 |
| `reversetail` | 5′ adapter tail for reverse primers | Illumina P7 |

### Panel Mode

| Parameter | Description | Default |
|---|---|---|
| `multiplex` | Generate two overlapping, multiplex-compatible amplicon pools | `true` |
| `homology` | Design consensus primers from shared sequences across input files | `false` |

### Input / Output Paths

| Parameter | Description | Repeatable? |
|---|---|---|
| `target_path` | Path to an individual GenBank or FASTA target file | Yes |
| `target_primers` | Path to an existing primer/probe list to incorporate | No |
| `folder_path` | Directory of target files (subdirectories included) | No |
| `folder_out` | Output directory for results | No |
| `genome_path` | Directory of reference genome FASTA files (subdirectories included) | No |

> **Note:** When both `target_path` and `folder_path` are provided, PCRpanel processes the **union** of all listed files plus all files discovered in the folder.

### Full Configuration Examples

<details>
<summary><strong>Windows example</strong></summary>

```ini
# Target sequences
target_path=C:\PCRpanel\test\NG_013019.gb
target_path=C:\PCRpanel\test\NG_011731.gb
target_path=C:\PCRpanel\test\NG_008847.gb
target_path=C:\PCRpanel\test\NC_000002.gb

# Existing primers to incorporate
target_primers=C:\PCRpanel\test\primers.txt

# Panel mode
homology=false
multiplex=true

# Amplicon size
minPCR=250
maxPCR=500

# Primer parameters
minLen=18
maxLen=24
minTm=60
maxTm=62

# End constraints
3end=w
5end=

# Illumina adapters
forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

</details>

<details>
<summary><strong>Linux / macOS example (batch processing)</strong></summary>

```ini
# Process all files in a directory
folder_path=/data/genes/
folder_out=/data/report/

# Reference genome for specificity checks
genome_path=/data/t2t/

# Existing primers
target_primers=/data/primers/primers.txt

# Panel mode
homology=false
multiplex=true

# Amplicon size
minPCR=250
maxPCR=500

# Primer parameters
minLen=18
maxLen=24
minTm=60
maxTm=62

# End constraints
3end=w
5end=

# Illumina adapters
forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

</details>

### Output Behaviour

| Configuration | Output Location |
|---|---|
| `folder_out` specified | The output directory is created, its contents are erased, and the results are written to the folder named 'folder_out'. |
| `folder_out` omitted | An output directory is created automatically |

---

## Input Formats

### GenBank (recommended)

**Extensions:** `*.gb`, `*.gbff`

PCRpanel parses RefSeqGene / RefSeq records (typically downloaded with `rettype=gbwithparts`). Each file should contain:

- A **FEATURES** table with gene annotations (exons, CDS, mRNA, etc.)
- An **ORIGIN** section with the nucleotide sequence

> **Fallback:** If a GenBank file lacks an ORIGIN section, PCRpanel automatically falls back to FASTA parsing.

### FASTA

Standard FASTA format with one or more records per file. Each record is treated as a **separate, full-length target region** — exon-level targeting requires GenBank annotations.

### Primer List (`target_primers`)

A plain-text file containing existing primers or probes to incorporate into the panel design. PCRpanel will use these as fixed anchors and design additional primers to fill coverage gaps.

---

## Target Region Detection

PCRpanel designs primers **only within identified target coordinates**. Targets are derived automatically from GenBank annotations or, in the absence of annotations, from the full sequence.

### Exon Recognition Priority

Many RefSeq/GenBank records lack explicit `exon` features. PCRpanel applies a tiered recognition hierarchy:

| Priority | Feature Type | Status | Notes |
|---|---|---|---|
| 1 | `exon` | ✅ Active | Used directly when present |
| 2 | `mRNA`, `ncRNA`, `rRNA`, `tRNA` | ⏸ Suspended | `join(…)` blocks would be extracted as exon coordinates |
| 3 | `CDS` | ⏸ Suspended | `join(…)` blocks used when transcript features are absent |
| 4 | *Full-sequence fallback* | ✅ Active | Entire sequence treated as one contiguous target |

> **Current behaviour:** Only priority 1 (explicit `exon`) and priority 4 (full-sequence fallback) are active. Priorities 2 and 3 are reserved for future activation.

### GenBank Coordinate Conventions

Coordinates are **1-based, inclusive** — for example, `5049..5095` spans 47 bp and includes both endpoints.

**Supported location qualifiers:**

| Qualifier | Behaviour |
|---|---|
| `join(…)` | Multi-exon features — each interval is extracted as a separate target |
| `order(…)` | Treated identically to `join` |
| `complement(…)` | Reverse strand; coordinates are extracted normally |
| `<` / `>` (partial) | Partial-bound symbols are stripped; numeric positions are used as-is |

**Examples from a FEATURES table:**

```
# Multi-exon mRNA
mRNA    join(5049..5095,26596..26643,30691..30834,46554..46688)

# Reverse-strand CDS
CDS     complement(join(330..403,1050..1120))
```

### Internal Coordinate Conversion

GenBank 1-based inclusive intervals are converted to 0-based half-open coordinates for internal processing:

```
GenBank:  a..b        (1-based, inclusive)
Internal: [a−1, b)    (0-based, exclusive end)

Example:  5049..5095  →  [5048, 5095)    length = 47 bp
```

---

## Output Files

PCRpanel generates the following output for each target:

| File | Contents |
|---|---|
| **Primer report** | Forward and reverse primer sequences, Tm, GC%, length, linguistic complexity, and amplicon coordinates |
| **Pool assignments** | Multiplex pool allocation (Pool A / Pool B) for each primer pair, ensuring no overlapping amplicons share a pool |
| **Tailed primers** | Full primer sequences including 5′ adapter tails, ready for ordering |
| **Coverage summary** | Target region coverage statistics and any gaps in amplicon tiling |

---

## Use Cases

### 1. Gene Panel Design

Design primers targeting specific exons across multiple genes — ideal for hereditary disease screening, cancer hotspot panels, and pharmacogenomics assays.

```ini
target_path=/data/BRCA1/NG_005905.gb
target_path=/data/BRCA2/NG_012772.gb
target_path=/data/TP53/NG_017013.gb

multiplex=true
minPCR=250
maxPCR=500
```

### 2. Whole-Genome Tiling (Viral Sequencing)

For complete genome coverage — for example, SARS-CoV-2 surveillance — the entire genome becomes the target. PCRpanel tiles amplicons across the full sequence without requiring exon annotations.

```ini
target_path=/data/NC_045512.2.gb
multiplex=true
minPCR=400
maxPCR=500
```

### 3. Consensus Primer Design (Homology Mode)

When you have multiple related sequences (e.g., alleles, strains, or gene family members), enable homology mode to design primers from conserved regions:

```ini
target_path=/data/strain_A.fa
target_path=/data/strain_B.fa
target_path=/data/strain_C.fa

homology=true
multiplex=false
```

### 4. Building on Existing Primers

Start from a validated primer set and let PCRpanel fill coverage gaps:

```ini
target_path=/data/target_gene.gb
target_primers=/data/validated_primers.txt

multiplex=true
minPCR=250
maxPCR=500
```

---

## Related Tools

| Tool | Description |
|---|---|
| **[NCBI RefSeq GenBank Downloader](https://github.com/rkalendar/genbanktools)** | Batch-download GenBank records by gene abbreviation or accession number — useful for assembling large target lists |

---

## Troubleshooting

### No exons or targets detected

- **Check annotations:** Verify the GenBank record contains a **FEATURES** table with `exon` entries and an **ORIGIN** section.
- **Check join syntax:** Multi-exon genes must use `join(…)` syntax in their feature locations.
- **Check encoding:** Ensure the file is saved as **UTF-8**.
- **Try FASTA:** As a fallback, convert to FASTA — the full sequence will be used as a single target.

### FASTA fallback

When FASTA input is provided, each record is treated as one full-length target region. For exon-level targeting, supply GenBank-formatted input with gene annotations.

### Java version errors

```bash
# Check version
java -version

# If multiple versions are installed, set JAVA_HOME explicitly
export JAVA_HOME=/path/to/jdk-25
export PATH="$JAVA_HOME/bin:$PATH"
```

### `OutOfMemoryError` with large genomes

Increase heap allocation:

```bash
java -Xms32g -Xmx128g -jar PCRpanel.jar config.file
```

See [Memory Configuration](#memory-configuration) for genome-size-specific recommendations.

### Primers not found for some targets

- **Widen Tm range:** Try `minTm=58` and `maxTm=64`.
- **Increase primer length range:** Try `minLen=17` and `maxLen=28`.
- **Increase amplicon size range:** Allow larger amplicons with a higher `maxPCR`.
- **Lower complexity threshold:** Reduce `minLC` (e.g., to `70`) for AT-rich or repetitive regions.

---

## FAQ

**Q: Can PCRpanel design primers for non-human genomes?**
Yes. PCRpanel is organism-agnostic. Any GenBank or FASTA sequence can be used as input — viral, bacterial, plant, or animal.

**Q: What is linguistic complexity, and why filter on it?**
Linguistic complexity measures sequence diversity on a 0–100% scale. Low-complexity regions (e.g., `AAAAAAA` or `ATATATATAT`) make poor primer binding sites. The default threshold of 80% filters out these regions.

**Q: Can I use custom adapters instead of Illumina?**
Yes. Set `forwardtail` and `reversetail` to any nucleotide sequence. Leave them blank for no tails.

**Q: How does multiplex pooling work?**
When `multiplex=true`, PCRpanel splits primer pairs into **two pools** (A and B) such that amplicons within each pool do not overlap. This enables two-reaction multiplex PCR with complete target coverage.

**Q: Can I run PCRpanel without a reference genome?**
Yes. The `genome_path` parameter is optional. Without it, PCRpanel skips the genome-wide specificity check — useful for rapid panel prototyping or when working with organisms that lack a reference assembly.

---

## Citation

If you use PCRpanel in your research, please cite:

> Kalendar, R. (2025). PCRpanel: Custom Amplicon Panel Designer. Available at <https://primerdigital.com/tools/panel.html>



---

## License

This project is distributed under the terms of the [GNU GENERAL PUBLIC LICENSE](LICENSE.txt).

---

<p align="center">
  <em>PCRpanel — Designed for researchers, by researchers.</em>
</p>
