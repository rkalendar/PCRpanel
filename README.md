# PCRpanel — Custom Amplicon Panel Designer

PCRpanel designs **custom amplicon panels** for **NGS** and **Oxford Nanopore (ONT)**, including **ultra‑high multiplex tiling PCR** and standard singleplex/multiplex PCR tasks.

**Online version:** https://primerdigital.com/tools/panel.html  
**Author:** Ruslan Kalendar  
**Contact:** ruslan.kalendar@helsinki.fi

---

## Key capabilities
- Multiplex tiling PCR and conventional PCR primer design
- Works with **any sequence length** and **any number of targets**
- Target regions can be **exons, introns, promoters, or arbitrary coordinates**
- Can design primers starting from **existing primer/probe lists**
- Supports optional **5′/3′ tails** (adapters, UMIs, etc.) during design

---

## Availability and requirements
- **OS:** Platform independent  
- **Java:** 25 or higher

Java downloads:
- https://www.oracle.com/java/technologies/downloads/

How to set/change JAVA_HOME / PATH:
- https://www.java.com/en/download/help/path.html

---

## Installing Java 25 with Conda (optional)
1) Add conda-forge and set strict priority:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

2) Create an environment and install OpenJDK 25:
```bash
conda create -n java25 openjdk=25
conda activate java25
```

3) Verify:
```bash
java -version
```

---

## Running from the command line
The executable JAR is **`PCRpanel.jar`** in the **`dist`** directory. Copy it anywhere and run:

```bash
java -jar <PCRpanelPath>/dist/PCRpanel.jar <config.file>
```

Example (Windows):
```bash
java -jar C:\PCRpanel\dist\PCRpanel.jar C:\PCRpanel\test\config.file
```

---

## Basic usage (config file)
Parameters and input paths are provided via a plain text config file:

**config.file**
```ini
target_path=C:\PCRpanel\test\NG_008690.txt
target_path=C:\PCRpanel\test\NG_011731.txt
target_path=C:\PCRpanel\test\NG_013019.txt
target_path=C:\PCRpanel\test\NG_008847.txt
target_path=C:\PCRpanel\test\NC_000002.txt

target_primers=C:\PCRpanel\test\primers.txt

minPCR=250
maxPCR=500
minLen=18
maxLen=24
minTm=60
maxTm=62
minLC=78

3end=w
5end=

forwardtail=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
reversetail=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

---



### Batch processing: per-file vs folder input, and output folder control
You can provide **individual files** (repeat `target_path=` multiple times) **or** point to a **folder** that contains multiple targets.

**Input folder:**
```ini
folder_path=C:\PCRpanel\test\
```

**Output folder:**
```ini
folder_out=C:\PCRpanel\report\
```

If `folder_out` is **not specified**, PCRpanel writes results into the **same directory as each input target file** (i.e., the directory of the `target_path=` file).  
If you use `folder_path` and `folder_out` is not specified, results are written into that **input folder**.

> When both `target_path` and `folder_path` are provided, PCRpanel should treat the final input set as the **union** of all listed files plus all files discovered in `folder_path`.


## Input formats
PCRpanel accepts:
- **GenBank flat files** (recommended): `*.gb`, `*.gbff`, RefSeqGene/RefSeq records (often downloaded with `rettype=gbwithparts`)
- **FASTA**: one or multiple records per file

If a GenBank file does not contain an `ORIGIN` sequence section (or sequence extraction fails), PCRpanel can fall back to FASTA parsing when the input is FASTA-formatted.

---

## Target regions and exon recognition for PCR design
PCRpanel designs primers **only within target coordinates**. These target coordinates can be derived automatically from GenBank annotations or specified as arbitrary regions.

### Coordinates are extracted from GenBank FEATURES
GenBank coordinate rules that matter:
- Locations are **1-based, inclusive** (e.g., `5049..5095` includes both ends).
- Locations can be wrapped by:
  - `join(...)` for multi-exon features
  - `order(...)` (treated like join for target extraction)
  - `complement(...)` (strand information; coordinates are still extracted normally)
- Partial bounds can be written as `<123..456` or `123..>456` (PCRpanel should treat `<`/`>` as the numeric position).

### What PCRpanel should recognize as “exons”
Many RefSeq/GenBank records **do not** contain explicit `exon` features. Instead, exon blocks are encoded inside transcript features.

Recommended recognition priority for exon-like blocks:
1. **`exon`** features (if present)
2. Transcript features with joined locations:
   - **`mRNA`**, **`ncRNA`**, **`rRNA`**, **`tRNA`** (use their `join(...)` blocks as exons)
3. **`CDS`** (use its `join(...)` blocks as coding exons when transcript features are absent)
4. Fallback: if nothing is annotated, use the **full sequence** as one target region

This is the most robust strategy across RefSeqGene (NG_*), genomic (NC_*), and transcript records.

### Examples of exon-bearing locations
From the FEATURES table:

**mRNA with exons**:
```
mRNA            join(5049..5095,26596..26643,30691..30834,46554..46688)
```

**CDS on reverse strand**:
```
CDS             complement(join(330..403,1050..1120))
```

The extracted target blocks are the individual intervals inside `join(...)`.

### Coordinate conversion (implementation note)
If your internal representation uses 0-based half-open coordinates (common in Java arrays),
convert GenBank intervals `a..b` to:
- `start = a - 1`
- `end   = b`   (end exclusive)

For example, `5049..5095` → `[5048, 5095)` (length 47).

---

## Whole-genome tiling panels (no exons)
Example: SARS-CoV-2 (NC_045512.2)  
https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2

For whole-genome tiling, the “target region” is simply the entire genome length. PCRpanel can design tiled amplicons across the full sequence, independent of exon annotations.

---

## Troubleshooting notes
- If no exons/targets are detected from FEATURES, verify:
  - The record contains a **FEATURES** table and an **ORIGIN** sequence
  - Transcript features use `join(...)` (common for multi-exon genes)
- If you supply FASTA, PCRpanel should treat each FASTA record as a separate target (full-length region).

