# Promoter Region Extraction Pipeline  
**hg38 | 500 bp upstream | Strand-aware | BEDTools**

## What This Pipeline Does

This workflow converts raw human gene annotation data into **strand-aware promoter regions** (500 bp upstream of each transcription start site, TSS) suitable for:

- Motif discovery and transcription factor scanning
- Gene Ontology (GO) enrichment workflows
- ChIP-seq / ATAC-seq overlap analysis
- Regulatory genomics pipelines

The final output is a clean BED file containing biologically correct promoter intervals.

---

## Pipeline Overview

```text
human_gene_annotation.tsv.gz
        │
        ▼
 [Parse TSS coordinates]
        │
        ▼
 [Filter valid chromosomes]
        │
        ▼
 [Generate 500 bp upstream promoters]
        │
        ▼
 promoters_500bp.bed   ✓
```

---

# Quick Start

## 1. Create Environment

```bash
mamba create -n go_enrichment python=3.12 -y
mamba activate go_enrichment

mamba install -c bioconda bedtools samtools emboss -y
```

---

## 2. Prepare Genome Index

Index the hg38 reference genome and generate a BEDTools-compatible genome file.

```bash
samtools faidx hg38.fa

cut -f1,2 hg38.fa.fai > hg38.genome
```

This creates:

| File | Purpose |
|---|---|
| `hg38.fa.fai` | FASTA index |
| `hg38.genome` | Chromosome sizes for BEDTools |

---

## 3. Parse TSS Coordinates

Extract transcription start sites from the annotation table and convert them into BED format.

```bash
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"}
NR>1{
    if($8==-1 || $8=="" || $7=="") next

    chrom="chr"$5

    if(chrom=="chrMT")
        chrom="chrM"

    tss=$8 + 0

    strand = ($6==-1) ? "-" : "+"

    print chrom,
          tss,
          tss+1,
          chrom"@"tss"-"(tss+1)"|"$7,
          ".",
          strand
}' > genes_tss_clean.bed
```

---

## 4. Remove Unsupported Chromosomes

Keep only chromosomes present in the reference genome index.

```bash
grep -Fwf <(cut -f1 hg38.genome) \
    genes_tss_clean.bed \
    > genes_tss_final.bed
```

---

## 5. Generate Promoter Windows

Extend each TSS 500 bp upstream in a strand-aware manner.

```bash
bedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed
```

---

# Input Files

| Resource | Description |
|---|---|
| `hg38.fa` | Human reference genome (UCSC hg38) |
| `human_gene_annotation.tsv.gz` | Gene annotation table containing TSS positions |

---

# Output Files

| File | Description |
|---|---|
| `genes_tss_final.bed` | Cleaned TSS coordinates |
| `promoters_500bp.bed` | Final promoter intervals for downstream analysis |

---

# BED File Format

Example:

```text
chr1    11868    12369    chr1@11868-11869|DDX11L1    .    +
```

| Column | Meaning |
|---|---|
| 1 | Chromosome |
| 2 | Start coordinate |
| 3 | End coordinate |
| 4 | Feature ID |
| 5 | Score (`.` unused) |
| 6 | Strand (`+` or `-`) |

---

# Coordinate Interpretation

For a **positive-strand gene**:

```text
chr1    11868    12369
```

- TSS = 12368
- Promoter extends 500 bp toward smaller coordinates

```text
<=====[500 bp promoter]=====[TSS]>>>>>>>>>>>
```

---

For a **negative-strand gene**:

```text
chrM    4400    4901
```

- TSS = 4400
- Upstream direction is toward larger coordinates

```text
<<<<<<<<<<<<<<[TSS]=====[500 bp promoter]=====>
```

This strand-aware behavior is controlled by:

```bash
bedtools slop -s
```

---

# Example

## Before Extension

```text
chrM    4400    4401    chrM@4400-4401|MT-TQ    .    -
```

## After `bedtools slop -l 500 -s`

```text
chrM    4400    4901    chrM@4400-4401|MT-TQ    .    -
```

---

# Downstream Applications

The generated promoter BED file can be used directly in:

## Motif Discovery

Compatible with:

- MEME Suite
- HOMER
- JASPAR motif scanning
- FIMO

Example:

```bash
bedtools getfasta \
    -fi hg38.fa \
    -bed promoters_500bp.bed \
    -s \
    -name \
    > promoters.fa
```

---

## GO Enrichment Analysis

Map promoter-associated genes to biological functions.

Typical workflow:

```text
Promoter regions
      ↓
Motif enrichment
      ↓
Gene mapping
      ↓
GO / pathway enrichment
```

---

## Peak Overlap Analysis

Intersect promoters with regulatory peaks.

Examples:

### ChIP-seq overlap

```bash
bedtools intersect \
    -a promoters_500bp.bed \
    -b chipseq_peaks.bed \
    -wa -u
```

### ATAC-seq overlap

```bash
bedtools intersect \
    -a promoters_500bp.bed \
    -b atac_peaks.bed \
    -wa -u
```

---

# Notes

- Coordinates are automatically clipped to chromosome boundaries using `hg38.genome`
- Mitochondrial chromosome names are normalized:

```text
chrMT → chrM
```

- Invalid or missing TSS entries are discarded during parsing

---

# Final Output

```text
promoters_500bp.bed
```

This is the primary file used for all downstream regulatory genomics analyses.
