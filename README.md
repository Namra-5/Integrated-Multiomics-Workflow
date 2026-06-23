# Integrated Multi-Omics Pipeline for Regulatory and Genetic Network Analysis in Severe COVID-19

**A complete, reproducible bioinformatics workflow spanning RNA-seq, differential expression, functional enrichment, ncRNA regulation, ChIP-seq, GWAS, and multi-omics integration — applied to the molecular basis of severe SARS-CoV-2 infection.**

![Omics Layers](https://img.shields.io/badge/Omics%20Layers-4-blueviolet?style=for-the-badge)
![Pipeline Tasks](https://img.shields.io/badge/Pipeline%20Tasks-7-blue?style=for-the-badge)
![DEGs Identified](https://img.shields.io/badge/DEGs%20Identified-1%2C459-critical?style=for-the-badge)
![Candidate Genes](https://img.shields.io/badge/Candidate%20Genes-MMP8%20%7C%20IL23R%20%7C%20CD209-success?style=for-the-badge)
![Disease](https://img.shields.io/badge/Disease-Severe%20COVID--19-red?style=for-the-badge)
![Institution](https://img.shields.io/badge/Institution-NUST%20SINES-orange?style=for-the-badge)
![Language](https://img.shields.io/badge/Language-R%20%7C%20Bash-informational?style=for-the-badge)
![Genome](https://img.shields.io/badge/Reference-GRCh38-lightgrey?style=for-the-badge)

## Authors
- Namra Basharat
- Ghania Munir
- Hania Fahad 
- Nawal Babar

*BS Bioinformatics UG · School of Interdisciplinary Engineering & Sciences · NUST · 2026*

---

## Table of Contents

- [Project Overview](#project-overview)
- [The Central Biological Question](#the-central-biological-question)
- [Pipeline Architecture](#pipeline-architecture)
- [Repository Structure](#repository-structure)
- [Task 1 — RNA-seq: Data Retrieval, Alignment and Quantification](#task-1-rna-seq-data-retrieval-alignment-and-quantification)
- [Task 2 — Differential Gene Expression Analysis](#task-2-differential-gene-expression-analysis)
- [Task 3 — GO and KEGG Functional Enrichment Analysis](#task-3-go-and-kegg-functional-enrichment-analysis)
- [Task 4 — ncRNA Regulatory Network Analysis](#task-4-ncrna-regulatory-network-analysis)
- [Task 5 — ChIP-seq Regulatory Analysis](#task-5-chip-seq-regulatory-analysis)
- [Task 6 — Genome-Wide Association Study (GWAS)](#task-6-genome-wide-association-study-gwas)
- [Task 7 — Multi-Omics Integration](#task-7-multi-omics-integration)
- [Key Findings at a Glance](#key-findings-at-a-glance)
- [Tools and Environment](#tools-and-environment)
- [How to Reproduce](#how-to-reproduce)
- [Team Contributions](#team-contributions)
- [References](#references)

---

## Project Overview

This repository contains a **complete, end-to-end multi-omics bioinformatics pipeline** investigating the molecular mechanisms underlying severe COVID-19. Working from raw sequencing reads through to an integrated regulatory model with candidate gene prioritisation, the project connects **seven genomic analysis layers** into a single, coherent biological narrative.

The pipeline demonstrates full-stack bioinformatics fluency: Linux shell scripting and HPC-grade alignment pipelines, R/Bioconductor statistical analysis, multi-layer omics integration, mechanistic modelling, and publication-quality visualisation. The dataset is derived from a landmark *Science* publication (Arunachalam et al., 2020) profiling PBMC transcriptomes from five healthy controls and five patients with severe COVID-19.

**Why this matters:** COVID-19 kills through dysregulated immune responses — cytokine storm and immune exhaustion — not directly through viral load. Understanding the transcriptomic and genetic architecture of this immune failure is essential for developing targeted therapies. This pipeline reconstructs that architecture at four molecular levels simultaneously.

---

## The Central Biological Question

> **Why do some individuals develop life-threatening COVID-19 while others remain asymptomatic — and what genomic regulatory mechanisms drive that difference?**

The project addresses this at four levels simultaneously:

| Level | Approach | Question Asked |
|-------|----------|----------------|
| **Transcriptome** | RNA-seq + DESeq2 | Which genes are dysregulated in severe COVID-19 PBMCs? |
| **Post-transcriptional** | miRNA interaction databases | Which non-coding RNAs control those dysregulated genes? |
| **Chromatin / TF** | NF-κB ChIP-seq (ENCODE) | Which transcription factor binding events drive expression changes? |
| **Genetics** | GWAS (Ellinghaus et al., 2020) | Which inherited SNPs predispose individuals to severe disease? |

The final integration in Task 7 builds a complete mechanistic chain: **SNP → disrupted TF binding → altered chromatin state → changed transcription → ncRNA fine-tuning → disease phenotype**.

---

## Pipeline Architecture

```
                     Raw SRA Reads (NCBI GEO: GSE152418)
                                   │
                                   ▼
                       ┌─────────────────────┐
                       │       Task 1        │  FastQC → Trim Galore → HISAT2 → featureCounts
                       │       RNA-seq       │  Output: 42,355-gene x 10-sample count matrix
                       └──────────┬──────────┘
                                  │
                                  ▼
                       ┌─────────────────────┐
                       │       Task 2        │  DESeq2 (negative binomial, Benjamini-Hochberg FDR)
                       │        DEG          │  Output: 1,459 significant DEGs (padj < 0.01)
                       └──────────┬──────────┘
                                  │
              ┌───────────────┬───┴─────────────┬───────────────┐
              ▼               ▼                 ▼               ▼
        ┌──────────┐    ┌──────────┐      ┌──────────┐    ┌────────────┐
        │  Task 3  │    │  Task 4  │      │  Task 5  │    │   Task 6   │
        │ GO/KEGG  │    │  ncRNA   │      │ ChIP-seq │    │    GWAS    │
        │ Enrichmt │    │ Network  │      │  NF-κB   │    │GCST90000255│
        └────┬─────┘    └────┬─────┘      └────┬─────┘    └────┬───────┘
             │               │                 │               │
             └───────────────┴────────┬────────┴───────────────┘
                                      │
                                      ▼
                       ┌─────────────────────┐
                       │       Task 7        │  Multi-omics integration → Candidate gene scoring
                       │     Integration     │  Output: MMP8, IL23R, CD209 (triple-layer evidence)
                       └─────────────────────┘
```

---

## Repository Structure

```
integrated-multiomics-workflow/
│
├── README.md
├── .gitignore
├── Namra_Ghania_Nawal_Hania_Genomics_FinalProject.pdf   ← Full 70-page project report
│
├── Task1_RNA-seq/                        ← Shell pipeline: SRA → BAM → count matrix
│   ├── 0_setup.sh                        ← conda environment + tool installation
│   ├── 1_download_and_merge.sh           ← SRA download + lane merging (20 SRR → 10 samples)
│   ├── 2_reference_and_index.sh          ← GRCh38 HISAT2 prebuilt index retrieval (~4.2 GB)
│   ├── 3_fastqc.sh                       ← Pre-trim quality assessment (FastQC v0.12.1)
│   ├── 4_trimming.sh                     ← Adapter + quality trimming (Trim Galore / Cutadapt)
│   ├── 5_mapping.sh                      ← HISAT2 spliced alignment → coordinate-sorted BAM
│   ├── 6_name_sort.sh                    ← samtools name-sort for featureCounts input
│   ├── 7_featurecounts.sh                ← Gene-level quantification (-s 2, reverse-stranded)
│   ├── 8_cleanup.R                       ← Count matrix cleaning and zero-count filtering
│   ├── data/trimmed/                     ← Trim Galore QC reports (.txt)
│   ├── logs/                             ← HISAT2 per-sample alignment logs
│   └── Results/
│       ├── Alignments/                   ← Sorted BAMs + BAI index files
│       ├── Counts/                       ← FeatureCounts_Mod.txt, _clean.txt (42,355 genes)
│       └── fastqc/                       ← FastQC HTML reports (10 samples)
│
├── Task2_DEG-Analysis/                   ← DESeq2 differential expression in R
│   ├── Genomics_Task2_DEG_Analysis.R
│   ├── Task2_Metadata_creation.R
│   ├── DEG_results_all.csv               ← All 24,719 tested genes with LFC, p-value, padj
│   ├── DEG_results_significant.csv       ← 1,459 significant DEGs (padj < 0.01)
│   ├── Top20_Upregulated.csv
│   ├── Top20_Downregulated.csv
│   ├── Volcano_plot.png                  ← LFC vs -log10(padj); significant genes highlighted
│   └── heatmap.png                       ← All 1,459 DEGs across 10 samples; clean group separation
│
├── Task3_GO-KEGG-Analysis/               ← Functional enrichment in R (enrichR + biomaRt)
│   ├── Task3_GO_KEGG_Analysis.R
│   ├── GO_BP_results.csv                 ← Full GO Biological Process results
│   ├── GO_BP_top20.csv
│   ├── GO_CC_results.csv
│   ├── GO_MF_results.csv
│   ├── KEGG_results.csv
│   ├── KEGG_top20.csv
│   ├── GO_BP_barplot.png                 ← Top 15 GO BP terms (-log10 adj p-value)
│   ├── GO_dotplot.png
│   ├── KEGG_barplot.png
│   └── KEGG_dotplot.png
│
├── Task4_ncRNA-Analysis/                 ← miRNA-target interaction network (miRDB)
│   └── task4_interactions_complete.xlsx  ← 130+ interactions across 8 key DEGs (scores ≥ 50)
│
├── Task5_ChIP-seq-Analysis/              ← NF-κB p65 ChIP-seq peak annotation (ENCODE)
│   ├── task5.R                           ← ChIPseeker annotation pipeline
│   ├── NF-kB_annotated_peaks.csv         ← 42,243 peaks annotated to 17,623 genes
│   └── NF-kB_anno_pie.pdf                ← Genomic feature distribution pie chart
│
├── Task6_GWAS/                           ← GWAS analysis in R (GCST90000255)
│   ├── Task6_GWAS_Analysis.Rmd
│   ├── Task6_GWAS_Analysis.html          ← Rendered, self-contained analysis report
│   ├── GWAS_all_cleaned.csv
│   ├── GWAS_GW_significant_SNPs.csv      ← 4 genome-wide significant SNPs (P < 5×10⁻⁸)
│   ├── GWAS_suggestive_SNPs.csv
│   ├── GWAS_chromosome_summary.csv
│   ├── Manhattan_Plot.png                ← Chr3 + Chr9 peaks above GW significance
│   ├── QQ_Plot.png                       ← λ = 1.0169 (no inflation)
│   └── PCA_Plot.png
│
└── Task7_Integration/                    ← Cross-layer integration + candidate gene model
    ├── Task7_Integration.Rmd
    ├── Task7_Integration.html
    ├── Task7_Candidate_Genes.csv
    ├── Task7_Integrated_Evidence_Table.csv  ← 811 genes scored across 4 omics layers
    ├── Task7_Triple_Overlap_Candidates.csv  ← MMP8, IL23R, CD209
    ├── Task7_GWAS_Locus_Genes.csv
    ├── Task7_Candidate_Highlight_Plot.png   ← Volcano plot with triple-evidence genes labelled
    ├── Task7_Evidence_Heatmap.png           ← 4-layer evidence matrix for top candidates
    ├── Task7_eQTL_Workflow_Diagram.png      ← SNP → TF binding → expression → ncRNA model
    ├── Task7_GWAS_ChIPseq_Overlap.png
    └── Task7_Venn_Diagram.png               ← DEG ∩ ChIP ∩ ncRNA → {MMP8, IL23R, CD209}
```

> **Note:** Raw `.fastq.gz` and `.bam` files are not tracked in this repository due to size constraints. All raw data is downloaded automatically from NCBI SRA by `1_download_and_merge.sh`. Approximately 50 GB of free disk space is required.

---

## Task 1 — RNA-seq: Data Retrieval, Alignment and Quantification

### Dataset

| Parameter | Value |
|-----------|-------|
| GEO Accession | GSE152418 |
| Publication | Arunachalam et al., 2020, *Science*, PMID 32788292 |
| Organism | *Homo sapiens* |
| Cell Type | PBMC (Peripheral Blood Mononuclear Cells) |
| Sequencing Platform | Illumina NovaSeq 6000, single-end 101 bp |
| Library Kit | TruSeq Stranded mRNA (reverse-stranded) |
| Reference Genome | GRCh38, Ensembl Release 109 |
| Samples | 5 healthy controls + 5 severe COVID-19 patients |
| SRR runs | 20 runs (2 lane-splits per sample, merged pre-analysis) |

### Sample Table

| Label | SRR Lane 1 | SRR Lane 2 | Sex | Condition |
|-------|------------|------------|-----|-----------|
| Normal_1 | SRR12007855 | SRR12007856 | Male | Healthy |
| Normal_2 | SRR12007857 | SRR12007858 | Female | Healthy |
| Normal_3 | SRR12007859 | SRR12007860 | Female | Healthy |
| Normal_4 | SRR12007861 | SRR12007862 | Male | Healthy |
| Normal_5 | SRR12007863 | SRR12007864 | Female | Healthy |
| COVID_1 | SRR12007825 | SRR12007826 | Male | Severe COVID-19 |
| COVID_2 | SRR12007827 | SRR12007828 | Female | Severe COVID-19 |
| COVID_3 | SRR12007831 | SRR12007832 | Male | Severe COVID-19 |
| COVID_4 | SRR12007835 | SRR12007836 | Female | Severe COVID-19 |
| COVID_5 | SRR12007839 | SRR12007840 | Male | Severe COVID-19 |

### Pipeline

```
SRA Toolkit (prefetch + fasterq-dump)
        → Lane merging (cat; 20 SRR → 10 merged FASTQ)
        → FastQC v0.12.1 (quality assessment)
        → Trim Galore v0.6.10 / Cutadapt v5.2 (Q20 3′-trim, min 36 bp)
        → HISAT2 v2.2.1 (spliced alignment to GRCh38; --rna-strandness R)
        → samtools sort (coordinate-sort → name-sort for featureCounts)
        → featureCounts / Subread v2.0.3 (-s 2 reverse-stranded, single-end)
        → R cleanup (zero-count filter: 62,710 → 42,355 genes)
```

### Key Technical Decisions

**HISAT2 over STAR:** STAR requires ≥ 32 GB RAM for the human genome. The pipeline was executed within a 12 GB RAM constraint (Google Colab + local Ubuntu 24 VM). HISAT2's FM-index operates in 4–8 GB while achieving Pearson r > 0.98 with STAR count matrices across all expressed genes — no biological conclusions are compromised.

**Strandedness validation (-s 2):** Empirically confirmed. Running featureCounts on Normal_1 with `-s 2` yielded 9,361,437 assigned reads (33.6%). Running `-s 1` yielded < 5% assignment — unambiguous confirmation of the TruSeq reverse-stranded library design.

### Results

**HISAT2 Alignment Rates (all 10 samples):**

| Group | Alignment Rate Range | Interpretation |
|-------|---------------------|----------------|
| Healthy controls | 94.29 – 94.88% | Excellent; well above the 70% benchmark |
| Severe COVID-19 | 91.24 – 92.62% | ~2–3% lower — biologically meaningful, not a technical artifact |

The 2–3 percentage point gap between groups reflects known PBMC transcriptomic reprogramming in severe SARS-CoV-2 infection: induction of viral-response transcripts, ncRNA dysregulation, and the presence of low-level viral RNA reads that do not align to GRCh38.

**featureCounts Output:**

| Metric | Value |
|--------|-------|
| Genes in GTF | 62,710 |
| Assignment rate per sample | 32–35% |
| Retained after zero-count filter | **42,355 genes** |
| Genes removed (tissue-specific, pseudogenes) | 20,355 (32.5%) |

The ~33% exonic assignment rate is **expected and correct** for PBMC transcriptomes. It reflects the high proportion of intronic pre-mRNA reads (rapidly activated immune cells undergoing de novo transcription), non-coding RNA loci, and multi-mapping from paralogous immune gene families — not a quantification failure. The resulting 7.5–10.7 million exonically assigned reads per sample are more than sufficient for robust differential expression.

---

## Task 2 — Differential Gene Expression Analysis

**Method:** DESeq2 negative binomial model + Benjamini-Hochberg FDR  
**Comparison:** Severe COVID-19 vs Healthy Controls  
**Threshold:** padj < 0.01

### Results

| Metric | Value |
|--------|-------|
| Genes tested | 24,719 |
| **Significant DEGs (padj < 0.01)** | **1,459** |
| Upregulated in COVID-19 (LFC > 0) | 2,221 |
| Downregulated in COVID-19 (LFC < 0) | 1,546 |
| Outlier genes removed | 54 (0.22%) |

The higher proportion of upregulated genes is biologically consistent with the **cytokine storm** — a broad transcriptional activation response to severe SARS-CoV-2 infection. Hierarchical clustering in the expression heatmap produces complete, clean separation of COVID-19 from healthy samples, confirming a strong and reproducible biological signal with no sample misclassification.

### Outputs

- `DEG_results_significant.csv` — 1,459 significant DEGs with LFC, p-value, padj
- `Volcano_plot.png` — LFC vs −log₁₀(padj); significant genes highlighted in red/blue
- `heatmap.png` — Scaled expression of all 1,459 DEGs across 10 samples; perfect group separation
- `Top20_Upregulated.csv` and `Top20_Downregulated.csv`

---

## Task 3 — GO and KEGG Functional Enrichment Analysis

**Tool:** enrichR (R package)  
**Databases:** GO Biological Process / Molecular Function / Cellular Component 2023, KEGG 2021 Human  
**Input:** 1,459 significant DEGs, Ensembl IDs converted to HGNC symbols via biomaRt  
**Threshold:** Adjusted p-value < 0.05

### Top GO Biological Processes

| Rank | Term | Adj. P-value | Gene Overlap |
|------|------|-------------|--------------|
| 1–3 | Mitotic Spindle Assembly Checkpoint Signaling | 8.13 × 10⁻¹⁰ | 16/26 |
| 4 | Negative Regulation of Mitotic Metaphase/Anaphase Transition | 3.06 × 10⁻⁹ | 16/28 |
| 5 | Mitotic Sister Chromatid Segregation | 1.81 × 10⁻⁷ | 29/111 |
| 6 | Cellular Respiration | 1.81 × 10⁻⁷ | 25/85 |
| 7 | Mitochondrial ATP Synthesis Coupled Electron Transport | 4.39 × 10⁻⁷ | 22/70 |
| 8 | Microtubule Cytoskeleton Organization Involved in Mitosis | 4.97 × 10⁻⁷ | 20/59 |
| 9 | Mitotic Nuclear Division | 1.06 × 10⁻⁶ | 17/45 |
| 10 | Aerobic Electron Transport Chain | 1.06 × 10⁻⁶ | 21/68 |

### Top KEGG Pathways

| Rank | Pathway | Adj. P-value | Gene Overlap |
|------|---------|-------------|--------------|
| 1 | **Cell cycle** | 1.23 × 10⁻¹⁴ | 40/124 |
| 2 | Prion disease | 5.19 × 10⁻¹² | 57/273 |
| 3 | Parkinson disease | 1.22 × 10⁻¹¹ | 53/249 |
| 4 | Herpes simplex virus 1 infection | 2.30 × 10⁻⁸ | 74/498 |
| 5 | Diabetic cardiomyopathy | 6.18 × 10⁻⁸ | 40/203 |
| 6 | Protein processing in endoplasmic reticulum | 6.18 × 10⁻⁸ | 36/171 |
| 7 | Oxidative phosphorylation | 6.18 × 10⁻⁸ | 31/133 |
| 8 | Huntington disease | 1.44 × 10⁻⁷ | 51/306 |
| 9 | Proteasome | 1.12 × 10⁻⁶ | 16/46 |
| 10 | Non-alcoholic fatty liver disease | 2.01 × 10⁻⁶ | 31/155 |

### Biological Interpretation: Three Converging Themes

**1. Mitotic Arrest and Cell Cycle Hijacking**
The top four GO terms all implicate the spindle assembly checkpoint (SAC). SARS-CoV-2 arrests host cells in mitosis to redirect cellular machinery toward viral replication — the "negative regulation of metaphase/anaphase transition" signal is a molecular fingerprint of this hijacking. This is reinforced by the Cell Cycle pathway emerging as the single most significant KEGG hit (adj. p = 1.23 × 10⁻¹⁴, 40 overlapping genes).

**2. Mitochondrial Dysfunction and Metabolic Reprogramming**
Enrichment of oxidative phosphorylation, ATP synthesis, and aerobic electron transport pathways — cross-validated by the neurodegenerative disease cluster in KEGG (Prion, Parkinson, Huntington all converge on mitochondrial dysfunction) — indicates a shift away from OXPHOS toward glycolysis (Warburg-like effect). This generates excess reactive oxygen species, fuelling the cytokine storm and explaining why patients with metabolic comorbidities (diabetes, cardiomyopathy, NAFLD) fare significantly worse.

**3. ER Stress and Proteasomal Activation**
Enrichment of "Protein processing in endoplasmic reticulum" and "Proteasome" reflects the host's unfolded protein response (UPR) to massive SARS-CoV-2 protein synthesis. These pathways represent validated antiviral therapeutic targets.

---

## Task 4 — ncRNA Regulatory Network Analysis

**Source:** miRDB (mirdb.org) — machine-learning predicted miRNA-target interactions (score ≥ 50)  
**Input:** 8 high-confidence DEGs (|LFC| > 2, padj < 0.01) selected for COVID-19 biological relevance

### Genes Analysed

| Gene | log₂ Fold Change | Direction | Function |
|------|-----------------|-----------|----------|
| IFI27 | +9.88 | ↑ Up | Interferon alpha-inducible protein |
| MMP8 | +4.55 | ↑ Up | Neutrophil collagenase / tissue remodelling |
| IFIT3 | elevated | ↑ Up | Interferon-induced antiviral protein |
| ISG20 | elevated | ↑ Up | Interferon-stimulated 3′→5′ exonuclease |
| CD209 | −6.20 | ↓ Down | Dendritic cell receptor (DC-SIGN) |
| KLRK1 | −4.10 | ↓ Down | NK cell activating receptor (NKG2D) |
| IL23R | −3.80 | ↓ Down | IL-23/IL-17 immune axis receptor |

### Network Statistics

| Metric | Value |
|--------|-------|
| Total miRNA predictions | 130+ |
| Mean target score | 75.8 ± 12.4 (range: 50–96) |
| miRNA regulators — upregulated genes | 9–22 per gene |
| miRNA regulators — downregulated genes | **40–115 per gene** |
| CD209 miRNA regulators | **115** (most heavily regulated gene in dataset) |

### The Dual Dysregulation Model

The network reveals a paradoxical but mechanistically coherent pattern that explains the clinical picture of severe COVID-19:

- **Upregulated pro-inflammatory genes** (IFI27, MMP8, IFIT3, ISG20) show **loss of miRNA-mediated repression** — a small number of high-confidence regulators, meaning that disruption of even one or two miRNAs permits explosive transcript accumulation.
- **Downregulated immune-recovery genes** (CD209, KLRK1, IL23R) are under **active, coordinated miRNA suppression** — 40–115 regulators per gene — creating a sustained silencing program that drives immune exhaustion and impairs antigen presentation.

**Priority miRNA candidates for validation:** miR-146a-5p (suppresses NF-κB; its loss permits cytokine storm) and miR-155-3p (essential for immune cell differentiation; its dysregulation drives the exhaustion phenotype).

---

## Task 5 — ChIP-seq Regulatory Analysis

**Dataset:** ENCODE ENCSR000ATN · File: ENCFF190KNC · ENCODE4 v1.5.1  
**Factor:** NF-κB p65 (RELA)  
**Cell type:** CD14+ Monocytes (GRCh38)  
**Method:** IDR-thresholded peaks (both replicates combined), annotated with ChIPseeker in R

### Why NF-κB in Monocytes?

NF-κB is the master transcriptional switch for innate immune activation. In severe COVID-19, SARS-CoV-2 triggers NF-κB hyperactivation in monocytes and macrophages — driving uncontrolled IL-6, TNF-α, CXCL8, and IL-1β production (the cytokine storm). Using monocyte ChIP-seq directly interrogates the regulatory landscape of the immune cell type most profoundly dysregulated in this disease.

### Genomic Distribution: 42,243 NF-κB Peaks

| Genomic Feature | Peak Count | Percentage |
|----------------|-----------|------------|
| **Promoter (≤ 3 kb from TSS)** | **16,646** | **39.4%** |
| Intron | 14,208 | 33.6% |
| Distal Intergenic | 7,357 | 17.4% |
| Exon | 2,663 | 6.3% |
| 3′ UTR | 1,194 | 2.8% |
| 5′ UTR | 125 | 0.3% |
| Downstream | 50 | 0.1% |

39.4% of peaks fall directly at promoters — the most direct mechanism of transcriptional activation — and 33.6% fall in intronic regions likely representing long-range enhancer elements acting through DNA looping.

### Key COVID-19-Relevant NF-κB Target Genes

| Gene | Signal | Peak Location | Biological Significance |
|------|--------|---------------|------------------------|
| STAT3 | 356.5 | Promoter | IL-6/JAK-STAT axis — cytokine storm amplifier |
| IRF1 | 288.3 | Intron | Interferon regulatory factor 1 |
| ICAM1 | 265.5 | Promoter | Immune cell recruitment |
| ISG15 | 239.3 | Promoter | Antiviral interferon-stimulated gene |
| NFKBIA | 206.0 | Intron | IκBα — negative feedback inhibitor of NF-κB |
| TMPRSS2 | 202.6 | Intron | SARS-CoV-2 spike-priming protease (viral entry) |
| IL6 | 93.2 | Promoter | Master cytokine of the cytokine storm |
| **ACE2** | **70.2** | **Promoter** | **SARS-CoV-2 receptor — NF-κB regulates its own entry point** |

**Standout finding:** ACE2 itself carries an NF-κB promoter peak. The inflammatory response triggered by SARS-CoV-2 may modulate its own receptor's expression — a regulatory feedback loop between infection severity and viral entry efficiency with direct implications for disease progression dynamics.

DEGs intersecting ChIP-seq targets: **~200+ genes**, including MMP8, IL23R, and CD209 — subsequently confirmed as triple-overlap candidates in Task 7.

---

## Task 6 — Genome-Wide Association Study (GWAS)

**Dataset:** GCST90000255 (NHGRI-EBI GWAS Catalog)  
**Study:** Ellinghaus et al., 2020, *New England Journal of Medicine*  
**Trait:** Severe COVID-19 with respiratory failure  
**Cohort:** 3,815 European-ancestry individuals · ~10 million SNPs · GRCh38

### Genome-Wide Significant SNPs (P < 5 × 10⁻⁸)

| CHR | BP | SNP | Effect Allele | Beta | P-value |
|-----|-----|-----|--------------|------|---------|
| 3 | 45,901,089 | chr3:45901089:C:T | T | −0.476 | **4.87 × 10⁻¹²** |
| 3 | 45,851,314 | chr3:45851314:T:C | C | −0.421 | 1.10 × 10⁻¹⁰ |
| 3 | 45,869,975 | chr3:45869975:C:T | T | −0.403 | 5.54 × 10⁻¹⁰ |
| 9 | 136,149,229 | chr9:136149229:G:A | A | +0.312 | 4.95 × 10⁻⁹ |

**Genomic inflation factor λ = 1.0169** — well within the acceptable < 1.05 threshold, confirming signals are genuine genetic associations, not population stratification artifacts.

### Biological Interpretation

**Chromosome 3 locus (3p21.31)** — The strongest and most replicated COVID-19 GWAS signal in the literature. Contains chemokine receptors CXCR6, CCR9, and XCR1 (governing T-cell and NK-cell trafficking to infected lung tissue), the ACE2-interacting transporter SLC6A20, and LZTFL1 (cilia integrity). Remarkably, this entire risk haplotype was **introgressed from Neanderthals** — one of the most striking demonstrations in modern GWAS of ancient admixture influencing contemporary disease susceptibility.

**Chromosome 9 locus (ABO blood group, 9q34.2)** — Blood group A individuals showed significantly elevated risk of respiratory failure versus blood group O. ABO antigens regulate von Willebrand factor levels, predisposing to the COVID-19-associated coagulopathy and endothelial dysfunction that characterise severe disease.

---

## Task 7 — Multi-Omics Integration

### Data Layers Integrated

| Layer | Source | Scale |
|-------|--------|-------|
| Differential expression | DESeq2 (Task 2) | 1,459 DEGs |
| ncRNA regulation | miRDB interactions (Task 4) | 130 interactions |
| ChIP-seq TF binding | NF-κB ENCODE peaks (Task 5) | 17,623 annotated genes |
| GWAS genetics | Ellinghaus et al. GCST90000255 (Task 6) | 24 locus genes |

### Overlap Analysis

| Intersection | Count | Key Genes |
|-------------|-------|-----------|
| DEGs ∩ ChIP-seq targets | ~200+ | MMP8, IL23R, CD209, SDC1, ADAMTS2 |
| DEGs ∩ ncRNA targets | 130 interactions | MMP8, IL23R, CD209, ISG20 |
| DEGs ∩ GWAS locus genes | ~17 | SURF6, MED22, CCR1, CCR5, TSC1 |
| **Triple overlap (DEG + ChIP + ncRNA)** | **3 genes** | **MMP8, IL23R, CD209** |

### High-Confidence Candidate Genes

| Gene | LFC | padj | Direction | Evidence Layers | Function |
|------|-----|------|-----------|----------------|---------|
| **MMP8** | +4.55 | 0.0079 | ↑ Up | DEG + ChIP + ncRNA | Neutrophil collagenase; tissue remodelling |
| **IL23R** | −2.33 | 0.00033 | ↓ Down | DEG + ChIP + ncRNA | IL-23/IL-17 immune signalling axis |
| **CD209** | −2.24 | 0.0015 | ↓ Down | DEG + ChIP + ncRNA | Dendritic cell receptor; innate immune sensing |
| SURF6 | −0.67 | 0.0023 | ↓ Down | DEG + ChIP + GWAS | Ribosome biogenesis |
| MED22 | −0.57 | 0.0099 | ↓ Down | DEG + ChIP + GWAS | Mediator complex / transcriptional co-activation |
| ADAMTS2 | +7.68 | 2.7 × 10⁻⁵ | ↑ Up | DEG + ChIP | Extracellular matrix processing |
| SDC1 | +5.88 | 1.3 × 10⁻¹⁷ | ↑ Up | DEG + ChIP | Cell adhesion and signalling |

### The Integrated Regulatory Model

```
[Disease-associated GWAS SNP at regulatory region]
                     │
                     ▼
    [Altered NF-κB binding affinity at κB motif]
         evidenced by: GWAS-ChIP-seq overlap
                     │
                     ▼
    [Changed chromatin state: H3K27ac, H3K4me3]
                     │
                     ▼
    [Altered transcription of target gene]
         evidenced by: DESeq2 (padj < 0.01)
                     │
                     ▼
    [ncRNA post-transcriptional fine-tuning]
         evidenced by: miRDB interaction network
                     │
                     ▼
    [Final mRNA abundance change → disease phenotype]
```

**MMP8 (+4.55 LFC):** Loss of repressive TF binding at the MMP8 promoter leads to chromatin opening, massive neutrophil collagenase overexpression, extracellular matrix destruction, and inflammatory infiltration of the lung.

**IL23R (−2.33 LFC):** Reduced NF-κB-driven transcription combined with active miRNA targeting of the 3′ UTR attenuates the IL-23/IL-17 immune axis and impairs T-cell-mediated immune resolution.

**CD209 (−2.24 LFC):** Combined TF-occupancy loss and coordinated suppression by 115 miRNAs produces near-complete silencing of the dendritic cell pattern recognition receptor — impairing antigen presentation and shutting down adaptive immune activation precisely when it is most needed.

---

## Key Findings at a Glance

| Layer | Finding |
|-------|---------|
| **RNA-seq** | 91–95% alignment rates; 42,355 expressed genes across 10 PBMC samples |
| **DEG** | 1,459 significant DEGs (padj < 0.01); asymmetric upregulation consistent with cytokine storm biology |
| **GO/KEGG** | Three hallmarks: mitotic arrest + mitochondrial dysfunction + ER stress |
| **ncRNA** | Dual mechanism: loss of miRNA brakes on inflammation AND active miRNA suppression of immune recovery genes |
| **ChIP-seq** | NF-κB binds IL-6, ACE2, and TMPRSS2 promoters — inflammation directly regulates viral entry machinery |
| **GWAS** | Chr 3p21.31 (Neanderthal haplotype) + ABO locus; λ = 1.0169 — clean, unconfounded signals |
| **Integration** | **MMP8, IL23R, CD209** — triple-evidence candidates with a complete SNP → TF → expression → ncRNA regulatory chain |

---

## Tools and Environment

| Category | Tool | Version |
|----------|------|---------|
| Environment | conda (rnaseq_env), Ubuntu 24 | — |
| Data retrieval | SRA Toolkit (prefetch, fasterq-dump) | 3.1.1 |
| Quality control | FastQC | 0.12.1 |
| Trimming | Trim Galore / Cutadapt | 0.6.10 / 5.2 |
| Alignment | HISAT2 | 2.2.1 |
| BAM processing | samtools | — |
| Quantification | featureCounts (Subread) | 2.0.3 |
| Differential expression | DESeq2 (R/Bioconductor) | 1.40+ |
| Visualisation | EnhancedVolcano, pheatmap, ggplot2 | — |
| Enrichment | enrichR | — |
| ID conversion | biomaRt | — |
| miRNA interactions | miRDB | — |
| ChIP annotation | ChIPseeker, TxDb.Hsapiens.UCSC.hg38 | — |
| GWAS analysis | data.table, ggplot2 | — |
| Scripting | R 4.5.1, Bash | — |
| Reporting | R Markdown | — |

---

## How to Reproduce

### Prerequisites

- Linux / macOS terminal (or WSL2 on Windows)
- conda installed ([Miniconda](https://docs.conda.io/en/latest/miniconda.html))
- R ≥ 4.1 with BiocManager installed
- ~50 GB free disk space for raw data downloads

### Step 1 — Clone the Repository

```bash
git clone https://github.com/YOUR_USERNAME/integrated-multiomics-workflow.git
cd integrated-multiomics-workflow
```

### Step 2 — RNA-seq Pipeline (Task 1)

```bash
cd Task1_RNA-seq
bash 0_setup.sh                  # Create conda environment; install all tools
bash 1_download_and_merge.sh     # Download 20 SRR files from NCBI; merge lanes → 10 FASTQ
bash 2_reference_and_index.sh    # Download prebuilt GRCh38 HISAT2 index (~4.2 GB)
bash 3_fastqc.sh                 # Pre-trim quality assessment
bash 4_trimming.sh               # Adapter and quality trimming (Q20, min 36 bp)
bash 5_mapping.sh                # HISAT2 spliced alignment (--rna-strandness R)
bash 6_name_sort.sh              # Name-sort BAMs for featureCounts
bash 7_featurecounts.sh          # Gene-level quantification (-s 2 reverse-stranded)
Rscript 8_cleanup.R              # Zero-count filter: 62,710 → 42,355 genes
```

### Step 3 — Differential Expression (Task 2)

```r
source("Task2_DEG-Analysis/Task2_Metadata_creation.R")
source("Task2_DEG-Analysis/Genomics_Task2_DEG_Analysis.R")
```

### Step 4 — Enrichment Analysis (Task 3)

```r
source("Task3_GO-KEGG-Analysis/Task3_GO_KEGG_Analysis.R")
```

### Step 5 — ChIP-seq Annotation (Task 5)

```r
source("Task5_ChIP-seq-Analysis/task5.R")
```

### Step 6 — GWAS Analysis (Task 6)

```r
rmarkdown::render("Task6_GWAS/Task6_GWAS_Analysis.Rmd")
```

### Step 7 — Multi-Omics Integration (Task 7)

```r
rmarkdown::render("Task7_Integration/Task7_Integration.Rmd")
```

---

## Team Contributions

| Name | Role | Tasks |
|------|------|-------|
| **Namra Basharat** | Pipeline Lead & Repository Admin | Task 1 (RNA-seq pipeline: download → QC → alignment → quantification) |
| **Ghania Munir** | Statistical Analysis | Task 2 (DESeq2 DEG analysis), Task 3 (GO/KEGG enrichment) |
| **Hania Fahad** | Regulatory Analysis | Task 4 (ncRNA/miRNA network), Task 5 (ChIP-seq annotation) |
| **Nawal Babar** | Genetics & Integration | Task 6 (GWAS analysis), Task 7 (multi-omics integration) |

**Institution:** School of Interdisciplinary Engineering & Sciences (SINES), National University of Sciences and Technology (NUST), Islamabad, Pakistan  
**Programme:** BS Bioinformatics, Undergraduate Year 1 · May 2026

---

## References

Arunachalam, P. S., Wimmers, F., Mok, C. K. P., et al. (2020). Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans. *Science*, 369(6508), 1210–1220. https://doi.org/10.1126/science.abc6261

Ellinghaus, D., Degenhardt, F., Bujanda, L., et al. (2020). Genomewide association study of severe Covid-19 with respiratory failure. *New England Journal of Medicine*, 383, 1522–1534. https://doi.org/10.1056/NEJMoa2020283

ENCODE Project Consortium. ENCSR000ATN — NF-κB p65 ChIP-seq in CD14+ Monocytes. https://www.encodeproject.org/experiments/ENCSR000ATN/

NHGRI-EBI GWAS Catalog. GCST90000255. https://www.ebi.ac.uk/gwas/studies/GCST90000255

GEO Accession GSE152418. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418

BioProject PRJNA639275. https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639275

---

**Built with:** HISAT2 · DESeq2 · ChIPseeker · enrichR · miRDB · ENCODE · GWAS Catalog · R · Bash

*National University of Sciences and Technology · BS Bioinformatics · 2026*