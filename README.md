# Integrated Multi-Omics Workflow

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-red)](#citation)
[![R](https://img.shields.io/badge/R-4.x-blue)](https://www.r-project.org/)
[![Shell](https://img.shields.io/badge/Shell-Bash-lightgrey)](https://www.gnu.org/software/bash/)

---

## Paper

**Multi-Omics Integration Across Regulatory Layers Reveals Candidate Drivers of Severe COVID-19 Immunopathology**
Namra Basharat, Ghania, Hania Fahad
School of Interdisciplinary Engineering and Sciences (SINES), National University of Sciences and Technology (NUST), Islamabad, Pakistan
*Preprint, bioRxiv, 2026.* DOI: *to be added upon posting*

## Table of Contents

- [Overview](#overview)
- [Central Question](#central-question)
- [Pipeline Architecture](#pipeline-architecture)
- [Repository Structure](#repository-structure)
- [Dataset](#dataset)
- [Results by Stage](#results-by-stage)
- [Candidate Genes](#candidate-genes)
- [Tools and Versions](#tools-and-versions)
- [How to Reproduce](#how-to-reproduce)
- [Data Availability](#data-availability)
- [Authors and Contributors](#authors-and-contributors)
- [Citation](#citation)
- [License](#license)

---

## Overview

Most people who contract SARS-CoV-2 recover without incident; a smaller fraction develop life-threatening disease. This repository investigates the regulatory basis of that divergence by integrating four independent molecular layers on a 10-sample PBMC RNA-seq subset (5 healthy controls, 5 severe COVID-19 patients) drawn from [GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418), originally profiled by Arunachalam et al. (2020, *Science*).

**A note on interpretation:** the ChIP-seq layer in this study profiles **CTCF**, a constitutive chromatin-architectural protein — not an inflammatory transcriptional activator such as NF-κB. Overlap between a CTCF peak and a differentially expressed gene indicates that gene sits within an architecturally organized chromatin domain; it is not, by itself, evidence that CTCF drives the observed expression change. This distinction is discussed explicitly in the paper's Discussion and Limitations sections and is preserved throughout this repository.

## Central Question

The project addresses this at four regulatory levels simultaneously:

| Level | Approach | Question Asked |
|---|---|---|
| Transcriptome | RNA-seq + DESeq2 | Which genes are dysregulated in severe COVID-19 PBMCs? |
| Chromatin architecture | CTCF ChIP-seq (ENCODE) | Which genes sit within CTCF-organized chromatin domains? |
| Post-transcriptional | miRNA target-network analysis | Which microRNAs regulate the most dysregulated genes? |
| Genetic risk | GWAS reanalysis (Ellinghaus et al.) | Which inherited variants predispose to severe disease? |

## Pipeline Architecture

                 Raw SRA Reads (NCBI GEO: GSE152418)
                               │
                               ▼
                   ┌───────────────────────┐
                   │        Stage 1        │  FastQC → Trim Galore → HISAT2 → featureCounts
                   │  01_rnaseq_processing  │  Output: 42,355-gene × 10-sample count matrix
                   └───────────┬───────────┘
                               │
                               ▼
                   ┌───────────────────────┐
                   │        Stage 2        │  DESeq2 v1.52.0 (negative binomial, BH FDR)
                   │02_differential_expression│ Output: 1,459 significant DEGs (padj < 0.01)
                   └───────────┬───────────┘
                               │
    ┌──────────────────┬──────┴───────────┬──────────────────┐
    ▼                  ▼                  ▼                  ▼
    ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌─────────────┐
│ Stage 3 │ │ Stage 4 │ │ Stage 5 │ │ Stage 6 │
│ GO/KEGG │ │ ncRNA │ │ CTCF │ │ GWAS │
│03_functional│ │04_mirna_ │ │05_ctcf_ │ │06_gwas_ │
│_enrichment │ │network │ │chipseq │ │analysis │
└──────┬──────┘ └──────┬──────┘ └──────┬──────┘ └──────┬──────┘
│ │ │ │
└─────────────────┴─────────┬────────┴──────────────────┘
│
▼
┌───────────────────────┐
│ Stage 7 │ Multi-omics integration → Candidate gene scoring
│07_multiomics_integration│ Output: MMP8, IL23R, CD209 (triple-layer evidence)
└───────────────────────┘


## Repository Structure

Integrated-Multiomics-Workflow/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitattributes
├── .gitignore
│
├── 01_rnaseq_processing/ # SRA download → QC → alignment → quantification
│ ├── 0_setup.sh
│ ├── 1_download_and_merge.sh
│ ├── 2_reference_and_index.sh
│ ├── 3_fastqc.sh
│ ├── 4_trimming.sh
│ ├── 5_mapping.sh
│ ├── 6_name_sort.sh
│ ├── 7_featurecounts.sh
│ ├── 8_cleanup.R
│ ├── data/trimmed/ # Trim Galore QC reports
│ ├── logs/ # HISAT2 per-sample alignment logs
│ ├── reference/ # Genome/annotation (not tracked — see below)
│ └── Results/
│ ├── Alignments/ # Sorted BAM index files
│ ├── Counts/ # Raw and cleaned featureCounts matrices
│ └── fastqc/ # Per-sample FastQC reports
│
├── 02_differential_expression/ # DESeq2 differential expression
│ ├── deg_analysis.R
│ ├── DEG_results_all.csv
│ ├── DEG_results_significant.csv
│ ├── Top20_Upregulated.csv
│ ├── Top20_Downregulated.csv
│ ├── Volcano_plot.png
│ └── heatmap.png
│
├── 03_functional_enrichment/ # GO / KEGG enrichment (enrichR)
│ ├── functional_enrichment.R
│ ├── GO_BP_results.csv / GO_BP_top20.csv
│ ├── GO_CC_results.csv / GO_MF_results.csv
│ ├── KEGG_results.csv / KEGG_top20.csv
│ └── GO_BP_barplot.png / GO_dotplot.png / KEGG_barplot.png / KEGG_dotplot.png
│
├── 04_mirna_network/ # miRDB miRNA target-network analysis
│ └── mirna_interactions_complete.xlsx
│
├── 05_ctcf_chipseq/ # CTCF ChIP-seq peak annotation (ChIPseeker)
│ ├── ctcf_chipseq_annotation.R
│ ├── CTCF_annotated_peaks.csv
│ └── CTCF_anno_pie.pdf
│
├── 06_gwas_analysis/ # GWAS summary-statistics reanalysis
│ ├── gwas_analysis.Rmd / gwas_analysis.html
│ ├── GWAS_all_cleaned.csv
│ ├── GWAS_GW_significant_SNPs.csv
│ ├── GWAS_suggestive_SNPs.csv
│ ├── GWAS_chromosome_summary.csv
│ └── Manhattan_Plot.png / QQ_Plot.png / PCA_Plot.png
│
└── 07_multiomics_integration/ # Cross-layer integration and candidate gene scoring
├── multiomics_integration.Rmd / multiomics_integration.html
├── Candidate_Genes.csv
├── Triple_Overlap_Candidates.csv
├── GWAS_Locus_Genes.csv
├── Integrated_Evidence_Table.csv
├── hypergeometric_results.rds
├── Candidate_Highlight_Plot.png
├── Evidence_Heatmap.png
├── GWAS_ChIPseq_Overlap.png
├── eQTL_Workflow_Diagram.png
└── Venn_Diagram.png


**Note:** raw `.fastq.gz`, `.bam`, and reference genome/annotation files are not tracked in this repository due to size (~50 GB). All raw sequencing data is retrieved automatically by `01_rnaseq_processing/1_download_and_merge.sh`; the reference genome index is retrieved by `2_reference_and_index.sh`.

## Dataset

| Parameter | Value |
|---|---|
| GEO Accession | GSE152418 |
| Source publication | Arunachalam et al., 2020, *Science*, PMID 32788292 |
| Cell type | PBMC (peripheral blood mononuclear cells) |
| Platform | Illumina NovaSeq 6000, single-end 101 bp |
| Library kit | TruSeq Stranded mRNA (reverse-stranded) |
| Reference genome | GRCh38, Ensembl release 109 |
| Samples analyzed | 5 healthy controls + 5 severe COVID-19 patients |
| SRA runs used | 20 (2 lane-split runs merged per sample) |

## Results by Stage

**Stage 1 — RNA-seq processing** (`01_rnaseq_processing/`)
HISAT2 alignment rates: 94.29–94.88% (healthy controls), 91.24–92.62% (severe COVID-19). 62,710 → 42,355 genes retained after zero-count filtering (32.5% removed).

**Stage 2 — Differential expression** (`02_differential_expression/`)
24,719 genes tested (DESeq2 v1.52.0, `~condition` design, Benjamini-Hochberg FDR). **1,459 significant DEGs** at padj < 0.01. At the relaxed padj < 0.1 threshold: 2,221 nominally upregulated, 1,546 nominally downregulated.

**Stage 3 — Functional enrichment** (`03_functional_enrichment/`)
Top GO Biological Process term: mitotic spindle assembly checkpoint signaling (adj. p = 8.13×10⁻¹⁰, 16/26 genes). Top KEGG pathway: Cell cycle (adj. p = 1.23×10⁻¹⁴, 40/124 genes). Two dominant biological themes emerge: cell-cycle/mitotic-checkpoint disruption and mitochondrial oxidative-phosphorylation impairment.

**Stage 4 — ncRNA regulatory network** (`04_mirna_network/`)
Eight high-effect-size DEGs selected (|log₂FC| > 2, padj < 0.01); IGLC2 excluded (no miRDB targets returned), **seven genes carried forward**: IFI27, MMP8, IFIT3, ISG20 (upregulated) and CD209, KLRK1, IL23R (downregulated). 129 total gene-miRNA interactions (miRDB score ≥ 50, mean score 79.1 ± 10.6). The four upregulated genes were linked to 68 unique miRNAs combined (~17/gene; IFI27 alone to just 9, at high confidence). The three downregulated genes were each linked to exactly 20 unique miRNAs (the top-20 filtering cap), 60 unique miRNAs combined, with no overlap between their target sets.

**Stage 5 — CTCF ChIP-seq** (`05_ctcf_chipseq/`)
42,243 IDR-thresholded peaks (ENCODE [ENCSR000ATN](https://www.encodeproject.org/experiments/ENCSR000ATN/), file ENCFF190KNC, CD14+ monocytes) annotated to 17,623 genes via ChIPseeker. Genomic distribution: 39.4% within ±3 kb of a TSS, 33.6% intronic, 17.4% distal intergenic, 6.3% exonic, 2.8% 3′ UTR, 0.3% 5′ UTR, 0.1% downstream. COVID-19-relevant genes carrying CTCF peaks include IL6, CXCL2, NFKBIA, ACE2, and TMPRSS2.

**Stage 6 — GWAS reanalysis** (`06_gwas_analysis/`)
[GCST90000255](https://www.ebi.ac.uk/gwas/studies/GCST90000255) (Ellinghaus et al., 2020, *NEJM*), 8,328,612 autosomal variants after quality-control filtering, 3,815 individuals of European ancestry. Two genome-wide-significant loci (P < 5×10⁻⁸):
- **3p21.31**: lead SNP chr3:45,834,967 (G>GA), β = −0.5727, P = 1.15×10⁻¹⁰ (27 variants significant at this locus)
- **9q34.2 (ABO locus)**: chr9:133,263,862 (A>C), β = 0.2814, P = 4.95×10⁻⁸

Genomic inflation factor λ = 1.0169 (well-controlled, < 1.05 threshold).

**Stage 7 — Multi-omics integration** (`07_multiomics_integration/`)
Hypergeometric enrichment test (background: 24,719 DESeq2-tested genes): DEG–CTCF overlap = 787 observed vs. 626.5 expected by chance (P = 2.09×10⁻¹⁸), against a background-restricted CTCF gene set of 10,614. DEG–GWAS-locus overlap = 2 genes (SURF6, MED22) vs. 1.36 expected (P = 0.40, underpowered given a locus gene set of n = 23).

## Candidate Genes

| Gene | log₂FC | padj | Direction | Evidence Layers | Function |
|---|---|---|---|---|---|
| **MMP8** | +4.55 | 0.0079 | ↑ Up | DEG + CTCF + ncRNA | Neutrophil collagenase; matrix remodeling |
| **IL23R** | −2.33 | 0.00033 | ↓ Down | DEG + CTCF + ncRNA | IL-23/IL-17 immune signaling axis |
| **CD209** | −2.24 | 0.0015 | ↓ Down | DEG + CTCF + ncRNA | Dendritic-cell receptor; antigen presentation |
| SURF6 | −0.67 | 0.0023 | ↓ Down | DEG + CTCF + GWAS proximity | Ribosome biogenesis |
| MED22 | −0.57 | 0.0099 | ↓ Down | DEG + CTCF + GWAS proximity | Mediator complex / transcription |
| ADAMTS2 | +7.68 | 2.7×10⁻⁵ | ↑ Up | DEG + CTCF | Extracellular matrix processing |
| SDC1 | +5.88 | 1.3×10⁻¹⁷ | ↑ Up | DEG + CTCF | Cell adhesion / signaling |

MMP8, IL23R, and CD209 are the three genes supported by all three of the differential-expression, CTCF-occupancy, and ncRNA-targeting layers — the study's primary candidates for future experimental follow-up. SURF6 and MED22 are supported by DEG, CTCF, and proximity to the 9q34.2 GWAS locus, but this proximity association did not reach statistical significance and should be read as descriptive rather than confirmed enrichment.

## Tools and Versions

| Category | Tool | Version |
|---|---|---|
| Data retrieval | SRA Toolkit (prefetch, fasterq-dump) | 3.1.1 |
| Quality control | FastQC | 0.12.1 |
| Trimming | Trim Galore / Cutadapt | 0.6.10 / 5.2 |
| Alignment | HISAT2 | 2.2.1 |
| Quantification | featureCounts (Subread) | 2.0.3 |
| Differential expression | DESeq2 (R/Bioconductor) | 1.52.0 |
| Enrichment | enrichR | — |
| ID conversion | biomaRt | — |
| miRNA targets | miRDB | v6.0 |
| ChIP-seq annotation | ChIPseeker | — |
| GWAS analysis | R (`data.table::fread`, `phyper`) | — |
| Visualization | EnhancedVolcano, pheatmap, ggplot2 | — |

## How to Reproduce

```bash
git clone https://github.com/Namra-5/Integrated-Multiomics-Workflow.git
cd Integrated-Multiomics-Workflow
```

Each numbered folder (`01_rnaseq_processing/` through `07_multiomics_integration/`) contains its own scripts, expected inputs, and outputs, and can generally be run in numeric order since later stages consume earlier stages' outputs. Raw `.fastq.gz` and `.bam` files are not tracked in this repository due to size; `01_rnaseq_processing` scripts download raw sequencing data automatically from NCBI SRA (~50 GB free disk space required).

## Data Availability

| Layer | Accession | Source |
|---|---|---|
| RNA-seq | [GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418) | Gene Expression Omnibus |
| CTCF ChIP-seq | [ENCSR000ATN](https://www.encodeproject.org/experiments/ENCSR000ATN/) | ENCODE (CD14+ monocytes) |
| GWAS summary statistics | [GCST90000255](https://www.ebi.ac.uk/gwas/studies/GCST90000255) | NHGRI-EBI GWAS Catalog |

## Authors and Contributors

**Authors** (manuscript byline)
- **Namra Basharat** — RNA-seq processing, alignment, quantification; manuscript structuring
- **Ghania** — GWAS reanalysis, multi-omics integration
- **Hania Fahad** — Differential expression, pathway enrichment, manuscript discussion

**Contributor**
- **Nawal Babar** — performed the miRNA-target prediction and CTCF ChIP-seq analysis presented in this study (acknowledged in the manuscript; not a listed author)

Institution: School of Interdisciplinary Engineering and Sciences (SINES), National University of Sciences and Technology (NUST), Islamabad, Pakistan — BS Bioinformatics.

## Citation

If you use this code, please cite the preprint:

> Basharat N, Ghania, Fahad H. Multi-Omics Integration Across Regulatory Layers Reveals Candidate Drivers of Severe COVID-19 Immunopathology. *bioRxiv*. 2026. doi: *to be added*

A machine-readable citation is available in [`CITATION.cff`](CITATION.cff).

## License

Code in this repository is released under the [MIT License](LICENSE). The associated manuscript is separately licensed under CC BY 4.0 on bioRxiv.

## Contact

Corresponding author: Namra Basharat — nbashrat.bsbi23sines@student.nust.edu.pk
