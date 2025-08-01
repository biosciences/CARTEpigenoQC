---
title: 'CARTEpigenoQC: A Quality Control Toolkit for CAR-T Single-Cell Epigenomic Data'
tags:
  - R
  - bioinformatics
  - pathogen genomics
  - reproducibility
authors:
  - name: Kaitao Lai
    affiliation: 1
    orcid: 0000-0002-9420-9352
affiliations:
  - name: The University of Sydney
    index: 1
date: "2025-07-28"
bibliography: paper.bib
---

# Summary

CARTEpigenoQC is an R-based toolkit designed to streamline quality control (QC) for single-cell epigenomic datasets involving Chimeric Antigen Receptor (CAR)-engineered T cells. With the growing application of scATAC-seq, scCUT&Tag, and scBS-seq to characterize CAR-T cell states [@satpathy2019], it has become critical to perform customized QC that not only addresses standard metrics like FRiP (Fraction of Reads in Peaks) and TSS enrichment, but also directly detects signal from CAR vector insertion sites.

CARTEpigenoQC supports both 10x Genomics and non-10x data formats and produces HTML and PNG summary outputs suited for exploratory analysis and regulatory-grade preclinical reporting. It is intended to assist researchers, core facilities, and translational immunologists in ensuring the validity of single-cell epigenomic profiling of engineered T cells [@finck2022; @robinson2023].

# Statement of Need

Despite the rise in single-cell chromatin profiling technologies in immunotherapy, current QC pipelines do not provide functionality to assess CAR vector integration in scATAC-seq or related datasets. Existing tools like ArchR [@granja2021], Signac [@stuart2021], or SnapATAC [@fang2021] focus on general chromatin metrics but lack transgene-specific support. CARTEpigenoQC fills this gap by integrating standard QC (e.g., FRiP, duplication, peak overlap) with direct inspection of vector signal such as $\text{EF1-}\alpha$, PGK1, or TRAC locus insertions that are common in CAR constructs [@eyquem2017].

This tool was inspired by the increasing use of vector-based gene therapies and the demand for preclinical quality pipelines [@maude2018]. It is lightweight, fast, and does not require complex dependencies. It generates clean, reproducible Markdown-based reports suitable for clinical handoff, regulatory review, or publication.

# Features

- FRiP calculation per cell using fragment overlaps with peak regions.
- Visualization of CAR vector site coverage across single cells.
- Markdown-based report with interactive and static plots.
- Compatible with 10x Genomics 'fragments.tsv.gz', ArchR fragments, and BED-formatted CAR insertion sites.
- Designed for reproducibility and integration into preclinical immunotherapy pipelines.

# Example Usage

CARTEpigenoQC includes a 'run_qc.R' script and report template. When supplied with fragment files and a BED file of CAR insertion coordinates, it calculates QC metrics and generates a self-contained HTML report.

# Example Outputs

## FRiP Score Distribution

This histogram visualizes the distribution of FRiP (Fraction of Reads in Peaks) across all cells. A peak around 0.2–0.4 is expected for high-quality scATAC-seq data [@chen2019]. A right-skewed distribution is typical, and cells with FRiP < 0.2 may be filtered out in downstream analyses.

![FRiP Histogram](figures/example_frip_plot.png)

## Top Cells by FRiP

This table lists the top 20 cells with the highest FRiP scores. These cells often exhibit clearer chromatin signatures and are more informative for downstream applications like motif enrichment or regulatory network analysis.

![Top FRiP Table](figures/frip_top_table.png)

## CAR Insertion Site Coverage Heatmap

The heatmap shows per-cell read coverage across CAR vector insertion sites (e.g., EF1a, PGK1, TRAC). Each row represents a cell barcode, and each column represents a CAR target site. This allows visualization of which cells are likely to carry vector insertions and where signal enrichment occurs, helping to verify construct delivery or identify mosaicism in editing [@torres2022].

![CAR Insertion Heatmap](figures/car_coverage_heatmap.png)

## Top CAR Sites by Total Coverage

This table summarizes the top 10 CAR insertion sites ranked by the total number of reads mapped to each site across all cells. It helps identify the most active or consistently covered transgene targets, which is critical in validating vector design and insertion stability [@frangieh2021].

![Top CAR Sites](figures/car_site_table.png)

# Repository and Installation

The source code is hosted on GitHub:  
https://github.com/biosciences/CARTEpigenoQC

# Acknowledgements

The author thanks collaborators at the University of Sydney for insights into CAR-T clinical pipelines, and for pipeline development guidance.

# References
