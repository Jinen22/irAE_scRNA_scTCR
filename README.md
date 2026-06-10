# Distinct Inflammatory Cytotoxic T Lymphocyte Populations Mediate PD-1 Blockade Induced Immune-Related Adverse Events in Multiple Organs

This repository contains the custom analysis code used in the study **"Distinct Inflammatory Cytotoxic T Lymphocyte Populations Mediate PD-1 Blockade Induced Immune-Related Adverse Events in Multiple Organs"** (Cancer Research, 2026).  
The scripts reproduce all major figures and analyses reported in the manuscript.

## Abstract

Immune checkpoint blockade-induced immune-related adverse events (irAEs) hamper the application of this revolutionary anti-tumor therapeutic strategy. Here, we explored the mechanisms driving irAEs by profiling the immune ecosystem of major irAE-affected organs at the single-cell scale. The analysis identified three populations of cytotoxic T lymphocytes that mediate anti-tumor immunity (CTL1) or that induce irAE in the gut (CTLirAE-I) or in multiple other organs (CTLirAE-II). Interleukin-JAK1 signaling was specifically activated in the CTLirAE-II population upon PD-1 blockade. Targeting JAK1 remarkably relieved the irAEs in the heart and lung, without compromising the anti-tumor efficacy. Tracking TCR sequence and transcriptome showed that CTLirAE-II and CTL1 populations originated from lymph node progenitor cells, while the CTLirAE-I population was derived from tissue-resident memory T cells. Moreover, irAEs could be monitored by assessing the CTLirAE-II population in circulation. In conclusion, this study elucidates the landscape of cellular changes in irAEs across multiple organs following immunotherapy and proposes strategies for relieving irAE symptoms and facilitating diagnosis.

**Significance:** Dissecting PD-1 blockade-induced immune-related adverse events across multiple organs at a single-cell scale elucidates regulators of pathogenic progression and clonal evolution, providing strategies for diagnosis and treatment without impairing anti-tumor efficacy.

## Data Availability

The processed single-cell RNA-seq data reported in this study have been deposited in the **OMIX** database (China National Center for Bioinformation) under accession number **OMIX016674** and are publicly available at [https://ngdc.cncb.ac.cn/omix](https://ngdc.cncb.ac.cn/omix).  
Raw FASTQ files generated in this study have been deposited in the **Genome Sequence Archive (GSA)** at the **National Genomics Data Center (NGDC)** under BioProject accession **PRJCA035557** (scRNA-seq: CRA040650; scTCR-seq: CRA041426) and are publicly accessible at [https://ngdc.cncb.ac.cn/gsa](https://ngdc.cncb.ac.cn/gsa).  

## Code Availability

The custom analysis code central to this work is publicly available in this GitHub repository:  
[https://github.com/YourUsername/YourRepo](https://github.com/YourUsername/YourRepo) *(replace with your actual repository URL)*


- **`Figure_2.R`** – Generates all panels of Figure 2, including UMAP visualizations colored by cell type, tissue, and treatment group; a marker gene dotplot; and log2 fold-change barplots of immune cell proportions across tissues.
- Additional scripts for other figures/tables are named accordingly and follow the same structure.

## System Requirements

The code has been tested under **R version 4.3.0** (or higher) and requires the following R packages:

| Package     | Version  | Usage                     |
|-------------|----------|---------------------------|
| Seurat      | ≥ 4.3.0  | Single-cell data handling |
| ggplot2     | ≥ 3.4.0  | Plotting                  |
| dplyr       | ≥ 1.1.0  | Data manipulation         |
| patchwork   | ≥ 1.1.0  | Combining plots           |
| ggsci       | ≥ 3.0.0  | Color palettes            |
| ggpubr      | (any)    | Theme elements            |
| reshape2    | (any)    | Data reshaping            |

You can install all required packages by running:

```r
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "ggsci", "ggpubr", "reshape2"))
```
Note: Some packages may require additional system dependencies (e.g., hdf5r, SeuratObject). Please refer to the official documentation of each package for installation instructions.


How to Run
```
Clone this repository
git clone https://github.com/YourUsername/YourRepo.git
cd YourRepo
```
Download the input data

Download the processed Seurat object (Mouse_anit-PD1_irAEs_scRNA_immune_cell_label.rds) from the OMIX database (accession OMIX008845) or directly from the source provided in the data availability statement.

Place this file into the data/ folder within the repository. Create the folder if it does not exist.

Open R or RStudio and set the working directory to the repository root. You can also open Figure_2.R directly in RStudio, which will automatically set the working directory via setwd(dirname(rstudioapi::getActiveDocumentContext()$path)).

Run the script
Execute Figure_2.R entirely or section by section. All output figures (PDFs) will be saved in the current working directory.

For other analysis scripts, follow the same procedure – ensure the required data files are present in the data/ folder before running them.

### Citation
If you use the data or code from this repository, please cite our manuscript:

Liu, Xiaowei, et al. "Distinct Inflammatory Cytotoxic T Lymphocyte Populations Mediate PD-1 Blockade Induced Immune-Related Adverse Events in Multiple Organs." Cancer Research (2026).

Contact
For questions regarding the code, please open an issue in this repository or contact the corresponding author at [email address].
For data access requests, please use the data availability statement information above.
