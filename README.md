# RNA-Seq Differential Expression Analysis Workflow (DESeq2 & Downstream)

This repository contains a set of R scripts implementing an RNA-Seq differential expression analysis workflow. It utilizes DESeq2 for the core analysis and includes scripts for downstream functional enrichment (Gene Ontology) and visualization (PCA, Volcano plots, MA plots, heatmaps).

**Developed for:** Analysis of time-course RNA-Seq data with multiple cell lines and treatments

## IMPORTANT DISCLAIMER: Data not included

**This repository contains ONLY the R scripts for the analysis pipeline.**

The raw and processed data used in the original analysis are **NOT** included in this repository as they pertain to unpublished work.

**Users wishing to run this pipeline MUST provide their own input data** in the format specified below and place it in the appropriate directories (`data/`, `data/gene_sets/`). The scripts will likely fail if the correct input data is not provided in the expected locations and format.

## Workflow overview

The analysis is performed in two main stages:

1.  **Core differential expression analysis (`01_deseq2_analysis.R`):**
    *   Loads raw count data (user-provided).
    *   Performs basic filtering and QC plots.
    *   Runs the standard DESeq2 workflow (`estimateSizeFactors`, `estimateDispersions`, `nbinomWaldTest`).
    *   Performs variance stabilization (`vst`).
    *   Generates PCA plots for QC.
    *   Conducts **all pairwise comparisons** between defined experimental groups using `lfcShrink` (ashr).
    *   Generates standard plots (MA, Volcano, p-value histograms, heatmaps) for each comparison.
    *   Performs overlap analysis using UpSet plots.
    *   Saves intermediate results (CSV per comparison, RData objects) required for the next stage.

2.  **Downstream analysis and visualization (`02_downstream_analysis.R`):**
    *   Loads intermediate results (CSVs, RData) generated by Script 01.
    *   Generates enhanced Volcano plots with gene labeling (using ggplot2 and ggrepel).
    *   Performs Gene Ontology (GO) enrichment analysis for Biological Processes (BP) using `clusterProfiler`.
    *   Generates dot plots for GO results.
    *   Generates heatmaps (using pheatmap) for significant genes within the top enriched GO terms.

## Expected directory structure

The scripts expect the following directory structure within the project root:
```
your_project_root/
├── data/ # USER MUST CREATE & POPULATE
│ ├── raw.RData # USER-PROVIDED: SummarizedExperiment object named 'x'
│ └── gene_sets/ # USER MUST CREATE & POPULATE
│ └── Genes_of_interest.txt # USER-PROVIDED: Text file with gene symbols
├── scripts/ # Contains the analysis scripts
│ ├── 01_deseq2_analysis.R
│ ├── 02_downstream_analysis.R
│ └── utils_deseq_functions.R # Helper functions for script 01
├── results/ # CREATED BY SCRIPTS: Output directory
│ ├── qc/ # QC plots (histograms, PCA, preproc)
│ ├── deseq2_object.RData # Saved 'dds' object
│ ├── deseq2_reslist.RData # Saved results list, vsd, etc.
│ ├── deseq2_sizefactors.txt
│ ├── tables/ # Aggregated result tables
│ ├── pairwise_results/ # Per-comparison CSVs and basic plots from script 1
│ ├── overlap_analysis/ # UpSet plots and related heatmaps from script 1
│ └── downstream_analysis/ # Volcano plots, GO analysis, heatmaps from script 2
├── your_project.Rproj # Optional: RStudio project file
├── README.md # This file
└── .gitignore # Tells Git which files/folders to ignore (e.g., data/, results/)
```
## Dependencies

*   **R:** Version >= 4.1.0 recommended.
*   **Bioconductor:** Version >= 3.14 recommended.
*   **R Packages:**
    *   From CRAN: `tidyverse` (includes `ggplot2`, `dplyr`, `tidyr`, etc.), `pheatmap`, `RColorBrewer`, `ggrepel`, `GGally`, `ggpointdensity`, `reshape2`, `parallel`, `matrixStats`, `scales`, `ggvenn`
    *   From Bioconductor: `DESeq2`, `SummarizedExperiment`, `edgeR`, `vsn`, `ashr`, `UpSetR`, `BiocParallel`, `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `pcaExplorer`

### Install packages using:

```R
# For CRAN packages
install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel", "GGally", "ggpointdensity", "reshape2", "matrixStats", "scales", "ggvenn"))
```

### For Bioconductor packages:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "SummarizedExperiment", "edgeR", "vsn", "ashr", "UpSetR", "BiocParallel", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "pcaExplorer"))
```

## Setup instructions
Clone Repository:
```bash
git clone https://github.com/tlboi/rnaseqanalysis.git
```

### Install dependencies: 
Ensure R and all required packages (see above) are installed.

### Create data directories: 
Manually create the data/ and data/gene_sets/ directories inside the cloned project folder.

### Add input data:
Place your raw count data file, saved as an R .RData file containing a SummarizedExperiment object named x, into the data/ directory.

### The SummarizedExperiment object x must contain:
- Raw counts accessible via counts(x).
- Sample metadata accessible via colData(x), including columns named Group (a factor defining the experimental groups for comparison), Time, CellLine, and Treatment. The Group column is essential for pairwise comparisons.
- Gene metadata accessible via rowData(x), including columns named Genetype (used to identify ERCC spike-ins if present), Symbol (for gene labeling), Description, Ensembltrans, Ensemblgene, basepairs, Chr.

### Genes of interest (optional):
Place you gene of interest (one gene symbol per line) into a text file named Genes_of_interest.txt within the data/gene_sets/ directory.

### Running the analysis
Execute the scripts sequentially from the project's root directory. It is recommended to use Rscript from the terminal or run the code interactively within an R session (e.g., in RStudio opened via the .Rproj file).

#### Run DESeq2 analysis:
```bash
Rscript scripts/01_deseq2_analysis.R
```
This script performs the core DE analysis and generates necessary intermediate files in the results/ directory.

#### Run downstream analysis:
```bash
Rscript scripts/02_downstream_analysis.R
```
This script uses the output from the first script to perform GO analysis and generate final plots.

### Script descriptions
#### scripts/01_deseq2_analysis.R:
The main script for differential expression analysis. Loads data, runs DESeq2, performs all pairwise comparisons, generates QC and basic result plots (PCA, MA, heatmaps), and saves results. Relies on utils_deseq_functions.R.

#### scripts/02_downstream_analysis.R: 
Performs downstream analysis using the results from script 01. Generates publication-style Volcano plots, runs GO enrichment analysis (clusterProfiler), and creates GO term-specific heatmaps.

#### scripts/utils_deseq_functions.R: 
Contains various helper functions (primarily for plotting and QC) used by 01_deseq2_analysis.R.

### Output description
The analysis generates various output files organized within the results/ directory:

- results/qc/: Quality control plots (count distributions, PCA, sample distances, dispersion, etc.).
- results/pairwise_results/: Contains one CSV file per comparison with full DESeq2+ashr results, plus basic MA, Volcano, p-value histograms, and significant gene heatmaps for each. Also includes count plots for specific genes of interest.
- results/tables/: Aggregated tables including combined statistics across all comparisons and normalized/stabilized counts.
- results/overlap_analysis/: UpSet plots showing overlaps between DE gene sets across comparisons, and example heatmaps of commonly regulated genes.
- results/downstream_analysis/: Contains subdirectories for each comparison, holding the enhanced Volcano plots, GO enrichment results (tables and dot plots), and GO term-specific heatmaps generated by script 02.
- results/deseq2_object.RData: The main DESeqDataSet object (dds) after running the analysis.
- results/deseq2_reslist.RData: Contains the list of results data frames (reslist_deseq), aggregated tables (restot, mytot), and the variance-stabilized data (vsd).
- results/deseq2_sizefactors.txt: Table of calculated size factors per sample.# RNA-Seq Differential Expression Analysis Workflow (DESeq2 & Downstream)

