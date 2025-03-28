# Differential Expression Analysis using DESeq2

# --- 1. Setup ---

# Libraries
library(DESeq2)
library(SummarizedExperiment) # For rowData, colData etc.
library(edgeR)              # For cpm()
library(vsn)                # For vst()
library(ggplot2)
library(ggrepel)
library(GGally)             # For ggpairs
library(ggpointdensity)
library(pheatmap)
library(RColorBrewer)
library(ashr)               # For LFC shrinkage
library(UpSetR)             # For upset plots
library(BiocParallel)       # For parallel processing
library(stringr)            # For string manipulation (used in original)
# library(tidyr)            # Needed by utils_deseq_functions
# library(reshape2)         # Needed by utils_deseq_functions
# library(pcaExplorer)      # Needed by utils_deseq_functions
# library(matrixStats)      # Needed by utils_deseq_functions
# library(scales)           # Needed by utils_deseq_functions


# Parameters
FDR_THRESHOLD <- 0.05 # Significance threshold for plots/heatmaps/overlaps
LFC_THRESHOLD <- 1.0  # Log2 fold change threshold for plots/heatmaps/overlaps
N_TOP_PCA <- 500      # Number of top variable genes for PCA
N_CPU <- max(1, detectCores() - 1) # Number of cores for parallel processing

# Genes of interest for plotting counts
# genes_of_interest <- c("", "", "", "", "") 

# Input/Output Files and Directories
# Assuming script is run from project root or using RStudio Project
data_dir <- "data"
results_dir <- "results"

raw_rdata_file <- file.path(data_dir, "raw.RData") # Expects a SummarizedExperiment 'x'
functions_file <- "scripts/utils_deseq_functions.R"

# Create output directories if they don't exist
dir.create(results_dir, showWarnings = FALSE)
qc_dir <- file.path(results_dir, "qc")
pairwise_dir <- file.path(results_dir, "pairwise_results")
tables_dir <- file.path(results_dir, "tables")
overlap_dir <- file.path(results_dir, "overlap_analysis")
dir.create(qc_dir, showWarnings = FALSE)
dir.create(pairwise_dir, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)
dir.create(overlap_dir, showWarnings = FALSE)

# Source utility functions
if(file.exists(functions_file)) {
  source(functions_file)
} else {
  stop("Utility functions file not found: ", functions_file)
}

# --- 2. Load Data ---

if(file.exists(raw_rdata_file)) {
  load(raw_rdata_file) # Should load object 'x' (SummarizedExperiment or similar)
  cat("Loaded data from:", raw_rdata_file, "\n")
  # Basic check of the loaded object
  if(!exists("x") || !inherits(x, "SummarizedExperiment")) {
    stop("'raw.RData' must contain a SummarizedExperiment object named 'x'.")
  }
  # Ensure rowData has 'Genetype' and colData has 'Group'
  if (!"Genetype" %in% colnames(rowData(x))) stop("'Genetype' column not found in rowData(x).")
  if (!"Group" %in% colnames(colData(x))) stop("'Group' column not found in colData(x).")
  if (!"Time" %in% colnames(colData(x))) stop("'Time' column not found in colData(x).")
  if (!"CellLine" %in% colnames(colData(x))) stop("'CellLine' column not found in colData(x).")
  if (!"Treatment" %in% colnames(colData(x))) stop("'Treatment' column not found in colData(x).")
  # Ensure Group is a factor for DESeq2 contrasts
  colData(x)$Group <- as.factor(colData(x)$Group)
  
} else {
  stop("Raw data file not found: ", raw_rdata_file)
}

cat("Dimensions of raw count data (genes x samples):", dim(x), "\n")
print(head(colData(x)))

# --- 3. Initial QC Plots (Raw Counts) ---

# Distribution of log10 total CPM per gene
pdf(file = file.path(qc_dir, "hist_log10_cpm_total.pdf"), paper = "a4")
hist(log10(rowSums(cpm(counts(x))) + 1), breaks = 50, main = "Log10 Total CPM per Gene", xlab = "Log10(Total CPM + 1)")
dev.off()

# Distribution of log10 total counts per gene
pdf(file = file.path(qc_dir, "hist_log10_counts_total.pdf"), paper = "a4")
hist(log10(rowSums(counts(x)) + 1), breaks = 50, main = "Log10 Total Raw Counts per Gene", xlab = "Log10(Total Counts + 1)")
dev.off()

# --- 4. DESeq2 Analysis ---

# Filter low count genes (require >= 5 counts in at least 3 samples)
# Using 3 samples as a minimum group size heuristic - adjust if needed
min_count <- 5
min_samples <- 3
keep_genes <- rowSums(counts(x) >= min_count) >= min_samples
cat("Filtering: Keeping", sum(keep_genes), "out of", nrow(x), "genes.\n")
y_filt <- x[keep_genes, ]

# Identify ERCC spike-ins
is_ercc <- rowData(y_filt)[, "Genetype"] == "ERCC-SPIKEIN"
cat("Number of ERCC spike-ins remaining after filtering:", sum(is_ercc), "\n")

# Check if sufficient non-ERCC genes remain for size factor estimation
if(sum(!is_ercc) < 10) {
  warning("Very few non-ERCC genes remain after filtering. Size factor estimation might be unreliable.")
  control_genes <- NULL # Use all genes if too few non-ERCC
} else {
  control_genes <- which(!is_ercc)
}

# Create DESeqDataSet (design includes the 'Group' factor for pairwise comparisons)
# Ensure 'Group' combines factors of interest if not already done
# Example: colData(y_filt)$Group <- factor(paste(colData(y_filt)$CellLine, colData(y_filt)$Treatment, colData(y_filt)$Time, sep="_"))
if(length(design(y_filt)) == 0 || !"Group" %in% all.vars(design(y_filt))) {
  cat("Setting design formula to ~ Group\n")
  design(y_filt) <- ~ Group
} else {
  cat("Using existing design formula:", as.character(design(y_filt)), "\n")
}

dds <- DESeqDataSet(y_filt, design = design(y_filt))

# Estimate size factors using non-ERCC genes if possible
cat("Estimating size factors...\n")
dds <- estimateSizeFactors(dds, controlGenes = control_genes)
print("Size Factors:")
print(sizeFactors(dds))
write.table(x = as.data.frame(sizeFactors(dds)),
            file = file.path(results_dir, "deseq2_sizefactors.txt"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA) # Use col.names=NA for better reading in spreadsheets

# Estimate dispersions
cat("Estimating dispersions...\n")
dds <- estimateDispersions(dds)

# Perform Wald test
cat("Running Wald test...\n")
dds <- nbinomWaldTest(dds)

# Save the processed DESeqDataSet object
deseq2_rdata_file <- file.path(results_dir, "deseq2_object.RData")
save(list = "dds", file = deseq2_rdata_file)
cat("Saved DESeqDataSet object 'dds' to:", deseq2_rdata_file, "\n")


# --- 5. Variance Stabilization and PCA ---

cat("Performing variance stabilizing transformation (VST)...\n")
vsd <- vst(dds, blind = TRUE) # Use blind=TRUE for QC, blind=FALSE if using for downstream DE testing input (not typical)

# PCA Plot - Top N Genes
pca_plot_ntop <- make.pca.ggplot2(vsd, ntop = N_TOP_PCA, pcs = c(1, 2)) + ggtitle(paste("PCA - Top", N_TOP_PCA, "genes"))
ggsave(file.path(qc_dir, paste0("pca_PC1_2_ntop", N_TOP_PCA, ".pdf")), pca_plot_ntop)

pca_plot_ntop_34 <- make.pca.ggplot2(vsd, ntop = N_TOP_PCA, pcs = c(3, 4)) + ggtitle(paste("PCA - Top", N_TOP_PCA, "genes"))
ggsave(file.path(qc_dir, paste0("pca_PC3_4_ntop", N_TOP_PCA, ".pdf")), pca_plot_ntop_34)

# Pairs Plot - Top N Genes (PCs 1-4)
# Need to extract legend from one plot first
pca_data_ntop <- make.pca.ggplot2(vsd, ntop = N_TOP_PCA)$data # Get data used in PCA plot
pca_legend_ntop <- GGally::grab_legend(pca_plot_ntop)

pca_pairs_ntop <- ggpairs(pca_data_ntop,
                          columns = 1:4, # PC1 to PC4
                          switch = "y",
                          ggplot2::aes(shape = Time, fill = CellLine, colour = Treatment),
                          upper = "blank",
                          legend = pca_legend_ntop,
                          axisLabels = "internal") +
  theme_bw() +
  scale_shape_manual(values = 21:24) +
  scale_color_manual(values = c("grey50", "black")) # Match pca_plot_ntop colors
ggsave(file.path(qc_dir, paste0("pca_pairs_1to4_ntop", N_TOP_PCA, ".pdf")), pca_pairs_ntop, width = 10, height = 10)

# PCA Plot - All Genes
pca_plot_all <- make.pca.ggplot2(vsd, ntop = NULL, pcs = c(1, 2)) + ggtitle("PCA - All genes")
ggsave(file.path(qc_dir, "pca_PC1_2_all.pdf"), pca_plot_all)

pca_plot_all_34 <- make.pca.ggplot2(vsd, ntop = NULL, pcs = c(3, 4)) + ggtitle("PCA - All genes")
ggsave(file.path(qc_dir, "pca_PC3_4_all.pdf"), pca_plot_all_34)

# Pairs Plot - All Genes (PCs 1-4)
pca_data_all <- make.pca.ggplot2(vsd, ntop = NULL)$data
pca_legend_all <- GGally::grab_legend(pca_plot_all)

pca_pairs_all <- ggpairs(pca_data_all,
                         columns = 1:4,
                         switch = "y",
                         ggplot2::aes(shape = Time, fill = CellLine, colour = Treatment),
                         upper = "blank",
                         legend = pca_legend_all,
                         axisLabels = "internal") +
  theme_bw() +
  scale_shape_manual(values = 21:24) +
  scale_color_manual(values = c("grey50", "black"))
ggsave(file.path(qc_dir, "pca_pairs_1to4_all.pdf"), pca_pairs_all, width = 10, height = 10)


# --- 6. All Pairwise Comparisons ---

reslist_deseq <- list()
all_groups <- levels(colData(dds)$Group)

# Prepare all unique pairs of groups for comparison
group_pairs <- combn(all_groups, 2, simplify = FALSE)

# Function to perform differential expression analysis for a pair of groups
run_deseq2_pair <- function(pair, dds_obj) {
  group1 <- pair[1]
  group2 <- pair[2]
  comp_name <- paste(group1, "vs", group2, sep = "_")
  cat("Running comparison:", comp_name, "\n")
  
  # Use tryCatch for robustness
  result <- tryCatch({
    # Use lfcShrink with ashr for moderated LFCs
    res_shrink <- lfcShrink(dds = dds_obj,
                            contrast = c("Group", group1, group2),
                            type = "ashr",
                            quiet = TRUE)
    
    # Add metadata (Symbol, Genetype) from rowData
    res_df <- as.data.frame(res_shrink)
    res_df <- cbind(rowData(dds_obj)[rownames(res_df), c("Symbol", "Genetype")], res_df) # Assuming Symbol/Genetype are in rowData
    
    # Save individual results to CSV
    output_filename <- file.path(pairwise_dir, paste0(comp_name, ".csv"))
    write.csv(res_df, file = output_filename, row.names = TRUE) # Save rownames (gene IDs)
    
    return(res_df) # Return the data frame directly
    
  }, error = function(e) {
    cat("ERROR processing comparison:", comp_name, "\n", e$message, "\n")
    return(NULL) # Return NULL on error
  })
  return(result)
}

# Set up parallel backend
register(MulticoreParam(workers = N_CPU))
cat("Running DESeq2 pairwise comparisons using", N_CPU, "cores...\n")

# Run DESeq2 analysis in parallel
reslist_deseq_raw <- bplapply(group_pairs, run_deseq2_pair, dds_obj = dds)

# Filter out NULL results from failed comparisons and name the list
reslist_deseq <- reslist_deseq_raw[!sapply(reslist_deseq_raw, is.null)]
names(reslist_deseq) <- sapply(group_pairs[!sapply(reslist_deseq_raw, is.null)], function(pair) paste(pair, collapse = "_vs_"))

cat("Finished", length(reslist_deseq), "pairwise comparisons.\n")

# --- 7. Generate Plots for Each Comparison ---

cat("Generating plots for each comparison...\n")
for (lab in names(reslist_deseq)) {
  cat("Plotting:", lab, "\n")
  res_df <- reslist_deseq[[lab]] # Already a data frame
  
  if(is.null(res_df) || nrow(res_df) == 0) {
    cat("Skipping plots for", lab, "- results are empty.\n")
    next
  }
  
  # Make histogram of p-values
  # Use padj for histogram to assess distribution after adjustment
  p_hist_data <- res_df[!is.na(res_df$padj), ]
  if(nrow(p_hist_data) > 0) {
    p_hist <- ggplot(p_hist_data, aes(x = padj)) +
      geom_histogram(aes(y = ..density..), bins = 50, fill = "grey", color="black") +
      theme_bw() +
      ggtitle(paste("Adj. P-value Distribution:", lab)) +
      xlab("Adjusted P-value (padj)")
    ggsave(file.path(pairwise_dir, paste0(lab, "_hist_padj.pdf")), p_hist)
  } else {
    cat("No non-NA adjusted p-values for histogram:", lab, "\n")
  }
  
  
  # Make MA plot
  # Ensure function handles potential missing columns gracefully if needed
  tryCatch({
    p_ma <- ma_ggplot2(res_df, lfc = LFC_THRESHOLD, fdr = FDR_THRESHOLD, comp_label = lab)
    ggsave(file.path(pairwise_dir, paste0(lab, "_MA_plot_LFC", LFC_THRESHOLD, "_FDR", FDR_THRESHOLD, ".pdf")), p_ma)
  }, error = function(e) { cat("Error generating MA plot for", lab, ":", e$message, "\n") })
  
  
  # Make volcano plot
  tryCatch({
    p_volcano <- volcano_ggplot2(res_df, fdr_threshold = FDR_THRESHOLD, lfc_threshold = LFC_THRESHOLD, topN = 20, title = lab, use_padj = TRUE)
    ggsave(file.path(pairwise_dir, paste0(lab, "_Volcano_plot_LFC", LFC_THRESHOLD, "_FDR", FDR_THRESHOLD, ".pdf")), p_volcano)
  }, error = function(e) { cat("Error generating Volcano plot for", lab, ":", e$message, "\n") })
  
  
  # Make heatmap of significant genes (using VST data)
  sig_genes_idx <- which(!is.na(res_df$padj) & res_df$padj <= FDR_THRESHOLD & abs(res_df$log2FoldChange) >= LFC_THRESHOLD)
  num_sig <- length(sig_genes_idx)
  cat("Number of significant genes (FDR<", FDR_THRESHOLD, ", |LFC|>", LFC_THRESHOLD, "):", num_sig, "\n")
  
  if (num_sig > 1) { # Need at least 2 genes for heatmap scaling
    sig_gene_ids <- rownames(res_df)[sig_genes_idx]
    heatmap_mat <- assay(vsd)[sig_gene_ids, ]
    
    # Handle cases with too many genes for visualization
    show_rownames_heatmap <- num_sig <= 100 # Show names only if manageable
    kmeans_k_heatmap <- NA
    main_title <- paste(lab, "\nSig. Genes:", num_sig, "(FDR<", FDR_THRESHOLD, ", |LFC|>", LFC_THRESHOLD, ")")
    if (num_sig > 2000) {
      kmeans_k_heatmap <- 8 # Use kmeans clustering for very large heatmaps
      main_title <- paste(main_title, "- K-means (k=8)")
      show_rownames_heatmap <- FALSE
    } else if (num_sig > 100) {
      show_rownames_heatmap <- FALSE
    }
    
    # Get annotation data and colors
    coldata_heatmap <- as.data.frame(colData(vsd)[, c("CellLine", "Time", "Treatment")])
    # Define annotation colors manually or use a helper function
    ann_colors <- list(
      Time = brewer.pal(length(levels(vsd$Time)), "BuPu"),
      CellLine = brewer.pal(length(levels(vsd$CellLine)), "Accent"),
      Treatment = c("grey", "black") # Assuming 2 levels
    )
    names(ann_colors$Time) <- levels(vsd$Time)
    names(ann_colors$CellLine) <- levels(vsd$CellLine)
    names(ann_colors$Treatment) <- levels(vsd$Treatment)
    
    # Generate heatmap
    heatmap_filename <- file.path(pairwise_dir, paste0(lab, "_heatmap_sig_genes_FDR", FDR_THRESHOLD, ".pdf"))
    pheatmap(heatmap_mat,
             scale = "row",
             kmeans_k = kmeans_k_heatmap,
             main = main_title,
             annotation_col = coldata_heatmap,
             annotation_colors = ann_colors,
             show_colnames = FALSE,
             show_rownames = show_rownames_heatmap,
             cluster_rows = TRUE,
             cluster_cols = FALSE, # Often better to keep sample order logical
             fontsize_row = ifelse(show_rownames_heatmap, 6, 4),
             filename = heatmap_filename)
    
  } else {
    cat("Skipping heatmap for", lab, "- fewer than 2 significant genes found.\n")
  }
}


# --- 8. Aggregate Results and Save Tables ---

cat("Aggregating results...\n")

# Extract key stats (LFC, pval, padj) for combined table
reslist_stats <- lapply(reslist_deseq, function(df) {
  df_sub <- df[, c("log2FoldChange", "pvalue", "padj")]
  colnames(df_sub) <- c("lfc", "p", "fdr")
  return(df_sub)
})

# Combine into a single data frame with comparison names as part of column names
restot <- do.call(cbind, reslist_stats)

# Combine gene annotation, baseMean, and aggregated stats
# Use rowData from the *original* unfiltered 'x' to include all gene info if available
# If using 'dds', annotations might be subsetted. Let's use 'dds' for consistency with analysis.
gene_anno <- as.data.frame(rowData(dds)[, c("Symbol", "Description", "Ensembltrans", "Ensemblgene", "Genetype", "basepairs", "Chr")])

# Get baseMean from dds object
base_means <- data.frame(baseMean = rowMeans(counts(dds, normalized=TRUE))) # Or use results(dds)$baseMean if running default results() once

# Combine annotation, baseMean, and stats
# Ensure rownames match for merging
all_genes <- rownames(dds)
mytot <- cbind(gene_anno[all_genes, ], base_means[all_genes, , drop=FALSE], restot[all_genes, ])

# Save aggregated stats table
agg_stats_file <- file.path(tables_dir, "DESeq2_all_comparisons_stats_gencode_hg38.txt")
write.table(x = mytot, file = agg_stats_file, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cat("Saved aggregated stats to:", agg_stats_file, "\n")

# Save normalized & variance stabilized counts
vsd_mat <- assay(vsd)
vsd_df_out <- data.frame(GeneID = rownames(vsd_mat), Symbol = rowData(dds)[rownames(vsd_mat), "Symbol"], vsd_mat)
norm_counts_file <- file.path(tables_dir, "DESeq2_vst_normalized_counts_gencode_hg38.txt")
write.table(x = vsd_df_out, file = norm_counts_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Saved VST normalized counts to:", norm_counts_file, "\n")

# Save key objects for potential downstream use (e.g., script 02)
reslist_rdata_file <- file.path(results_dir, "deseq2_reslist.RData")
save(list = c("reslist_deseq", "restot", "mytot", "vsd"), file = reslist_rdata_file)
cat("Saved aggregated results list, tables, and VST object to:", reslist_rdata_file, "\n")

# --- 9. Plot Counts for Genes of Interest ---

cat("Plotting counts for specific genes:", paste(genes_of_interest, collapse=", "), "\n")
plot_goi_dir <- file.path(pairwise_dir, "Genes_of_Interest_Plots") # Put in subfolder
dir.create(plot_goi_dir, showWarnings = FALSE)

for (g in genes_of_interest) {
  # Check if gene exists in the results (using Symbol from rowData)
  target_gene_id <- rownames(dds)[rowData(dds)$Symbol == g]
  if (length(target_gene_id) == 1) {
    p <- plotCounts2(dds, gene = target_gene_id, groupby = "Time", colorby = "Treatment") +
      facet_wrap(~CellLine, scales = "free_y") +
      theme_bw() +
      ggtitle(paste("Normalized Counts:", g))
    ggsave(file.path(plot_goi_dir, paste0(g, "_gene_counts.pdf")), p, width=8, height=6)
  } else if (length(target_gene_id) > 1) {
    cat("Warning: Multiple gene IDs found for Symbol", g, ". Plotting first one:", target_gene_id[1], "\n")
    # Plot the first one, or consider alternative handling
    p <- plotCounts2(dds, gene = target_gene_id[1], groupby = "Time", colorby = "Treatment") +
      facet_wrap(~CellLine, scales = "free_y") +
      theme_bw() +
      ggtitle(paste("Normalized Counts:", g, "(ID:", target_gene_id[1], ")"))
    ggsave(file.path(plot_goi_dir, paste0(g, "_", target_gene_id[1], "_gene_counts.pdf")), p, width=8, height=6)
  } else {
    cat("Warning: Gene Symbol", g, "not found in results. Skipping count plot.\n")
  }
}

# --- 10. Standard QC Plots (using utility function) ---

cat("Generating standard QC plots...\n")
qc_preproc_dir <- file.path(qc_dir, "preproc_plots")
deseq2.preprocess.plots(dds, output_dir = qc_preproc_dir, lab = "TB_RNAseq")
cat("Standard QC plots saved in:", qc_preproc_dir, "\n")

# --- 11. Overlap Analysis (UpSet Plots) ---

cat("Performing overlap analysis using UpSet plots...\n")

# Extract DE gene lists (using rownames - Gene IDs)
all_de_genes_list <- lapply(reslist_deseq, function(df) {
  rownames(df)[which(!is.na(df$padj) & df$padj <= FDR_THRESHOLD & abs(df$log2FoldChange) >= LFC_THRESHOLD)]
})
down_de_genes_list <- lapply(reslist_deseq, function(df) {
  rownames(df)[which(!is.na(df$padj) & df$padj <= FDR_THRESHOLD & df$log2FoldChange <= -LFC_THRESHOLD)]
})
up_de_genes_list <- lapply(reslist_deseq, function(df) {
  rownames(df)[which(!is.na(df$padj) & df$padj <= FDR_THRESHOLD & df$log2FoldChange >= LFC_THRESHOLD)]
})

# Filter out empty lists to avoid UpSetR errors
all_de_genes_list <- all_de_genes_list[sapply(all_de_genes_list, length) > 0]
down_de_genes_list <- down_de_genes_list[sapply(down_de_genes_list, length) > 0]
up_de_genes_list <- up_de_genes_list[sapply(up_de_genes_list, length) > 0]

# Function to safely generate UpSet plot
safe_upset <- function(data_list, filename, nsets_max = 12) {
  if(length(data_list) == 0) {
    cat("Skipping UpSet plot:", filename, "- No data sets with DE genes.\n")
    return()
  }
  # Convert list to UpSetR input format using helper function
  upset_data <- upset_fromList(data_list)
  if(nrow(upset_data) == 0 || ncol(upset_data) == 0) {
    cat("Skipping UpSet plot:", filename, "- Conversion from list failed or resulted in empty data.\n")
    return()
  }
  nsets_plot <- min(nsets_max, ncol(upset_data))
  
  pdf(filename, paper = "a4", width = 11, height = 8)
  tryCatch({
    print( # Upset needs explicit print inside loops/functions sometimes
      upset(upset_data, nsets = nsets_plot, order.by = "freq", mb.ratio = c(0.6, 0.4), text.scale = 1.2)
    )
  }, error = function(e) {
    cat("Error generating UpSet plot:", filename, "\n", e$message, "\n")
    plot.new(); title(main="Error generating UpSet plot", sub=e$message) # Placeholder plot on error
  })
  dev.off()
}

# Generate UpSet plots
safe_upset(all_de_genes_list, filename = file.path(overlap_dir, paste0("upset_overlap_all_DE_FDR", FDR_THRESHOLD, "_LFC", LFC_THRESHOLD, ".pdf")))
safe_upset(up_de_genes_list, filename = file.path(overlap_dir, paste0("upset_overlap_up_DE_FDR", FDR_THRESHOLD, "_LFC", LFC_THRESHOLD, ".pdf")))
safe_upset(down_de_genes_list, filename = file.path(overlap_dir, paste0("upset_overlap_down_DE_FDR", FDR_THRESHOLD, "_LFC", LFC_THRESHOLD, ".pdf")))


# --- 12. Specific Overlap Example (Down-regulated 3h vs 0h) ---

cat("Analyzing common down-regulated genes in '3h_vs_0h' comparisons...\n")

# Identify relevant comparisons (adjust pattern if needed)
pattern_3h_vs_0h <- "3h_vs_0h$" # Assuming group names end like this
comps_3h_vs_0h_down <- grep(pattern_3h_vs_0h, names(down_de_genes_list), value = TRUE)

if(length(comps_3h_vs_0h_down) > 1) {
  # Get the list of down-regulated genes for these comparisons
  down_list_subset <- down_de_genes_list[comps_3h_vs_0h_down]
  
  # Find common genes using intersect or matrix approach
  upset_subset_data <- upset_fromList(down_list_subset)
  common_genes_down <- rownames(upset_subset_data)[rowSums(upset_subset_data) == length(down_list_subset)]
  num_common <- length(common_genes_down)
  cat("Found", num_common, "genes commonly down-regulated in:", paste(comps_3h_vs_0h_down, collapse=", "), "\n")
  
  # Prepare annotation colors (reuse from section 7)
  coldata_heatmap <- as.data.frame(colData(vsd)[, c("CellLine", "Time", "Treatment")])
  ann_colors <- list(
    Time = brewer.pal(length(levels(vsd$Time)), "BuPu"),
    CellLine = brewer.pal(length(levels(vsd$CellLine)), "Accent"),
    Treatment = c("grey", "black") # Assuming 2 levels
  )
  names(ann_colors$Time) <- levels(vsd$Time)
  names(ann_colors$CellLine) <- levels(vsd$CellLine)
  names(ann_colors$Treatment) <- levels(vsd$Treatment)
  
  # Heatmap of common down-regulated genes across all samples
  if(num_common > 1) {
    pheatmap(assay(vsd)[common_genes_down, ],
             scale = "row",
             kmeans_k = NA,
             annotation_col = coldata_heatmap,
             annotation_colors = ann_colors,
             show_colnames = FALSE,
             show_rownames = num_common <= 100, # Show names only if reasonable
             cluster_cols = FALSE,
             main = paste("Common Down-regulated Genes in 3h vs 0h Comparisons (N =", num_common, ")"),
             filename = file.path(overlap_dir, "heatmap_common_down_3h_vs_0h.pdf"))
  }
  
  # Optional: Heatmaps per cell line for their specific 3h vs 0h down-regulated genes
  # for (cur_cl in levels(coldata_heatmap$CellLine)) {
  #   cat("Heatmap for CellLine:", cur_cl, "\n")
  #   cur_comp_name <- grep(paste0(cur_cl, ".*", pattern_3h_vs_0h), names(down_list_subset), value = TRUE)
  #   if(length(cur_comp_name) == 1) {
  #        cur_down_genes <- down_list_subset[[cur_comp_name]]
  #        num_cur_down <- length(cur_down_genes)
  #        if(num_cur_down > 1) {
  #           samples_to_include <- grep(cur_cl, colnames(vsd)) # Find columns matching cell line
  #           pheatmap(assay(vsd)[cur_down_genes, samples_to_include],
  #                    scale = "row",
  #                    kmeans_k = NA,
  #                    annotation_col = coldata_heatmap[samples_to_include, c("Time", "Treatment"), drop=FALSE],
  #                    annotation_colors = ann_colors[c("Time", "Treatment")],
  #                    show_colnames = FALSE,
  #                    show_rownames = num_cur_down <= 100,
  #                    cluster_cols = FALSE,
  #                    main = paste(cur_cl, ": Down-regulated in 3h vs 0h (N =", num_cur_down, ")"),
  #                    filename = file.path(overlap_dir, paste0("heatmap_", cur_cl, "_down_3h_vs_0h.pdf")))
  #        }
  #   }
  # }
  
  # Venn Diagram for the overlap (if using ggvenn)
  if (length(down_list_subset) >= 2 && length(down_list_subset) <= 6) { # ggvenn works best for 2-6 sets
    library(ggvenn)
    p_venn <- ggvenn(down_list_subset, set_name_size = 3) + ggtitle("Overlap of Down-regulated Genes (3h vs 0h)")
    ggsave(file.path(overlap_dir, "venn_down_3h_vs_0h.pdf"), p_venn)
  }
  
} else {
  cat("Skipping specific '3h vs 0h' common down-regulated analysis - fewer than 2 relevant comparisons found.\n")
}


# --- 13. Final Cleanup / Zipping (Optional) ---

# Optional: Zip results folder (can be system dependent)
# zip_file <- file.path(results_dir, "deseq2_analysis_results.zip")
# zip_command <- paste("zip -r", shQuote(zip_file), shQuote(basename(results_dir)))
# cat("Attempting to zip results folder...\n")
# system(paste("cd", shQuote(dirname(results_dir)), "&&", zip_command)) # Navigate, then zip

cat("--- DESeq2 Analysis Script Finished ---\n")
