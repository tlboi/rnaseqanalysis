# --- START OF FILE scripts/utils_deseq_functions.R ---

######################################################################################
# Utility Functions for DESeq2 RNA-Seq Analysis
# Assumes necessary libraries (ggplot2, DESeq2, pheatmap, etc.) are loaded
# by the calling script.
######################################################################################


#' Generate Standard Preprocessing Plots for DESeq2 Analysis
#'
#' @param y A DESeqDataSet object (post DESeq() run).
#' @param output_dir Directory path to save the plots.
#' @param lab Label for filenames (e.g., "rnaseq_project").
#' @param varStabTransform Method for variance stabilization ("vst" or "rlog").
#'
#' @return Invisible NULL. Plots are saved to files.
#' @details Generates plots for size factors, zero counts, data density,
#'          sample distances, PCA, PC correlations, variance stabilization,
#'          dispersion, sparsity, and Cook's distance.
deseq2.preprocess.plots <- function(y, output_dir = ".", lab = "rnaseq", varStabTransform = "vst") {
  #input must be class "DESeqDataSet" of 'DESeq' processed data
  
  if (class(y) != "DESeqDataSet") {
    stop("Input must be class 'DESeqDataSet' of 'DESeq' processed data")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Required libraries (ensure they are loaded by the main script)
  # library(tidyr)
  # library(DESeq2)
  # library(ggplot2)
  # library(pheatmap)
  # library(vsn)
  # library(pcaExplorer)
  # library(reshape2)
  
  # --- Plotting Functions ---
  
  # Plot sizeFactors
  pdf(file.path(output_dir, paste("sizefactors", lab, "pdf", sep = ".")), paper = "a4")
  par(mar = c(5, 12, 4, 2))
  barplot(sizeFactors(y), horiz = TRUE, las = 2, main = "SizeFactors")
  dev.off()
  
  # Plot fraction of genes with zero count
  pdf(file.path(output_dir, paste("fraction_zero_count_by_sample", lab, "pdf", sep = ".")), paper = "a4")
  fraction_zero_count <- apply(counts(y), 2, function(v) 100 * (sum(v == 0) / length(v)))
  par(mar = c(5, 12, 4, 2))
  barplot(fraction_zero_count, horiz = TRUE, las = 2, xlab = "Percent", main = "Percent zero count genes")
  dev.off()
  
  # Density plot of log2 raw counts+1
  # Ensure tidyr::gather is available if using older tidyr, otherwise use pivot_longer
  count_data_long <- tidyr::pivot_longer(data.frame(log2(counts(y) + 1)), cols = everything(), names_to = "key", values_to = "value")
  p <- ggplot(count_data_long, aes(x = value, colour = key)) + geom_density() + ggtitle("log2+1 transformed raw counts") + theme(legend.position = "none")
  ggsave(filename = file.path(output_dir, paste("density_log2_rawcounts", lab, "pdf", sep = ".")), plot = p)
  
  # Perform variance stabilization
  if (varStabTransform == "rlog") {
    yt <- rlog(y, blind = TRUE) # Use rlog directly
  } else {
    yt <- vst(y, blind = TRUE) # Use vst directly
  }
  
  # Extract colData
  coldata <- as.data.frame(colData(y))
  # Clean up coldata if needed (remove columns that are all identical or problematic)
  # coldata <- coldata[, sapply(coldata, function(x) length(unique(x)) > 1)]
  
  # Heatmap of sample distances (using correlation on transformed data)
  pdf(file.path(output_dir, paste("heatmap_sample2sample_dist", lab, "pdf", sep = ".")), paper = "a4")
  mat <- as.matrix(cor(assay(yt)))
  pheatmap::pheatmap(mat, annotation_col = coldata[, sapply(coldata, is.factor), drop = FALSE]) # Show factor annotations
  dev.off()
  
  # PCA plots for each factor in colData
  mytypes <- colnames(coldata)[sapply(coldata, function(x) is.factor(x) || is.character(x))] # Plot factors/characters
  for (mytype in mytypes) {
    if (length(unique(coldata[[mytype]])) > 1 && length(unique(coldata[[mytype]])) < nrow(coldata)) { # Avoid plotting unique IDs or constant columns
      tryCatch({
        p <- DESeq2::plotPCA(yt, intgroup = mytype, ntop = 500) + ggtitle(paste("PCA by", mytype, "(Top 500 genes)"))
        ggsave(filename = file.path(output_dir, paste("PCA", mytype, lab, "ntop500", "pdf", sep = ".")), plot = p)
      }, error = function(e) {
        cat("Could not generate PCA for:", mytype, "\nError:", e$message, "\n")
      })
    }
  }
  
  # PC and covariate correlations
  if (nrow(assay(yt)) > ncol(assay(yt))) { # Ensure samples < features for prcomp on t(assay(yt))
    pcaobj <- prcomp(t(assay(yt)))
    p <- pcaExplorer::pcascree(pcaobj, type = "pev", title = "PCA Scree Plot (Variance Explained)")
    ggsave(file.path(output_dir, paste("pcascree_plot", lab, "pdf", sep = ".")), plot = p)
    
    # Determine correlation of covariates with PCs
    # Select meaningful covariates (numeric or factors with fewer levels than samples)
    covars2plot_idx <- sapply(coldata, function(v) is.numeric(v) || ( (is.factor(v) || is.character(v)) && length(unique(v)) < length(v) && length(unique(v)) > 1) )
    if(sum(covars2plot_idx) > 0) {
      # Convert character columns to factors for correlation analysis
      coldata_factors <- coldata[, covars2plot_idx, drop=FALSE]
      for(cn in names(coldata_factors)) {
        if(is.character(coldata_factors[[cn]])) {
          coldata_factors[[cn]] <- as.factor(coldata_factors[[cn]])
        }
      }
      
      # Calculate correlations (handle potential errors)
      pcs_to_correlate <- 1:min(6, ncol(pcaobj$x))
      tryCatch({
        res_pc <- pcaExplorer::correlatePCs(pcaobj, coldata_factors, pcs = pcs_to_correlate)
        
        # Extract variance explained
        pca_var_expl <- round(100 * (pcaobj$sdev^2) / sum(pcaobj$sdev^2), 1)
        rownames(res_pc) <- paste(paste0("PC", pcs_to_correlate), paste0(pca_var_expl[pcs_to_correlate], "%"), sep = " : ")
        
        # Make barplots of covariate correlations
        res_pc_dat <- reshape2::melt(res_pc, varnames = c("PC", "Covariate"), value.name = "pval")
        res_pc_dat$log10pval <- -log10(res_pc_dat$pval)
        res_pc_dat$log10pval[is.infinite(res_pc_dat$log10pval)] <- max(res_pc_dat$log10pval[is.finite(res_pc_dat$log10pval)], 300) # Cap infinite values
        
        p <- ggplot(data = res_pc_dat, aes(x = Covariate, y = log10pval)) +
          geom_bar(stat = "identity") +
          facet_wrap(~PC) + # Use facet_wrap for clarity
          theme_bw() +
          coord_flip() +
          ylab("-log10(p-value)") +
          ggtitle("Correlation between PCs and Covariates")
        ggsave(file.path(output_dir, paste("pc_covar_cor_barplots", lab, "pdf", sep = ".")), plot = p, width=10, height=7)
      }, error = function(e) {
        cat("Could not perform PC correlation analysis.\nError:", e$message, "\n")
      })
    } else {
      cat("No suitable covariates found for PC correlation analysis.\n")
    }
  } else {
    cat("Skipping PCA Scree Plot and PC correlation: Number of genes is less than or equal to the number of samples.\n")
  }
  
  
  # Plot mean-sd relationship for different transformations
  pdf(file.path(output_dir, paste("variance_stabilizing_transformations", lab, "pdf", sep = ".")), paper = "a4", width = 12, height = 4)
  par(mfrow = c(1, 3))
  notAllZero <- (rowSums(DESeq2::counts(y)) > 0)
  if(sum(notAllZero) > 0) {
    vsn::meanSdPlot(log2(DESeq2::counts(y, normalized = TRUE)[notAllZero, ] + 1), main = "Log2(normalized counts + 1)")
    # meanSdPlot(assay(rld[notAllZero,])) # If using rlog
    vsn::meanSdPlot(assay(yt[notAllZero, ]), main = paste(varStabTransform, "(blind=TRUE)"))
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="") # Placeholder for third plot if needed
  } else {
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="", main="No non-zero genes found.")
  }
  dev.off()
  
  # Dispersion plot
  pdf(file.path(output_dir, paste("dispersion", lab, "pdf", sep = ".")), paper = "a4")
  plotDispEsts(y)
  dev.off()
  
  # Sparsity plot
  pdf(file.path(output_dir, paste("sparsity", lab, "pdf", sep = ".")), paper = "a4")
  plotSparsity(y)
  dev.off()
  
  # Plot Cook's distance boxplot
  cooks_data <- assays(y)[["cooks"]]
  if(!is.null(cooks_data)) {
    # Ensure tidyr::gather is available if using older tidyr, otherwise use pivot_longer
    cooks_long <- tidyr::pivot_longer(data.frame(cooks_data), cols = everything(), names_to = "key", values_to = "value")
    p <- ggplot(cooks_long, aes(y = value, x = key)) +
      scale_y_log10() +
      geom_boxplot(outlier.shape = NA) + # Avoid plotting outliers twice
      geom_jitter(width = 0.2, alpha = 0.3, size=0.5) + # Show individual points
      ggtitle("Cook's distance per sample") +
      theme_bw() +
      coord_flip() +
      xlab("Sample") + ylab("Cook's Distance (log10 scale)")
    ggsave(file.path(output_dir, paste("cooks_distance", lab, "pdf", sep = ".")), plot=p)
  } else {
    cat("Cook's distance not found in assays(y).\n")
  }
  invisible(NULL)
}

#-------------------------------------------------------------------------------------

#' Create PCA Plot using ggplot2
#'
#' @param x A DESeqTransform object, DESeqDataSet object, or numeric matrix (genes x samples).
#' @param ann Annotation data frame (samples x factors). If x is DESeq object, colData is used.
#' @param pcs Principal components to plot (numeric vector of length 2).
#' @param color_fac Column name in 'ann' to use for color aesthetic.
#' @param shape_fac Column name in 'ann' to use for shape aesthetic.
#' @param alpha_fac Column name in 'ann' to use for alpha aesthetic.
#' @param ntop Number of top variable genes to use for PCA. If NULL, uses all genes.
#'
#' @return A ggplot object.
make.pca.ggplot <- function(x, ann = NULL, pcs = c(1, 2), color_fac = NULL, shape_fac = NULL, alpha_fac = NULL, size = 3, ntop = NULL) {
  # library(ggplot2)
  # library(DESeq2)
  # library(SummarizedExperiment) # For colData, assay
  
  # Extract data matrix & annotation
  if (inherits(x, "DESeqTransform")) {
    mat <- assay(x)
    if (is.null(ann)) ann <- as.data.frame(colData(x))
  } else if (inherits(x, "DESeqDataSet")) {
    # Use variance stabilized data if available and ntop is specified, otherwise normalized counts
    if (!is.null(ntop)) {
      cat("Using variance stabilized data (VST) for PCA with ntop.\n")
      mat <- assay(vst(x, blind = TRUE))
    } else {
      cat("Using normalized counts for PCA (ntop=NULL).\n")
      mat <- counts(x, normalized = TRUE)
    }
    if (is.null(ann)) ann <- as.data.frame(colData(x))
  } else if (is.matrix(x) && is.numeric(x)) {
    mat <- x
    if (is.null(ann)) stop("Annotation 'ann' must be provided if 'x' is a matrix.")
    if (nrow(ann) != ncol(mat)) stop("Number of rows in 'ann' must match number of columns in 'mat'.")
  } else {
    stop("'x' must be a DESeqTransform, DESeqDataSet, or numeric matrix.")
  }
  
  # Use only top N genes with highest variance
  if (!is.null(ntop)) {
    rv <- matrixStats::rowVars(mat)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    mat <- mat[select, ]
  }
  
  # Perform PCA
  pca <- prcomp(t(mat))
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  # Create data frame for ggplot
  pcs.lab <- paste0("PC", pcs)
  if(ncol(pca$x) < max(pcs)) {
    stop(paste("Cannot plot PC", max(pcs), "- only", ncol(pca$x), "PCs available."))
  }
  dat <- data.frame(pca$x[, pcs, drop = FALSE], ann)
  colnames(dat)[1:2] <- pcs.lab # Ensure PC column names are correct
  
  # Build ggplot object
  aes_map <- aes_string(x = pcs.lab[1], y = pcs.lab[2], color = color_fac, shape = shape_fac, alpha = alpha_fac)
  g <- ggplot(dat, aes_map) +
    geom_point(size = size) +
    theme_bw() +
    xlab(paste0(pcs.lab[1], ": ", percentVar[pcs[1]], "% variance")) +
    ylab(paste0(pcs.lab[2], ": ", percentVar[pcs[2]], "% variance"))
  
  return(g)
}

#-------------------------------------------------------------------------------------

#' Create PCA Plot using ggplot2 (Tailored for specific 3-factor design)
#'
#' @param x A DESeqTransform object, typically from vst() or rlog().
#' @param pcs Principal components to plot (numeric vector of length 2).
#' @param ntop Number of top variable genes to use for PCA. If NULL, uses all genes.
#' @param size Point size.
#'
#' @return A ggplot object. Assumes colData has 'Time', 'CellLine', 'Treatment'.
make.pca.ggplot2 <- function(x, pcs = c(1, 2), size = 3, ntop = NULL) {
  # library(ggplot2)
  # library(DESeq2)
  # library(SummarizedExperiment) # For colData, assay
  
  # Extract data matrix & annotation
  if (inherits(x, "DESeqTransform")) {
    mat <- assay(x)
    ann <- as.data.frame(colData(x))
  } else {
    stop("'x' must be a DESeqTransform object (e.g., from vst or rlog).")
  }
  
  # Check for required columns
  required_cols <- c("Time", "CellLine", "Treatment")
  if (!all(required_cols %in% colnames(ann))) {
    stop("colData must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Use only top N genes with highest variance
  if (!is.null(ntop)) {
    rv <- matrixStats::rowVars(mat)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    mat <- mat[select, ]
  }
  
  # Perform PCA
  pca <- prcomp(t(mat))
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  # Create data frame for ggplot
  pcs.lab <- paste0("PC", pcs)
  if(ncol(pca$x) < max(pcs)) {
    stop(paste("Cannot plot PC", max(pcs), "- only", ncol(pca$x), "PCs available."))
  }
  dat <- data.frame(pca$x[, pcs, drop = FALSE], ann)
  colnames(dat)[1:2] <- pcs.lab # Ensure PC column names are correct
  
  # Plot requested principal components
  p <- ggplot(dat, aes_string(x = pcs.lab[1], y = pcs.lab[2], shape = "Time")) +
    geom_point(aes(fill = CellLine, color = Treatment), size = size) +
    scale_color_manual(values = c("grey50", "black")) + # Assuming 2 levels for Treatment
    scale_shape_manual(values = 21:24) + # Assuming up to 4 levels for Time
    guides(fill = guide_legend(override.aes = list(shape = 21))) + # Ensure fill legend uses a filled shape
    theme_bw() +
    xlab(paste0(pcs.lab[1], ": ", percentVar[pcs[1]], "% variance")) +
    ylab(paste0(pcs.lab[2], ": ", percentVar[pcs[2]], "% variance"))
  
  return(p)
}

#-------------------------------------------------------------------------------------

#' Create MA Plot using ggplot2 and ggpointdensity
#'
#' @param res A DESeqResults object or a data frame with 'baseMean', 'log2FoldChange', 'padj', 'Genetype' columns.
#' @param lfc Log2 fold change threshold for highlighting.
#' @param fdr Adjusted p-value (FDR) threshold for highlighting.
#' @param comp_label Optional label for the plot title/subtitle.
#'
#' @return A ggplot object.
ma_ggplot2 <- function(res, lfc = 1, fdr = 0.05, comp_label = "") {
  # require(ggplot2)
  # require(ggpointdensity) # Ensure this is installed and loaded
  
  res_df <- as.data.frame(res)
  
  # Ensure required columns exist
  req_cols <- c("baseMean", "log2FoldChange", "padj")
  if (!all(req_cols %in% colnames(res_df))) {
    stop("Input data frame must contain columns: ", paste(req_cols, collapse=", "))
  }
  # Add Genetype if missing, needed for ERCC highlighting
  if (!"Genetype" %in% colnames(res_df)) {
    res_df$Genetype <- "Unknown"
    warning("Column 'Genetype' not found. ERCC spike-ins will not be highlighted.")
  }
  
  
  # Add significance column
  res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < fdr & abs(res_df$log2FoldChange) > lfc, "Significant", "Not Significant")
  res_df$highlight_group <- ifelse(res_df$Genetype == "ERCC-SPIKEIN", "ERCC", res_df$significant)
  
  # Determine limits and points for highlighting
  plot_data <- subset(res_df, baseMean > 0) # Avoid log10(0)
  ercc_data <- subset(plot_data, highlight_group == "ERCC")
  sig_data <- subset(plot_data, highlight_group == "Significant")
  
  p <- ggplot(plot_data, aes(x = log10(baseMean), y = log2FoldChange)) +
    geom_pointdensity(size = 0.8, alpha = 0.7) + # Use point density for the background
    scale_color_viridis_c(name = "Density") + # Default viridis color scale for density
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    geom_hline(yintercept = c(-lfc, lfc), colour = "lightblue", linetype = "dotted") +
    geom_point(data = ercc_data, colour = "grey60", size = 1, alpha = 0.8) + # Overlay ERCC points
    geom_point(data = sig_data, colour = "red", size = 1, alpha = 0.8) +   # Overlay Significant points
    theme_bw() +
    labs(
      x = "Log10 Mean Normalized Counts (baseMean)",
      y = "Log2 Fold Change",
      title = "MA Plot",
      subtitle = comp_label
    ) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  # Add annotations for counts (optional)
  # n_sig <- nrow(sig_data)
  # n_ercc <- nrow(ercc_data)
  # p <- p + annotate("text", x = min(log10(plot_data$baseMean)), y = max(plot_data$log2FoldChange),
  #                   label = paste("Significant:", n_sig, "\nERCC:", n_ercc), hjust = 0, vjust = 1, size = 3)
  
  return(p)
}

#-------------------------------------------------------------------------------------

#' Plot Normalized Counts for a Gene using ggplot2
#'
#' @param y A DESeqDataSet object (must have run estimateSizeFactors).
#' @param gene Gene identifier (must be in rownames of y).
#' @param groupby Column name in colData(y) for x-axis grouping.
#' @param colorby Column name in colData(y) for point color.
#'
#' @return A ggplot object.
plotCounts2 <- function(y, gene, groupby, colorby) {
  # library(ggplot2)
  # library(DESeq2)
  # library(SummarizedExperiment) # For colData
  
  if (!gene %in% rownames(y)) {
    stop("Gene '", gene, "' not found in rownames of the DESeqDataSet object.")
  }
  if (!groupby %in% colnames(colData(y))) {
    stop("Group factor '", groupby, "' not found in colData.")
  }
  if (!colorby %in% colnames(colData(y))) {
    stop("Color factor '", colorby, "' not found in colData.")
  }
  
  # Use plotCounts internal mechanism to get data points
  plot_data <- DESeq2::plotCounts(y, gene = gene, intgroup = c(groupby, colorby), returnData = TRUE)
  
  # Create the plot
  p <- ggplot(plot_data, aes_string(x = groupby, y = "count", colour = colorby)) +
    geom_point(position = position_jitter(width = 0.1, height = 0), size = 2) +
    scale_y_log10(name = paste("Normalized Counts (log10 scale) -", gene)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) # Improved angle
  
  return(p)
}

#-------------------------------------------------------------------------------------

#' Create Volcano Plot using ggplot2 and ggrepel
#'
#' @param myresults A DESeqResults object or a data frame. Must contain 'log2FoldChange', 'pvalue', 'padj'.
#'                  Also uses 'Symbol' (or 'gene') for labeling.
#' @param fdr_threshold Adjusted p-value threshold for significance.
#' @param lfc_threshold Absolute log2 fold change threshold for significance.
#' @param topN Number of top significant genes to label (by smallest padj).
#' @param label_genes A character vector of specific gene symbols to label, regardless of rank.
#' @param title Optional plot title.
#' @param use_padj Use adjusted p-value ('padj') for y-axis? Default TRUE. If FALSE, uses 'pvalue'.
#'
#' @return A ggplot object.
volcano_ggplot2 <- function(myresults, fdr_threshold = 0.05, lfc_threshold = 1, topN = 10, label_genes = NULL, title = "Volcano Plot", use_padj = TRUE) {
  # library(ggplot2)
  # library(ggrepel)
  
  if (inherits(myresults, "DESeqResults")) {
    results_df <- as.data.frame(myresults)
    results_df$gene <- rownames(myresults) # Ensure gene name is a column
  } else if (is.data.frame(myresults)) {
    results_df <- myresults
    # Try to find a gene name column if 'gene' isn't present
    if (!"gene" %in% colnames(results_df)) {
      if("Symbol" %in% colnames(results_df)) {
        results_df$gene <- results_df$Symbol
        cat("Using 'Symbol' column for gene labels.\n")
      } else if (rownames(results_df)[1] != "1") { # Check if rownames look like gene names
        results_df$gene <- rownames(results_df)
        cat("Using rownames for gene labels.\n")
      } else {
        stop("Could not find a 'gene' or 'Symbol' column, or informative rownames for labels.")
      }
    }
  } else {
    stop("'myresults' must be a DESeqResults object or a data frame.")
  }
  
  # Check for required columns
  value_col <- ifelse(use_padj, "padj", "pvalue")
  req_cols <- c("log2FoldChange", value_col, "gene")
  if (!all(req_cols %in% colnames(results_df))) {
    stop("Input data frame must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # Prepare data: remove NA p-values, calculate -log10 p-value
  results_df <- results_df[!is.na(results_df[[value_col]]), ]
  results_df$log10_p <- -log10(results_df[[value_col]])
  # Handle infinite -log10(p) values by setting them slightly above the max finite value
  max_finite_log10p <- max(results_df$log10_p[is.finite(results_df$log10_p)], na.rm = TRUE)
  results_df$log10_p[is.infinite(results_df$log10_p)] <- max_finite_log10p * 1.05
  
  # Add significance status column (using padj for status, even if plotting pvalue)
  results_df$status <- "Not Significant"
  if("padj" %in% colnames(results_df)) {
    sig_idx <- !is.na(results_df$padj) & results_df$padj <= fdr_threshold & abs(results_df$log2FoldChange) >= lfc_threshold
    results_df$status[sig_idx & results_df$log2FoldChange > lfc_threshold] <- "Up-regulated"
    results_df$status[sig_idx & results_df$log2FoldChange < -lfc_threshold] <- "Down-regulated"
  } else {
    cat("Warning: 'padj' column not found. Significance status based only on p-value and LFC.\n")
    sig_idx <- !is.na(results_df$pvalue) & results_df$pvalue <= fdr_threshold & abs(results_df$log2FoldChange) >= lfc_threshold # Use fdr_threshold conceptually here
    results_df$status[sig_idx & results_df$log2FoldChange > lfc_threshold] <- "Up-regulated"
    results_df$status[sig_idx & results_df$log2FoldChange < -lfc_threshold] <- "Down-regulated"
  }
  results_df$status <- factor(results_df$status, levels = c("Down-regulated", "Not Significant", "Up-regulated"))
  
  
  # Determine genes to label
  significant_genes <- subset(results_df, status != "Not Significant")
  significant_genes <- significant_genes[order(significant_genes[[value_col]]), ] # Order by p-value/padj
  genes_to_label_topN <- head(significant_genes$gene, topN)
  genes_to_label_user <- intersect(label_genes, results_df$gene)
  all_genes_to_label <- unique(c(genes_to_label_topN, genes_to_label_user))
  label_data <- subset(results_df, gene %in% all_genes_to_label)
  
  
  # Define colors
  color_map <- c("Down-regulated" = "blue", "Not Significant" = "grey", "Up-regulated" = "red")
  
  
  # Build ggplot
  p <- ggplot(results_df, aes(x = log2FoldChange, y = log10_p)) +
    geom_point(aes(color = status), alpha = 0.6, size = 0.8) +
    scale_color_manual(name = "Significance", values = color_map) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), lty = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_threshold), lty = "dashed", color = "grey50") + # Line for FDR threshold
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = gene),
      size = 2.5,
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.3, "lines"),
      segment.color = 'grey50',
      max.overlaps = 20 # Increase max overlaps
    ) +
    theme_bw() +
    labs(
      x = "Log2 Fold Change",
      y = paste0("-Log10 ", ifelse(use_padj, "Adjusted P-value", "P-value")),
      title = title
    ) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  
  return(p)
}


#-------------------------------------------------------------------------------------

#' Convert Factor to Color Vector and Palette
#'
#' @param vec A factor or a vector coercible to a factor.
#'
#' @return A list containing 'col_vec' (color vector matching input) and 'col_pal' (named color palette).
fac2col <- function(vec) {
  if (!is.factor(vec)) vec <- as.factor(vec)
  levs <- levels(vec)
  cols <- scales::hue_pal()(length(levs)) # Use ggplot's default color palette generator
  # Alternative: cols <- rainbow(length(levs))
  names(cols) <- levs
  mycols <- cols[match(vec, levs)] # Map colors back to the original vector order
  return(list(col_vec = mycols, col_pal = cols))
}

#-------------------------------------------------------------------------------------

#' Create Binary Matrix from List for UpSet Plot
#'
#' Same as fromList() from 'UpSetR', but keeps element names as rownames.
#'
#' @param input A named list where each element is a vector of items.
#'
#' @return A data frame with rows=unique items, columns=list names, value=1 if item in list element, 0 otherwise.
upset_fromList <- function (input) {
  # Check input is a list
  if (!is.list(input)) stop("Input must be a list.")
  # Check list has names
  if (is.null(names(input))) stop("Input list must be named.")
  
  elements <- unique(unlist(input))
  if (length(elements) == 0) {
    cat("Warning: No elements found in the input list.\n")
    return(data.frame())
  }
  
  # Create matrix 
  mat <- sapply(input, function(set) {
    as.integer(elements %in% set)
  })
  rownames(mat) <- elements
  colnames(mat) <- names(input) 
  
  return(data.frame(mat))
}
