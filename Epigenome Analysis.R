# =============================================================================
# Epigenome-Transcriptome Integration Analysis
# =============================================================================
# 
# Project: Integrated analysis of DNA methylation (RRBS) and gene expression data
# Purpose: Correlate differential methylation patterns with gene expression changes
#          across multiple time points and genomic regions
# 
# Analysis includes:
# - Differential methylated cytosines (DMCs) analysis
# - Differentially expressed proteins (DEPs) analysis  
# - Gene Ontology and pathway enrichment analysis
# - Correlation analysis between methylation and expression
# - Comprehensive visualization suite
#
# Author: Lauren Urban
# =============================================================================

# SETUP AND DEPENDENCIES =====================================================

#' Package installer function
#' @param pkg Vector of package names to install and load
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Define required packages
packages <- c(
  "ggsci", "devtools", "hrbrthemes", "ggplot2", "ggpubr", "dplyr", "ggrepel",
  "tidyverse", "readxl", "readr", "tibble", "reshape2", "gprofiler2",
  "RColorBrewer", "scales", "sm", "viridis", "writexl", "GGally",
  "pheatmap", "ComplexHeatmap", "circlize", "colorspace", "heatmaply",
  "superheat", "lattice", "gplots", "paletteer", "kableExtra"
)

# Install and load packages
suppressMessages(ipak(packages))

# Install additional packages from GitHub if needed
if (!require("ggpubr")) {
  devtools::install_github("kassambara/ggpubr")
  library(ggpubr)
}

# Set global options
options(stringsAsFactors = FALSE, knitr.table.format = "html")
theme_set(theme_bw())

# UTILITY FUNCTIONS ===========================================================

#' Create custom color palette for heatmaps
#' @param n Number of colors to generate
#' @param type Type of palette ("diverging", "sequential")
create_color_palette <- function(n = 15, type = "diverging") {
  if (type == "diverging") {
    return(colorRampPalette(c('#B22526', 'gray98', '#145ec9'), 
                           bias = 1.23, space = "rgb")(n))
  } else {
    return(sequential_hcl(n, "Heat"))
  }
}

#' Load and prepare methylation data
#' @param file_path Path to methylation data file
#' @param sample_cols Vector of sample column names
load_methylation_data <- function(file_path, sample_cols = c("nsc_avg", "d1_avg", "wk1_avg")) {
  data <- read_delim(file_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  data <- data.frame(data)
  data <- na.omit(data)
  
  # Set gene names as row names
  rownames(data) <- data$GeneList
  
  # Select only sample columns for matrix conversion
  mat_data <- data[, sample_cols, drop = FALSE]
  
  return(list(data = data, matrix = as.matrix(mat_data)))
}

#' Filter DMCs by methylation difference and significance
#' @param data Input DMC data
#' @param meth_diff_col Column name for methylation difference
#' @param qval_col Column name for q-values
#' @param meth_threshold Methylation difference threshold (default: 20)
#' @param qval_threshold Q-value threshold (default: 0.01)
filter_significant_dmcs <- function(data, meth_diff_col, qval_col, 
                                   meth_threshold = 20, qval_threshold = 0.01) {
  hyper <- data %>%
    filter(get(meth_diff_col) > meth_threshold, get(qval_col) < qval_threshold) %>%
    mutate(type = "hyper")
  
  hypo <- data %>%
    filter(get(meth_diff_col) < -meth_threshold, get(qval_col) < qval_threshold) %>%
    mutate(type = "hypo")
  
  return(rbind(hyper, hypo))
}

# DATA PROCESSING FUNCTIONS ===================================================

#' Process DEP (Differentially Expressed Proteins) data
#' @param dep_file Path to DEP data file
process_dep_data <- function(dep_file) {
  deps_df <- read_delim(dep_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  deps_df <- data.frame(deps_df)
  deps_df <- na.omit(deps_df)
  
  return(deps_df)
}

#' Process unique DMCs for specific time points
#' @param promoter_file Path to promoter DMC file
#' @param gene_file Path to gene body DMC file
#' @param time_point Time point identifier ("Day1" or "Wk1")
process_unique_dmcs <- function(promoter_file, gene_file, time_point) {
  # Load data
  promoter_data <- read_delim(promoter_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  gene_data <- read_delim(gene_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Add location annotations
  promoter_data$loc <- "Promoters"
  gene_data$loc <- "Gene_bodies"
  
  # Combine datasets
  combined_data <- rbind(data.frame(promoter_data), data.frame(gene_data))
  
  return(combined_data)
}

#' Create correlation scatter plot for methylation vs expression
#' @param data Combined methylation and expression data
#' @param x_col X-axis variable (methylation)
#' @param y_col Y-axis variable (expression)
#' @param color_by Grouping variable for colors
#' @param facet_by Variable for faceting
create_correlation_plot <- function(data, x_col, y_col, color_by = NULL, facet_by = NULL) {
  # Create significance categories
  data <- data %>%
    mutate(significance = case_when(
      get(y_col) < -1 ~ 'downregulated',
      get(y_col) > 1 ~ 'upregulated',
      TRUE ~ 'not significant'
    ))
  
  # Get top genes for labeling
  if ("loc" %in% names(data)) {
    top_genes <- data %>%
      filter(loc == "Promoters", significance != 'not significant') %>%
      slice_max(abs(get(x_col)), n = 20) %>%
      pull(GeneList)
    
    data$labels <- ifelse(data$GeneList %in% top_genes, data$GeneList, "")
  }
  
  # Create base plot
  p <- ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = significance), alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "darkgray") +
    scale_color_manual(values = c("downregulated" = "cornflowerblue", 
                                 "not significant" = "grey", 
                                 "upregulated" = "firebrick")) +
    geom_vline(xintercept = c(-20, 20), linetype = "dotted", color = "black") +
    geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "black") +
    labs(x = "Methylation Difference (%)", 
         y = "Expression Log2 Fold Change",
         color = "Expression Status") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  
  # Add gene labels if available
  if ("labels" %in% names(data)) {
    p <- p + geom_text_repel(aes(label = labels), size = 3, max.overlaps = 20)
  }
  
  # Add faceting if specified
  if (!is.null(facet_by) && facet_by %in% names(data)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free")
  }
  
  # Add correlation statistics
  p <- p + stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.9)
  
  return(p)
}

# GO ENRICHMENT ANALYSIS FUNCTIONS ===========================================

#' Perform GO enrichment analysis using gprofiler2
#' @param gene_list Vector of gene symbols or IDs
#' @param organism Organism name (default: "hsapiens")
#' @param sources Data sources to query (default: all)
#' @param correction_method Multiple testing correction method
run_go_enrichment <- function(gene_list, organism = "hsapiens", 
                             sources = c("GO", "KEGG", "REAC"), 
                             correction_method = "g_SCS") {
  
  # Remove any NA or empty values
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  
  # Run GO enrichment
  gostres <- gost(
    query = gene_list,
    organism = organism,
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = correction_method,
    domain_scope = "annotated",
    sources = sources
  )
  
  return(gostres)
}

#' Create publication-ready GO plot
#' @param gostres Results from gost() function
#' @param highlight_terms Vector of terms to highlight
#' @param title Plot title
create_go_plot <- function(gostres, highlight_terms = NULL, title = "GO Enrichment") {
  
  # Create the base plot
  p <- gostplot(gostres, capped = TRUE, interactive = FALSE)
  
  # Add highlighting if specified
  if (!is.null(highlight_terms)) {
    p <- publish_gostplot(p, highlight_terms = highlight_terms)
  }
  
  return(p)
}

# ADVANCED VISUALIZATION FUNCTIONS ===========================================

#' Create violin plot for methylation data
#' @param data Input data frame
#' @param x_col X-axis variable (typically gene names)
#' @param y_col Y-axis variable (methylation values)
#' @param color_by Coloring variable
#' @param orientation Plot orientation ("vertical" or "horizontal")
create_violin_plot <- function(data, x_col, y_col, color_by = NULL, 
                              orientation = "vertical", title = "Violin Plot") {
  
  # Set default color if not specified
  if (is.null(color_by)) color_by <- x_col
  
  # Create base plot
  p <- ggviolin(data, x = x_col, y = y_col,
                add = "jitter",
                orientation = ifelse(orientation == "horizontal", "horiz", "vert"),
                color = color_by,
                font.label = list(size = 8, color = "black")) +
    labs(title = title, x = x_col, y = y_col) +
    theme_bw() +
    theme(legend.position = ifelse(length(unique(data[[color_by]])) > 10, "none", "right"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Create parallel coordinate plot for temporal data
#' @param data Input data frame
#' @param value_cols Vector of column names containing values
#' @param group_col Column name for grouping
#' @param title Plot title
create_parallel_plot <- function(data, value_cols, group_col, title = "Parallel Coordinate Plot") {
  
  # Prepare data
  plot_data <- data[, c(value_cols, group_col)]
  
  # Create plot
  p <- ggparcoord(plot_data,
                  columns = 1:length(value_cols),
                  groupColumn = group_col,
                  showPoints = TRUE,
                  scale = "globalminmax",
                  title = title,
                  alphaLines = 0.3) +
    scale_color_viridis(discrete = TRUE) +
    theme_ipsum() +
    theme(plot.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Create circular bar plot
#' @param data Input data frame
#' @param x_col X-axis variable (categories)
#' @param y_col Y-axis variable (values) 
#' @param fill_color Fill color for bars
#' @param title Plot title
create_circular_barplot <- function(data, x_col, y_col, fill_color = "#69b3a2", 
                                   title = "Circular Bar Plot") {
  
  # Arrange data by values
  plot_data <- data %>%
    arrange(desc(!!sym(y_col))) %>%
    mutate(!!sym(x_col) := factor(!!sym(x_col), levels = !!sym(x_col)))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = as.factor(!!sym(x_col)), y = !!sym(y_col))) +
    geom_bar(stat = "identity", fill = alpha(fill_color, 0.8)) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm")
    ) +
    coord_polar(start = 0) +
    labs(title = title)
  
  return(p)
}

#' Create lollipop plot
#' @param data Input data frame
#' @param x_col X-axis variable
#' @param y_col Y-axis variable
#' @param color_col Color variable
#' @param baseline Baseline value for segments
create_lollipop_plot <- function(data, x_col, y_col, color_col = NULL, 
                                baseline = 0, title = "Lollipop Plot") {
  
  if (is.null(color_col)) color_col <- x_col
  
  p <- ggplot(data, aes(x = !!sym(x_col), y = !!sym(y_col), color = !!sym(color_col))) +
    geom_segment(aes(x = !!sym(x_col), xend = !!sym(x_col), 
                     y = baseline, yend = !!sym(y_col))) +
    geom_point(size = 4, alpha = 0.7) +
    coord_flip() +
    theme_light() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    ) +
    labs(title = title, x = x_col, y = y_col)
  
  return(p)
}

# PCA ANALYSIS FUNCTIONS ==================================================

#' Perform PCA analysis on expression matrix
#' @param data_matrix Matrix with samples as columns, features as rows
#' @param scale Whether to scale the data
#' @param center Whether to center the data
perform_pca_analysis <- function(data_matrix, scale = TRUE, center = TRUE) {
  
  # Perform PCA (transpose so samples are rows)
  pca <- prcomp(t(data_matrix), scale = scale, center = center)
  
  # Calculate variance explained
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
  
  # Prepare results
  results <- list(
    pca = pca,
    variance_explained = pca.var.per,
    total_variance = pca.var
  )
  
  return(results)
}

#' Create PCA plot
#' @param pca_results Results from perform_pca_analysis()
#' @param pc_x PC for x-axis (default: 1)
#' @param pc_y PC for y-axis (default: 2)
#' @param sample_groups Optional grouping variable for samples
#' @param title Plot title
create_pca_plot <- function(pca_results, pc_x = 1, pc_y = 2, 
                           sample_groups = NULL, title = "PCA Analysis") {
  
  # Extract PCA data
  pca_data <- data.frame(
    Sample = rownames(pca_results$pca$x),
    X = pca_results$pca$x[, pc_x],
    Y = pca_results$pca$x[, pc_y]
  )
  
  # Add grouping if provided
  if (!is.null(sample_groups)) {
    pca_data$Group <- sample_groups
  }
  
  # Create plot
  p <- ggplot(data = pca_data, aes(x = X, y = Y)) +
    {if (!is.null(sample_groups)) geom_point(aes(color = Group), size = 3) 
     else geom_point(size = 3)} +
    geom_text_repel(aes(label = Sample), size = 3) +
    xlab(paste0("PC", pc_x, " - ", pca_results$variance_explained[pc_x], "%")) +
    ylab(paste0("PC", pc_y, " - ", pca_results$variance_explained[pc_y], "%")) +
    theme_bw() +
    labs(title = title)
  
  return(p)
}

#' Create PCA scree plot
#' @param pca_results Results from perform_pca_analysis()
#' @param n_components Number of components to show
create_scree_plot <- function(pca_results, n_components = 10) {
  
  # Prepare data
  scree_data <- data.frame(
    PC = 1:min(n_components, length(pca_results$variance_explained)),
    Variance = pca_results$variance_explained[1:min(n_components, length(pca_results$variance_explained))]
  )
  
  # Create plot
  p <- ggplot(scree_data, aes(x = factor(PC), y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(group = 1), color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(title = "PCA Scree Plot",
         x = "Principal Component",
         y = "Percent Variation") +
    theme_bw()
  
  return(p)
}

#' Get top contributing features for a PC
#' @param pca_results Results from perform_pca_analysis()
#' @param pc_number PC number to analyze
#' @param n_features Number of top features to return
get_top_pc_contributors <- function(pca_results, pc_number = 1, n_features = 10) {
  
  # Get loading scores
  loading_scores <- pca_results$pca$rotation[, pc_number]
  
  # Get magnitudes and rank
  feature_scores <- abs(loading_scores)
  feature_score_ranked <- sort(feature_scores, decreasing = TRUE)
  top_features <- names(feature_score_ranked[1:n_features])
  
  # Return results with scores and directions
  results <- data.frame(
    Feature = top_features,
    Loading_Score = loading_scores[top_features],
    Absolute_Score = feature_scores[top_features]
  )
  
  return(results)
}

# HOMER MOTIF ANALYSIS FUNCTIONS ==========================================

#' Process HOMER motif analysis results
#' @param homer_data Data frame with HOMER results
#' @param value_cols Vector of column names with significance values
#' @param term_col Column name containing motif terms
process_homer_data <- function(homer_data, value_cols, term_col = "terms") {
  
  # Convert to data frame and handle NAs
  homer_df <- data.frame(homer_data, na.strings = "NA")
  
  # Set row names and select value columns
  rownames(homer_df) <- homer_df[[term_col]]
  mat_data <- homer_df[, value_cols, drop = FALSE]
  
  # Convert to matrix
  mat_data <- as.matrix(mat_data)
  
  return(mat_data)
}

#' Create HOMER heatmap using superheat
#' @param homer_matrix Processed HOMER matrix
#' @param heat_limits Vector of min/max values for color scale
#' @param title Plot title
#' @param color_scheme Color scheme ("grey", "viridis", etc.)
create_homer_heatmap <- function(homer_matrix, heat_limits = NULL, 
                                title = "HOMER Motif Analysis", 
                                color_scheme = "grey") {
  
  # Set default limits if not provided
  if (is.null(heat_limits)) {
    heat_limits <- c(min(homer_matrix, na.rm = TRUE), 
                    max(homer_matrix, na.rm = TRUE))
  }
  
  # Create heatmap
  superheat(homer_matrix,
            left.label.text.size = 3,
            bottom.label.text.size = 4,
            heat.lim = heat_limits,
            heat.col.scheme = color_scheme,
            scale = FALSE,
            bottom.label.text.angle = 90,
            grid.hline.col = "white",
            grid.vline.col = "white",
            heat.na.col = "black",
            title = title)
}

# VISUALIZATION FUNCTIONS =====================================================

#' Create enhanced heatmap for expression data
#' @param matrix_data Matrix of expression values
#' @param title Plot title
#' @param cluster_rows Whether to cluster rows
#' @param cluster_cols Whether to cluster columns
#' @param k_means Number of k-means clusters
create_expression_heatmap <- function(matrix_data, title = "Expression Heatmap", 
                                    cluster_rows = TRUE, cluster_cols = TRUE, k_means = 3) {
  
  # Create color palette
  colors <- create_color_palette(15, "diverging")
  
  # Create heatmap
  pheatmap(matrix_data, 
           color = colors,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           k = k_means,
           cutree_rows = k_means,
           border_color = 'gray88',
           main = title,
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 10)
}

#' Create GO enrichment bar plot
#' @param go_data GO enrichment results
#' @param top_n Number of top terms to display
#' @param title Plot title
create_go_barplot <- function(go_data, top_n = 20, title = "GO Enrichment Analysis") {
  
  # Prepare data
  plot_data <- go_data %>%
    arrange(desc(Fold.Enrichment)) %>%
    slice_head(n = top_n) %>%
    mutate(Pathway = factor(Pathway, levels = Pathway))
  
  # Create color scale
  colors <- sequential_hcl(5, "Heat")
  
  # Create plot
  ggplot(plot_data, aes(x = reorder(Pathway, Fold.Enrichment), y = Fold.Enrichment)) +
    geom_col(aes(fill = Enrichment.FDR)) +
    scale_fill_gradientn(colors = colors, limits = c(0, 0.05), 
                        name = "FDR", guide = guide_colorbar(reverse = TRUE)) +
    coord_flip() +
    labs(title = title,
         x = "GO Terms",
         y = "Fold Enrichment") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 8))
}

#' Create volcano plot for differential analysis
#' @param data Differential analysis results
#' @param fc_col Fold change column name
#' @param pval_col P-value column name
#' @param fc_threshold Fold change threshold
#' @param pval_threshold P-value threshold
create_volcano_plot <- function(data, fc_col, pval_col, fc_threshold = 1, pval_threshold = 0.05) {
  
  # Prepare data
  plot_data <- data %>%
    mutate(
      neg_log_pval = -log10(get(pval_col)),
      significance = case_when(
        abs(get(fc_col)) > fc_threshold & get(pval_col) < pval_threshold & get(fc_col) > 0 ~ "Upregulated",
        abs(get(fc_col)) > fc_threshold & get(pval_col) < pval_threshold & get(fc_col) < 0 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Create plot
  ggplot(plot_data, aes_string(x = fc_col, y = "neg_log_pval")) +
    geom_point(aes(color = significance), alpha = 0.7) +
    scale_color_manual(values = c("Upregulated" = "red", 
                                 "Downregulated" = "blue", 
                                 "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    labs(x = "Log2 Fold Change",
         y = "-Log10 P-value",
         color = "Significance",
         title = "Volcano Plot") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# ANALYSIS WORKFLOW ===========================================================

#' Main analysis workflow
#' @param config_file Path to configuration file with file paths and parameters
run_integration_analysis <- function(config_file = "config.yaml") {
  
  cat("Starting Comprehensive Epigenome-Transcriptome Integration Analysis...\n")
  
  # Example workflow - adapt file paths as needed
  
  # 1. Load and process DEP data
  cat("Processing DEP data...\n")
  # deps_data <- process_dep_data("data/DEPs_data.txt")
  
  # 2. Load and process DMC data
  cat("Processing DMC data...\n")
  # dmcs_day1 <- process_unique_dmcs("data/promoters_day1.txt", "data/genes_day1.txt", "Day1")
  # dmcs_wk1 <- process_unique_dmcs("data/promoters_wk1.txt", "data/genes_wk1.txt", "Wk1")
  
  # 3. Integration analysis
  cat("Performing integration analysis...\n")
  # integrated_data <- inner_join(deps_data, dmcs_day1, by = "GeneList")
  
  # 4. Create visualizations
  cat("Generating comprehensive visualizations...\n")
  # correlation_plot <- create_correlation_plot(integrated_data, "methdiff_day1", "d1_avg", facet_by = "loc")
  # violin_plot <- create_violin_plot(integrated_data, "GeneList", "methdiff_day1")
  # expression_heatmap <- create_expression_heatmap(deps_matrix)
  
  # 5. GO enrichment analysis
  cat("Running GO enrichment analysis...\n")
  # significant_genes <- integrated_data %>% filter(abs(expression_fc) > 1) %>% pull(GeneList)
  # go_results <- run_go_enrichment(significant_genes)
  # go_plot <- create_go_plot(go_results)
  # go_barplot <- create_go_barplot(go_results$result)
  
  # 6. PCA analysis
  cat("Performing PCA analysis...\n")
  # pca_results <- perform_pca_analysis(expression_matrix)
  # pca_plot <- create_pca_plot(pca_results)
  # scree_plot <- create_scree_plot(pca_results)
  
  # 7. HOMER motif analysis
  cat("Processing HOMER motif analysis...\n")
  # homer_matrix <- process_homer_data(homer_data, c("day1", "week1"))
  # homer_heatmap <- create_homer_heatmap(homer_matrix)
  
  cat("Comprehensive analysis completed successfully!\n")
  
  # Return results list
  return(list(
    message = "Comprehensive analysis framework ready for your data",
    functions_available = c(
      # Data processing
      "process_dep_data", "process_unique_dmcs", "filter_significant_dmcs",
      # Basic visualizations
      "create_correlation_plot", "create_expression_heatmap", "create_volcano_plot",
      # Advanced visualizations  
      "create_violin_plot", "create_parallel_plot", "create_circular_barplot", "create_lollipop_plot",
      # GO enrichment
      "run_go_enrichment", "create_go_plot", "create_go_barplot",
      # PCA analysis
      "perform_pca_analysis", "create_pca_plot", "create_scree_plot", "get_top_pc_contributors",
      # HOMER analysis
      "process_homer_data", "create_homer_heatmap"
    )
  ))
}

#' Comprehensive pathway analysis workflow
#' @param gene_lists Named list of gene vectors for different conditions/timepoints
#' @param organism Organism for GO analysis
run_pathway_analysis <- function(gene_lists, organism = "hsapiens") {
  
  cat("Running comprehensive pathway analysis...\n")
  
  results <- list()
  
  for (condition in names(gene_lists)) {
    cat(paste("Analyzing", condition, "...\n"))
    
    # Run GO enrichment
    go_results <- run_go_enrichment(gene_lists[[condition]], organism = organism)
    
    # Create plots
    if (!is.null(go_results$result)) {
      go_plot <- create_go_plot(go_results)
      go_barplot <- create_go_barplot(go_results$result, title = paste("GO Enrichment -", condition))
      
      results[[condition]] <- list(
        go_results = go_results,
        go_plot = go_plot,
        go_barplot = go_barplot
      )
    }
  }
  
  return(results)
}

#' Multi-timepoint visualization workflow
#' @param data Integrated methylation and expression data
#' @param timepoints Vector of timepoint column names
#' @param gene_col Column name containing gene identifiers
run_temporal_analysis <- function(data, timepoints, gene_col = "GeneList") {
  
  cat("Creating temporal analysis visualizations...\n")
  
  results <- list()
  
  # Create violin plots for each timepoint
  for (tp in timepoints) {
    if (tp %in% names(data)) {
      violin_plot <- create_violin_plot(
        data = data,
        x_col = gene_col,
        y_col = tp,
        title = paste("Distribution at", tp)
      )
      results[[paste0("violin_", tp)]] <- violin_plot
    }
  }
  
  # Create parallel coordinate plot if multiple timepoints
  if (length(timepoints) > 1) {
    parallel_plot <- create_parallel_plot(
      data = data,
      value_cols = timepoints,
      group_col = gene_col,
      title = "Temporal Expression Patterns"
    )
    results[["parallel_plot"]] <- parallel_plot
  }
  
  # Create circular barplot for final timepoint
  if (length(timepoints) > 0) {
    final_tp <- timepoints[length(timepoints)]
    circular_plot <- create_circular_barplot(
      data = data,
      x_col = gene_col,
      y_col = final_tp,
      title = paste("Circular View -", final_tp)
    )
    results[["circular_plot"]] <- circular_plot
  }
  
  return(results)
}

# EXAMPLE USAGE ===============================================================

# # Set working directory
# # Run complete analysis
# results <- run_integration_analysis()
# 
# # BASIC ANALYSIS WORKFLOW:
# 
# # Load your data
# deps_data <- read_delim("data/your_deps_file.txt", delim = "\t")
# dmcs_data <- read_delim("data/your_dmcs_file.txt", delim = "\t")
# 
# # Create correlation plot
# correlation_plot <- create_correlation_plot(
#   data = your_integrated_data,
#   x_col = "methylation_diff",
#   y_col = "expression_fc",
#   facet_by = "genomic_region"
# )
# 
# # Save plot
# ggsave("results/correlation_plot.png", correlation_plot, width = 12, height = 8, dpi = 300)

# ADVANCED ANALYSIS EXAMPLES:

# # 1. GO ENRICHMENT ANALYSIS
# significant_genes <- c("EGFR", "TP53", "BRCA1", "MYC")  # Your gene list
# go_results <- run_go_enrichment(significant_genes, organism = "hsapiens")
# 
# # Create GO plots
# go_manhattan <- create_go_plot(go_results, highlight_terms = c("GO:0008283", "GO:0006915"))
# go_barplot <- create_go_barplot(go_results$result, top_n = 20)
# 
# # Save GO results
# ggsave("results/go_enrichment_plot.png", go_manhattan, width = 14, height = 10, dpi = 300)
# ggsave("results/go_barplot.png", go_barplot, width = 12, height = 8, dpi = 300)

# # 2. ADVANCED VISUALIZATIONS
# 
# # Violin plot for methylation distribution
# violin_plot <- create_violin_plot(
#   data = methylation_data,
#   x_col = "GeneList",
#   y_col = "methdiff_day1",
#   orientation = "horizontal",
#   title = "Methylation Distribution by Gene"
# )
# 
# # Parallel coordinate plot for temporal patterns
# parallel_plot <- create_parallel_plot(
#   data = temporal_data,
#   value_cols = c("nsc_methylation", "d1_methylation", "wk1_methylation"),
#   group_col = "Gene",
#   title = "Temporal Methylation Patterns"
# )
# 
# # Circular barplot
# circular_plot <- create_circular_barplot(
#   data = expression_data,
#   x_col = "GeneList",
#   y_col = "fold_change",
#   title = "Expression Changes - Circular View"
# )
# 
# # Lollipop plot
# lollipop_plot <- create_lollipop_plot(
#   data = differential_data,
#   x_col = "GeneList",
#   y_col = "log2FC",
#   color_col = "significance"
# )

# # 3. PCA ANALYSIS
# 
# # Prepare expression matrix (genes as rows, samples as columns)
# expression_matrix <- deps_data %>%
#   column_to_rownames("GeneList") %>%
#   select(nsc_avg, d1_avg, wk1_avg) %>%
#   as.matrix()
# 
# # Perform PCA
# pca_results <- perform_pca_analysis(expression_matrix, scale = TRUE)
# 
# # Create PCA visualizations
# pca_plot <- create_pca_plot(
#   pca_results = pca_results,
#   sample_groups = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
#   title = "PCA of Gene Expression"
# )
# 
# scree_plot <- create_scree_plot(pca_results, n_components = 10)
# 
# # Get top contributing genes for PC1
# top_contributors <- get_top_pc_contributors(pca_results, pc_number = 1, n_features = 20)
# print(top_contributors)
# 
# # Save PCA plots
# ggsave("results/pca_plot.png", pca_plot, width = 10, height = 8, dpi = 300)
# ggsave("results/scree_plot.png", scree_plot, width = 10, height = 6, dpi = 300)

# # 4. HOMER MOTIF ANALYSIS
# 
# # Process HOMER data
# homer_matrix_hypo <- process_homer_data(
#   homer_data = homer_hypo_data,
#   value_cols = c("Day1_pvalue", "Week1_pvalue"),
#   term_col = "Motif"
# )
# 
# homer_matrix_hyper <- process_homer_data(
#   homer_data = homer_hyper_data,
#   value_cols = c("Day1_pvalue", "Week1_pvalue", "Combined_pvalue"),
#   term_col = "Motif"
# )
# 
# # Create HOMER heatmaps
# homer_heatmap_hypo <- create_homer_heatmap(
#   homer_matrix = homer_matrix_hypo,
#   heat_limits = c(40, 510),
#   title = "HOMER Analysis - Hypomethylated Regions",
#   color_scheme = "grey"
# )
# 
# homer_heatmap_hyper <- create_homer_heatmap(
#   homer_matrix = homer_matrix_hyper,
#   heat_limits = c(5, 150),
#   title = "HOMER Analysis - Hypermethylated Regions",
#   color_scheme = "viridis"
# )

# # 5. COMPREHENSIVE WORKFLOW EXAMPLES
# 
# # Multi-condition pathway analysis
# gene_lists <- list(
#   "Day1_Hyper" = day1_hyper_genes,
#   "Day1_Hypo" = day1_hypo_genes,
#   "Week1_Hyper" = week1_hyper_genes,
#   "Week1_Hypo" = week1_hypo_genes
# )
# 
# pathway_results <- run_pathway_analysis(gene_lists, organism = "hsapiens")
# 
# # Temporal analysis workflow
# temporal_results <- run_temporal_analysis(
#   data = integrated_temporal_data,
#   timepoints = c("nsc_avg", "d1_avg", "wk1_avg"),
#   gene_col = "GeneList"
# )
# 
# # Save all temporal plots
# for (plot_name in names(temporal_results)) {
#   ggsave(paste0("results/", plot_name, ".png"), 
#          temporal_results[[plot_name]], 
#          width = 12, height = 8, dpi = 300)
# }

# # 6. PUBLICATION-READY FIGURE CREATION
# 
# # Create a multi-panel figure
# library(cowplot)
# 
# # Combine plots into publication figure
# figure_1 <- plot_grid(
#   correlation_plot, pca_plot,
#   violin_plot, go_barplot,
#   labels = c("A", "B", "C", "D"),
#   ncol = 2
# )
# 
# # Save publication figure
# ggsave("results/publication_figure_1.png", figure_1, 
#        width = 16, height = 12, dpi = 300)

# Print session info for reproducibility
print("Session Info:")
sessionInfo()