# Epigenome-Transcriptome Integration Analysis

A comprehensive R package for integrating DNA methylation (RRBS) and gene expression data to identify coordinated epigenetic and transcriptional changes across multiple time points and genomic regions.

## 🧬 Overview

This project provides a complete analytical framework for:
- **Differential Methylation Analysis**: Identification and characterization of differentially methylated cytosines (DMCs)
- **Gene Expression Analysis**: Processing of differentially expressed proteins/genes (DEPs)
- **Integration Analysis**: Correlation analysis between methylation and expression changes
- **GO/Pathway Enrichment**: Comprehensive pathway analysis using gprofiler2
- **Advanced Visualizations**: Violin plots, parallel coordinates, circular plots, lollipop plots
- **PCA Analysis**: Principal component analysis with contribution analysis
- **HOMER Motif Analysis**: Transcription factor binding site enrichment visualization
- **Publication-Ready Figures**: High-quality, customizable visualizations

## 🚀 Features

- ✅ Modular, reusable functions for bioinformatics analysis
- ✅ Comprehensive visualization suite with 10+ plot types
- ✅ Support for multiple time points and genomic regions (promoters, gene bodies)  
- ✅ Statistical correlation analysis with significance testing
- ✅ Automated GO/KEGG pathway enrichment analysis
- ✅ PCA analysis with variance decomposition and feature contribution
- ✅ HOMER motif analysis integration
- ✅ Advanced plot types (violin, parallel coordinate, circular, lollipop)
- ✅ Temporal analysis workflows
- ✅ Multi-condition comparative analysis
- ✅ Customizable color palettes and themes
- ✅ Publication-ready figure generation
- ✅ Reproducible analysis workflow

## 📋 Requirements

### R Version
- R ≥ 4.0.0

### Required Packages
```r
# CRAN packages
packages <- c(
  "ggsci", "devtools", "hrbrthemes", "ggplot2", "ggpubr", "dplyr", "ggrepel",
  "tidyverse", "readxl", "readr", "tibble", "reshape2", "gprofiler2",
  "RColorBrewer", "scales", "sm", "viridis", "writexl", "GGally",
  "pheatmap", "ComplexHeatmap", "circlize", "colorspace", "heatmaply",
  "superheat", "lattice", "gplots", "paletteer", "kableExtra"
)

# Install all packages
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize"))

# GitHub packages
devtools::install_github("kassambara/ggpubr")
```

## 🗂️ Project Structure

```
epigenome-transcriptome-analysis/
├── analysis.R                 # Main analysis script
├── README.md                  # This file
├── data/                      # Input data directory
│   ├── DEPs_data.txt         # Gene expression data
│   ├── DMCs_promoters.txt    # Promoter methylation data
│   └── DMCs_genes.txt        # Gene body methylation data
├── results/                   # Output directory
│   ├── plots/                # Generated visualizations
│   └── tables/               # Analysis results tables
└── config/                   # Configuration files
    └── analysis_config.yaml  # Analysis parameters
```

## 🔧 Usage

### Quick Start

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/epigenome-transcriptome-analysis.git
cd epigenome-transcriptome-analysis
```

2. **Prepare your data**
Place your data files in the `data/` directory:
- Gene expression data (DEPs): Tab-delimited file with GeneList and sample columns
- Methylation data (DMCs): Tab-delimited files with methylation differences and significance values

3. **Run the analysis**
```r
source("analysis.R")

# Run complete workflow
results <- run_integration_analysis()

# Or use individual functions
correlation_plot <- create_correlation_plot(
  data = your_data,
  x_col = "methylation_diff", 
  y_col = "expression_fc",
  facet_by = "genomic_region"
)
```

### Key Functions

#### Data Processing
- `process_dep_data()`: Load and process gene expression data
- `process_unique_dmcs()`: Process methylation data for specific time points
- `filter_significant_dmcs()`: Filter DMCs by significance thresholds
- `process_homer_data()`: Process HOMER motif analysis results

#### GO Enrichment & Pathway Analysis
- `run_go_enrichment()`: Perform GO/KEGG enrichment using gprofiler2
- `create_go_plot()`: Create publication-ready GO Manhattan plots
- `run_pathway_analysis()`: Multi-condition pathway analysis workflow

#### Basic Visualizations
- `create_correlation_plot()`: Methylation vs expression correlation plots
- `create_expression_heatmap()`: Clustered heatmaps with k-means clustering
- `create_go_barplot()`: Gene Ontology enrichment bar plots
- `create_volcano_plot()`: Volcano plots for differential analysis

#### Advanced Visualizations
- `create_violin_plot()`: Violin plots with jitter points for distribution analysis
- `create_parallel_plot()`: Parallel coordinate plots for temporal patterns
- `create_circular_barplot()`: Circular bar plots for categorical data
- `create_lollipop_plot()`: Lollipop plots for ranked data
- `create_homer_heatmap()`: HOMER motif analysis heatmaps

#### PCA Analysis
- `perform_pca_analysis()`: Principal component analysis with variance calculation
- `create_pca_plot()`: PCA scatter plots with sample grouping
- `create_scree_plot()`: Scree plots showing variance explained
- `get_top_pc_contributors()`: Identify top contributing features for each PC

#### Workflow Functions
- `run_integration_analysis()`: Complete analysis pipeline
- `run_temporal_analysis()`: Multi-timepoint visualization workflow
- `ipak()`: Package installation and loading utility

#### Utilities
- `create_color_palette()`: Generate custom color palettes
- `load_methylation_data()`: Standardized data loading

## 📊 Example Analyses

### 1. Basic Correlation Analysis
```r
# Create correlation plot between methylation and expression
plot <- create_correlation_plot(
  data = integrated_data,
  x_col = "methdiff_day1",
  y_col = "expression_log2fc",
  facet_by = "genomic_location"
)

ggsave("results/plots/correlation_analysis.png", plot, width = 12, height = 8)
```

### 2. GO Enrichment Analysis
```r
# Run GO enrichment analysis
significant_genes <- c("EGFR", "TP53", "BRCA1", "MYC")
go_results <- run_go_enrichment(significant_genes, organism = "hsapiens")

# Create visualization
go_plot <- create_go_plot(go_results, highlight_terms = c("GO:0008283"))
go_barplot <- create_go_barplot(go_results$result, top_n = 20)
```

### 3. Advanced Visualizations
```r
# Violin plot for methylation distribution
violin_plot <- create_violin_plot(
  data = methylation_data,
  x_col = "GeneList", 
  y_col = "methdiff_day1",
  orientation = "horizontal"
)

# Parallel coordinate plot for temporal patterns  
parallel_plot <- create_parallel_plot(
  data = temporal_data,
  value_cols = c("nsc_methylation", "d1_methylation", "wk1_methylation"),
  group_col = "Gene"
)

# Circular barplot
circular_plot <- create_circular_barplot(
  data = expression_data,
  x_col = "GeneList",
  y_col = "fold_change"
)
```

### 4. PCA Analysis  
```r
# Perform PCA on expression data
expression_matrix <- deps_data %>%
  column_to_rownames("GeneList") %>%
  select(nsc_avg, d1_avg, wk1_avg) %>%
  as.matrix()

pca_results <- perform_pca_analysis(expression_matrix, scale = TRUE)

# Create PCA visualizations
pca_plot <- create_pca_plot(pca_results, sample_groups = c("Control", "Treatment"))
scree_plot <- create_scree_plot(pca_results)

# Get top contributing genes
top_contributors <- get_top_pc_contributors(pca_results, pc_number = 1, n_features = 20)
```

### 5. HOMER Motif Analysis
```r
# Process HOMER results
homer_matrix <- process_homer_data(
  homer_data = homer_results,
  value_cols = c("Day1_pvalue", "Week1_pvalue"),
  term_col = "Motif"
)

# Create heatmap
homer_heatmap <- create_homer_heatmap(
  homer_matrix = homer_matrix,
  heat_limits = c(5, 150),
  title = "HOMER Motif Enrichment"
)
```

### 6. Multi-Condition Workflow
```r
# Analyze multiple conditions simultaneously
gene_lists <- list(
  "Day1_Hyper" = day1_hyper_genes,
  "Day1_Hypo" = day1_hypo_genes,
  "Week1_Hyper" = week1_hyper_genes,
  "Week1_Hypo" = week1_hypo_genes
)

pathway_results <- run_pathway_analysis(gene_lists)

# Temporal analysis across timepoints
temporal_results <- run_temporal_analysis(
  data = integrated_data,
  timepoints = c("nsc_avg", "d1_avg", "wk1_avg")
)
```

## 📈 Input Data Format

### Gene Expression Data (DEPs)
```
GeneList    nsc_avg    d1_avg    wk1_avg    
Gene1       2.5        1.2       -0.8
Gene2      -1.1        0.3        2.1
...
```

### Methylation Data (DMCs)
```
GeneList    start    end    methdiff_day1    qvalue_day1    location
Gene1       1000     1200   25.5             0.001          Promoters
Gene2       5000     5200   -18.2            0.005          Gene_bodies
...
```

## 🎨 Visualization Gallery

The package generates publication-ready visualizations including:

### Basic Plots
- **Correlation Plots**: Scatter plots showing methylation vs expression relationships with regression lines
- **Heatmaps**: Clustered expression patterns across time points with k-means clustering
- **Volcano Plots**: Differential analysis results with significance highlighting
- **Bar Plots**: Pathway and GO term enrichment visualizations

### Advanced Visualizations  
- **Violin Plots**: Distribution analysis with jitter points and statistical summaries
- **Parallel Coordinate Plots**: Multi-dimensional temporal pattern visualization
- **Circular Bar Plots**: Polar coordinate system for categorical data visualization
- **Lollipop Plots**: Ranked data with baseline comparisons
- **PCA Plots**: Principal component analysis with variance explanation and sample grouping
- **Scree Plots**: Variance explained by each principal component
- **HOMER Heatmaps**: Transcription factor binding motif enrichment patterns with customizable color schemes

### Specialized Analysis Plots
- **GO Manhattan Plots**: Genome-wide association style plots for pathway enrichment
- **Multi-panel Figures**: Combined publication-ready figures with multiple plot types
- **Temporal Series**: Time-course analysis visualizations
- **Comparative Analysis**: Side-by-side condition comparisons

All plots feature:
- Customizable color palettes (viridis, RColorBrewer, custom schemes)
- Publication-ready themes and formatting
- High-resolution output (300+ DPI)
- Flexible sizing and aspect ratios
- Statistical annotations and significance indicators


---
**Keywords**: bioinformatics, epigenetics, transcriptomics, DNA methylation, RRBS, gene expression, R, data integration, GO enrichment, PCA analysis, HOMER motifs, pathway analysis
