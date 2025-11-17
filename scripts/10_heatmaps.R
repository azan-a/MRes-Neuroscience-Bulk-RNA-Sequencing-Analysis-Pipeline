# Configuration-----------------------------------------------------------------
library(sessioninfo)
library(org.Rn.eg.db)
library(ComplexHeatmap)
library(SummarizedExperiment)


design_formula <- "singlefactor"


# Get functions.
source("scripts/plotting_utils.R")
# Get dds results.
dds <- readRDS(paste0("data/rds/", design_formula, "_dds.rds"))


# Parameters
n_top_genes <- 50
colour_steps <- 50


# r-log transformation----------------------------------------------------------
# Perform r-log transformation.
rld <- DESeq2::rlog(dds, blind = TRUE)

# Prepare for plotting----------------------------------------------------------

# Add annotations.
anno <- as.data.frame(colData(rld)[, c("treatment", "genotype")])
colnames(anno) <- c("Treatment", "Viral Modification")

annotation_colours <- list(
  Treatment = c("Control" = "#009E73", "Insulin" = "#56B4E9"),
  "Viral Modification" = c("GFP" = "#CC79A7", "Drp-1 CA" = "#0072B2")
)


# Extract top variable gene names.
top_variable_genes <- head(
  order(genefilter::rowVars(assay(rld)), decreasing = TRUE),
  n_top_genes
)


# Extract normalised counts for top variable genes.
mat <- assay(rld)[top_variable_genes, ]


# Map gene IDs to symbols.
gene_ids <- sub("GeneID:", "", rownames(mat))
gene_symbols <- mapIds(
  org.Rn.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
rownames(mat) <- gene_symbols


# Non-scaled heat map-----------------------------------------------------------
# Build non-scaled heat map.
nonscaled_heatmap <- pheatmap(
  mat,
  annotation_col = anno,
  name = " ",
  color = colorRampPalette(c("#F9F9F9", "#79ABE2", "#00366C"))(colour_steps),
  annotation_colors = annotation_colours
)


# Check plot.
print(nonscaled_heatmap)


# Save plot.
save_tiff(
  paste0("results/figures/", design_formula, "/non_scaled_heatmap.tiff"),
  height = 250,
  plot_code = function() {
    print(nonscaled_heatmap)
  }
)

# Z-score scaled heat map-------------------------------------------------------
# Calculate Z-score.
mat_scaled <- t(scale(t(mat)))


# Build Z-score scaled heat map
scaled_heatmap <- pheatmap(
  mat_scaled,
  annotation_col = anno,
  name = "Z-Score",
  color = colorRampPalette(c("#4A6FE3", "#E2E2E2", "#D33F6A"))(colour_steps),
  annotation_colors = annotation_colours
)


#Check plot.
print(scaled_heatmap)


# Save plot.
save_tiff(
  paste0("results/figures/", design_formula, "/scaled_heatmap.tiff"),
  height = 250,
  plot_code = function() {
    print(scaled_heatmap)
  }
)


#Sample-distance heat map-------------------------------------------------------
# Calculate sample distance.
sample_dist <- dist(t(assay(rld)))

# Convert into matrix
sample_dist_matrix <- as.matrix(sample_dist)

sample_dist_heatmap <- pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  color = colorRampPalette(c("#F9F9F9", "#79ABE2", "#00366C"))(colour_steps),
  name = " ",
  annotation_colors = annotation_colours
)


# Check plot.
print(sample_dist_heatmap)


# Save plot.
save_tiff(
  paste0("results/figures/", design_formula, "/sample_distance_heatmap.tiff"),
  height = 150,
  width = 200,
  plot_code = function() {
    print(sample_dist_heatmap)
  }
)


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/10_heat_maps.txt"
)
