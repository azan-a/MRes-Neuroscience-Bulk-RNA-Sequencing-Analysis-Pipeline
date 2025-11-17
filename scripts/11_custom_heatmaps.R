#Configurations-----------------------------------------------------------------
library(sessioninfo)
library(ComplexHeatmap)
library(org.Rn.eg.db)
library(SummarizedExperiment)


design_formula <- "multifactor"


# Get functions.
source("scripts/plotting_utils.R")
# Get dds results.
dds <- readRDS(paste0("data/rds/", design_formula, "_dds.rds"))
# Get custom gene list
custom_gene_list <- readLines("config/custom_gene_list.txt")


# Parameters
colour_steps <- 50

# r-log transformation----------------------------------------------------------
# Perform r-log transformation.
rld <- DESeq2::rlog(dds, blind = FALSE)


# QC plot.
save_tiff(
  paste0("results/qc_plots/", design_formula, "/mean-sd_plot_blind_false.tiff"),
  height = 150,
  plot_code = function() {
    vsn::meanSdPlot(assay(rld))
  }
)

# Prepare for plotting----------------------------------------------------------
# Add annotations.
anno <- as.data.frame(colData(rld)[, c("treatment", "genotype")])
colnames(anno) <- c("Treatment", "Viral Modification")

annotation_colours <- list(
  Treatment = c("Control" = "#009E73", "Insulin" = "#56B4E9"),
  "Viral Modification" = c("GFP" = "#CC79A7", "Drp-1 CA" = "#0072B2")
)

# Map to gene symbols.
gene_ids <- gsub("GeneID:", "", rownames(rld))
gene_symbols <- mapIds(
  org.Rn.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Replace rownames in rld with gene symbols.
rownames(rld) <- gene_symbols

# Subset assay.
mat <- assay(rld)[custom_gene_list, ]

# Compute row Z-scores.
mat_scaled <- t(scale(t(mat)))


# Custom heat map---------------------------------------------------------------
# Build Z-score scaled heat map
custom_heatmap <- pheatmap(
  mat_scaled,
  annotation_col = anno,
  name = "Z-Score",
  color = colorRampPalette(c("#4A6FE3", "#E2E2E2", "#D33F6A"))(colour_steps),
  annotation_colors = annotation_colours
)


#Check plot.
print(custom_heatmap)


#Save plot
save_tiff(
  paste0("results/figures/", design_formula, "/custom_heatmap.tiff"),
  height = 150,
  width = 200,
  plot_code = function() {
    print(custom_heatmap)
  }
)


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/11_custom_heat_maps.txt"
)
