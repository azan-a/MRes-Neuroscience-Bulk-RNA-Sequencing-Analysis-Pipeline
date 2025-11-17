#Configuration------------------------------------------------------------------
library(sessioninfo)
library(EnhancedVolcano)
library(org.Rn.eg.db)
library(dplyr)


# Get functions.
source("scripts/plotting_utils.R")


# Set parameters.
lfc_cutoff <- 0.26
alpha <- 0.05
n_top_genes <- 30


table_name <- "genotype_effect"
design_formula <- "multifactor"

# Get contrast results.
deg_table <- readRDS(paste0("data/rds/", table_name, ".rds"))


#Convert to gene symbol---------------------------------------------------------
# Convert deg_table into a dataframe.
deg_table <- as.data.frame(deg_table)

# Map NCBI gene ID to gene symbol.
gene_ids <- sub("GeneID:", "", rownames(deg_table))
gene_symbols <- mapIds(
  org.Rn.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multivals = "first"
)

# Use gene symbols where available, fall back on ID otherwise.
gene_symbols <- ifelse(
  is.na(gene_symbols),
  gene_ids,
  gene_symbols
)

rownames(deg_table) <- gene_symbols

# Colour-label genes------------------------------------------------------------
# For Non-significant genes.
keyvals <- rep("#808080", nrow(deg_table))

# For up-regulated genes.
keyvals[
  deg_table$log2FoldChange >= lfc_cutoff & deg_table$padj <= alpha
] <- "#E69F00"

# For down-regulated genes.
keyvals[
  deg_table$log2FoldChange <= -lfc_cutoff & deg_table$padj <= alpha
] <- "#56B4E9"


# Add labels if you want legends.
names(keyvals) <- rep("NS", nrow(deg_table))
names(keyvals)[
  deg_table$log2FoldChange >= lfc_cutoff & deg_table$padj <= alpha
] <- "Upregulated"
names(keyvals)[
  deg_table$log2FoldChange <= -lfc_cutoff & deg_table$padj <= alpha
] <- "Downregulated"


#Choose gene label--------------------------------------------------------------
# Either label top differentially expressed genes.
top_genes <- deg_table %>%
  dplyr::filter(
    !is.na(padj),
    padj <= alpha,
    abs(log2FoldChange) >= lfc_cutoff
  ) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = n_top_genes) %>%
  rownames(deg_table)

top_genes <- grep("^LOC", top_genes, invert = TRUE, value = TRUE)


# Or label a custom list of genes.
custom_gene_list <- readLines("config/custom_gene_list.txt")


#Volcano plot-------------------------------------------------------------------
# Build Volcano plot.
volcano_plot <- EnhancedVolcano(
  deg_table,
  lab = rownames(deg_table),
  x = "log2FoldChange",
  y = "padj",
  selectLab = top_genes, # Change for top_genes or custom_genes_list.
  ylab = bquote(~ -Log[10] ~ italic(P[adj])),
  title = NULL,
  labSize = 4,
  subtitle = NULL,
  caption = NULL,
  cutoffLineType = "blank",
  labFace = 'bold',
  colAlpha = 0.5,
  colCustom = keyvals,
  legendPosition = "none",
  drawConnectors = TRUE,
  colConnectors = "grey60",
  max.overlaps = Inf,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)


# Check plot before saving.
print(volcano_plot)


# Save plot.
save_tiff(
  paste0(
    "results/figures/",
    design_formula,
    "/",
    table_name,
    "_volcano_plot_top_genes.tiff"
  ),
  width = 200,
  height = 150,
  plot_code = function() {
    print(volcano_plot)
  }
)


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/12_volcano_plot.txt"
)
