# Configuration-----------------------------------------------------------------
library(sessioninfo)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(RNAseqQC)
library(apeglm)


# Get functions
source("scripts/deseq2_utils.R")
source("scripts/plotting_utils.R")


txi <- readRDS("data/rds/txi.rds")
sample_table <- readxl::read_excel("config/sample_table.xlsx")


# DESeq2 Analysis---------------------------------------------------------------
# Set factor levels.
sample_table$condition <- relevel(
  factor(sample_table$condition),
  ref = "GFP + Control"
)


# Build DESeq2 dataset.
dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~condition
)


# Run DESeq2 pipeline.
dds <- DESeq(dds)


#Pre-DGE QC plots---------------------------------------------------------------
save_tiff(
  "results/qc_plots/singlefactor/total_counts_plot.tiff",
  plot_code = function() {
    print(plot_total_counts(dds))
  }
)

save_tiff(
  "results/qc_plots/singlefactor/library_complexity.tiff",
  plot_code = function() {
    print(plot_library_complexity(dds))
  }
)

save_tiff(
  "results/qc_plots/singlefactor/gene_detection_plot.tiff",
  plot_code = function() {
    print(plot_gene_detection(dds))
  }
)

save_tiff(
  "results/qc_plots/singlefactor/dispersion_estimate_plot.tiff",
  height = 150,
  plot_code = function() {
    plotDispEsts(dds)
  }
)


# Calculate contrasts-----------------------------------------------------------
control_vs_drp1ca <- lfcShrink(
  dds,
  coef = "condition_Drp.1.CA...Control_vs_GFP...Control",
  type = "apeglm"
)

control_vs_insulin <- lfcShrink(
  dds,
  coef = "condition_GFP...Insulin_vs_GFP...Control",
  type = "apeglm"
)


dds$condition <- relevel(sample_table$condition, ref = "GFP + Insulin")
dds <- nbinomWaldTest(dds)
insulin_vs_drp1ca_insulin <- lfcShrink(
  dds,
  coef = "condition_Drp.1.CA...Insulin_vs_GFP...Insulin",
  type = "apeglm"
)

dds$condition <- relevel(sample_table$condition, ref = "Drp-1 CA + Control")
dds <- nbinomWaldTest(dds)
drp1ca_vs_drp1ca_insulin <- lfcShrink(
  dds,
  coef = "condition_Drp.1.CA...Insulin_vs_Drp.1.CA...Control",
  type = "apeglm"
)

# Post-DGE QC Plots-------------------------------------------------------------
save_tiff(
  "results/qc_plots/singlefactor/control_vs_drp1ca_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      control_vs_drp1ca$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)

save_tiff(
  "results/qc_plots/singlefactor/control_vs_insulin_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      control_vs_insulin$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)

save_tiff(
  "results/qc_plots/singlefactor/insulin_vs_drp1ca_insulin_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      insulin_vs_drp1ca_insulin$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)

save_tiff(
  "results/qc_plots/singlefactor/drp1ca_vs_drp1ca_insulin_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      drp1ca_vs_drp1ca_insulin$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)

# Save files--------------------------------------------------------------------
# Save dds for heat maps, PCA plots, and single gene plots.
saveRDS(dds, "data/rds/singlefactor_dds.rds")


# Obtain background genes for ORA.
background_genes <- sub("^GeneID:", "", rownames(dds))
saveRDS(
  background_genes,
  "data/rds/singlefactor_background_genes.rds"
)


# Save contrasts for MA plots and volcano plots.
saveRDS(control_vs_drp1ca, "data/rds/control_vs_drp1ca.rds")
saveRDS(control_vs_insulin, "data/rds/control_vs_insulin.rds")
saveRDS(insulin_vs_drp1ca_insulin, "data/rds/insulin_vs_drp1ca_insulin.rds")
saveRDS(drp1ca_vs_drp1ca_insulin, "data/rds/drp1ca_vs_drp1ca_insulin.rds")


# Save contrasts to excel.
save_excel(control_vs_drp1ca, "control_vs_drp1ca")
save_excel(control_vs_insulin, "control_vs_insulin")
save_excel(insulin_vs_drp1ca_insulin, "insulin_vs_drp1ca_insulin")
save_excel(drp1ca_vs_drp1ca_insulin, "drp1ca_vs_drp1ca_insulin")


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/08_deseq2_singlefactor.txt"
)
