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
sample_table$genotype <- factor(
  sample_table$genotype,
  levels = c("GFP", "Drp-1 CA")
)

sample_table$treatment <- factor(
  sample_table$treatment,
  levels = c("Control", "Insulin")
)


# Build DESeq2 dataset.
dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ genotype * treatment
)


# Run DESeq2 pipeline.
dds <- DESeq(dds)


#Pre-DGE QC plots---------------------------------------------------------------
save_tiff(
  "results/qc_plots/multifactor/total_counts_plot.tiff",
  plot_code = function() {
    print(plot_total_counts(dds))
  }
)

save_tiff(
  "results/qc_plots/multifactor/library_complexity.tiff",
  plot_code = function() {
    print(plot_library_complexity(dds))
  }
)

save_tiff(
  "results/qc_plots/multifactor/gene_detection_plot.tiff",
  plot_code = function() {
    print(plot_gene_detection(dds))
  }
)

save_tiff(
  "results/qc_plots/multifactor/dispersion_estimate_plot.tiff",
  height = 150,
  plot_code = function() {
    plotDispEsts(dds)
  }
)


# Calculate contrasts-----------------------------------------------------------
genotype_effect <- lfcShrink(
  dds,
  coef = "genotype_Drp.1.CA_vs_GFP",
  type = "apeglm"
)

treatment_effect <- lfcShrink(
  dds,
  coef = "treatment_Insulin_vs_Control",
  type = "apeglm"
)

interaction_effect <- lfcShrink(
  dds,
  coef = "genotypeDrp.1.CA.treatmentInsulin",
  type = "apeglm"
)


# Post-DGE QC Plots-------------------------------------------------------------
save_tiff(
  "results/qc_plots/multifactor/genotype_effect_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      genotype_effect$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)

save_tiff(
  "results/qc_plots/multifactor/treatment_effect_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      treatment_effect$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)

save_tiff(
  "results/qc_plots/multifactor/interaction_effect_p_value_histogram.tiff",
  plot_code = function() {
    hist(
      interaction_effect$pvalue,
      breaks = 20,
      col = "grey",
      main = "p-value Distribution",
      xlab = "p-value"
    )
  }
)


# Save files--------------------------------------------------------------------
# Save dds for heat maps, PCA plots, and single gene plots.
saveRDS(dds, "data/rds/multifactor_dds.rds")


# Obtain background genes for ORA.
background_genes <- sub("^GeneID:", "", rownames(dds))
saveRDS(
  background_genes,
  "data/rds/multifactor_background_genes.rds"
)


# Save contrasts for MA plots and volcano plots.
saveRDS(genotype_effect, "data/rds/genotype_effect.rds")
saveRDS(treatment_effect, "data/rds/treatment_effect.rds")
saveRDS(interaction_effect, "data/rds/interaction_effect.rds")


# Save contrasts to excel.
save_excel(genotype_effect, "genotype_effect")
save_excel(treatment_effect, "treatment_effect")
save_excel(interaction_effect, "interaction_effect")


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/08_deseq2_multifactor.txt"
)
