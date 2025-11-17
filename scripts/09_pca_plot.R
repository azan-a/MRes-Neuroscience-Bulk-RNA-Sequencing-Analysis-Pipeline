#Configuration------------------------------------------------------------------
library(sessioninfo)
library(DESeq2)
library(ggplot2)


design_formula <- "multifactor"


# Get functions.
source("scripts/plotting_utils.R")
# Get dds results.
dds <- readRDS(paste0("data/rds/", design_formula, "_dds.rds"))


# r-log transformation----------------------------------------------------------
# Calculate r-log transformation.
rld <- rlog(dds, blind = TRUE)


#QC plot.
save_tiff(
  paste0("results/qc_plots/", design_formula, "/mean-sd_plot_blind_true.tiff"),
  height = 150,
  plot_code = function() {
    vsn::meanSdPlot(assay(rld))
  }
)

# Prepare for plotting---------------------------------------------------------
# Extract PCA data,
pca_data <- plotPCA(
  rld,
  intgroup = c("treatment", "genotype"),
  returnData = TRUE
)
pc_percent <- round(100 * attr(pca_data, "percentVar"))


# Define palette.
pca_palette <- c("Control" = "#009E73", "Insulin" = "#56B4E9")


# PCA Plot----------------------------------------------------------------------
# Build PCA plot.
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, fill = treatment, shape = genotype)
) +
  geom_point(size = 3, stroke = 0.5, alpha = 0.7) +
  theme_light(base_size = 15) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  xlab(paste0("PC1: ", pc_percent[1], "% variance")) +
  ylab(paste0("PC2: ", pc_percent[2], "% variance")) +
  scale_y_continuous(limits = c(-10, 10), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(-15, 15), expand = expansion(mult = c(0, 0))) +
  coord_fixed() +
  scale_fill_manual(values = pca_palette) +
  scale_shape_manual(values = c("GFP" = 21, "Drp-1 CA" = 24)) +
  labs(fill = "Treatment", shape = "Viral modification") +
  guides(
    fill = guide_legend(override.aes = list(color = pca_palette))
  ) +
  remove_grid()


#Check plot.
print(pca_plot)


# Save plot.
save_tiff(
  paste0("results/figures/", design_formula, "/pca_plot.tiff"),
  plot_code = function() {
    print(pca_plot)
  }
)


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/09_pca_plot.txt"
)
