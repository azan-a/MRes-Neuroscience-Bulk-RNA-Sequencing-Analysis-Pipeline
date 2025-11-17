#Configuration------------------------------------------------------------------
library(sessioninfo)
library(ggplot2)

# Get functions.
source("scripts/plotting_utils.R")

# Set parameters.
lfc_cutoff <- 0.26
alpha <- 0.05

design_formula <- "multifactor"
contrast_name <- "treatment_effect"

# Get contrast results.
contrast_result <- readRDS(paste0("data/rds/", contrast_name, ".rds"))
contrast_df <- as.data.frame(contrast_result)
contrast_df <- contrast_df[!is.na(contrast_df$padj), ]


# Label genes-------------------------------------------------------------------
# For Non-significant genes.
keyvals <- rep("NS", nrow(contrast_df))

# For up-regulated genes.
keyvals[
  contrast_df$log2FoldChange >= lfc_cutoff & contrast_df$padj <= alpha
] <- "Upregulated"

# For down-regulated genes.
keyvals[
  contrast_df$log2FoldChange <= -lfc_cutoff & contrast_df$padj <= alpha
] <- "Downregulated"

#Add to contrast_df
contrast_df$colour_group <- keyvals


# MA plot-----------------------------------------------------------------------
# Build MA plot.
ma_plot <- ggplot(
  contrast_df,
  aes(x = baseMean, y = log2FoldChange, color = colour_group)
) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "NS" = "#808080",
      "Upregulated" = "#E69F00",
      "Downregulated" = "#56B4E9"
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = lfc_cutoff, linetype = "dotted", color = "black") +
  geom_hline(yintercept = -lfc_cutoff, linetype = "dotted", color = "black") +
  labs(
    x = expression("Mean of normalised counts (log"[10] * ")"),
    y = "Log2 Fold Change"
  ) +
  theme_classic() +
  theme(legend.position = "none")


# Check plot.
print(ma_plot)


# Save plot.
save_tiff(
  paste0(
    "results/figures/",
    design_formula,
    "/",
    contrast_name,
    "_ma_plot.tiff"
  ),
  plot_code = function() {
    print(ma_plot)
  }
)


# Save session info
writeLines(capture.output(session_info()), "results/sessioninfo/13_ma_plot.txt")
