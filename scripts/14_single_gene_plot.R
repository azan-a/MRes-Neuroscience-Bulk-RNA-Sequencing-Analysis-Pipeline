#Configuration------------------------------------------------------------------
library(sessioninfo)
library(ggplot2)
library(dplyr)


design_formula <- "multifactor"


# Get functions.
source("scripts/plotting_utils.R")
# Get dds results.
dds <- readRDS(paste0("data/rds/", design_formula, "_dds.rds"))

gene <- "GeneID:89829"
gene_name <- "Socs3"

# Prepare for plotting----------------------------------------------------------
# Extract counts for gene.
gene_of_interest <- DESeq2::plotCounts(
  dds,
  gene,
  returnData = TRUE,
  intgroup = c("genotype", "treatment")
)


# Obtain summary data.
summary_data <- gene_of_interest %>%
  group_by(genotype, treatment) %>%
  summarise(
    mean = mean(count),
    SEM = sd(count) / sqrt(n())
  )


#Add colour to bar chart.
bar_fill_colours <- c("Control" = "#009E73", "Insulin" = "#56B4E9")
bar_outline_colours <- c("Control" = "#005841", "Insulin" = "#2A6E98")


# Single gene plot--------------------------------------------------------------
# Build single gene plot.
single_gene_plot <- ggplot(
  summary_data,
  aes(x = genotype, y = mean, fill = treatment)
) +
  geom_bar(
    stat = "identity",
    width = 0.8,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  geom_errorbar(
    aes(ymin = mean - SEM, ymax = mean + SEM),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    data = gene_of_interest,
    aes(x = genotype, y = count, color = treatment),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 2,
    shape = 4,
    stroke = 1
  ) +
  labs(title = gene_name, x = "Viral modification", y = "Normalised Counts") +
  theme_light(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 2),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_line(colour = "black"),
    legend.position = "right"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = bar_fill_colours, name = "Treatment") +
  scale_color_manual(values = bar_outline_colours, name = "Treatment") +
  remove_grid() + #Comment before "+" to see how ggplot scales Y axis.
  coord_cartesian(ylim = c(0, 16000)) # Adjust accordingly.


# Check how single gene plot looks before saving.
print(single_gene_plot)


# Save plot as TIFF.
save_tiff(
  paste0("results/figures/", design_formula, "/", gene_name, "_plot.tiff"),
  width = 80,
  height = 80,
  plot_code = function() {
    print(single_gene_plot)
  }
)


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/14_single_gene_plot.txt"
)
