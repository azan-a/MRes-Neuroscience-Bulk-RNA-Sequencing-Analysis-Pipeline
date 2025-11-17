#Configuration------------------------------------------------------------------
library(sessioninfo)
library(pathview)
library(dplyr)


#Set parameters.
alpha <- 0.05
lfc_cutoff <- 0.26

design_formaula <- "singlefactor"
contrast_name <- "control_vs_drp1ca"

# Replace with KEGG pathway ID of interest
pathway_id <- "rno04930"


contrast_results <- readRDS(paste0("data/rds/", contrast_name, ".rds")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID") %>%
  dplyr::mutate(ncbi_gene_id = sub("GeneID:", "", GeneID)) %>%
  na.omit() %>%
  filter(padj <= alpha & abs(log2FoldChange) >= lfc_cutoff) %>%
  dplyr::select(ncbi_gene_id, log2FoldChange)


# Define inputs-----------------------------------------------------------------
gene_list <- contrast_results$log2FoldChange
names(gene_list) <- contrast_results$ncbi_gene_id

# Pathview----------------------------------------------------------------------
# Build pathview plot (automatically saves).
pathview(
  gene.data = gene_list,
  pathway.id = pathway_id,
  species = "rno",
  node.sum = "mean",
  low = list(gene = "#4A6FE3"),
  high = list(gene = "#D33F6A"),
  limit = list(gene = c(-1, 1)),
  bins = 20,
  kegg.dir = "data/pathview/",
  out.suffix = contrast_name,
)

# Move png file into figures directory.
png_file <- paste0(pathway_id, ".", contrast_name, ".png")
file.rename(
  png_file,
  file.path("results/figures/", design_formaula, "/", basename(png_file))
)


# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/17_pathview.txt"
)
