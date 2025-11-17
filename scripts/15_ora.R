#Configuration------------------------------------------------------------------
library(sessioninfo)
library(clusterProfiler)
library(org.Rn.eg.db)
library(dplyr)
library(writexl)

#Set parameters.
alpha <- 0.05
lfc_cutoff <- 0.26

contrast_name <- "drp1ca_vs_drp1ca_insulin"
design_formula <- "singlefactor"


background_genes <- readRDS(paste0(
  "data/rds/",
  design_formula,
  "_background_genes.rds"
))


# Extract significant DEGs from contrast results.
contrast_results <- readRDS(paste0("data/rds/", contrast_name, ".rds")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID") %>%
  dplyr::mutate(ncbi_gene_id = sub("GeneID:", "", GeneID)) %>%
  na.omit() %>%
  filter(padj <= alpha) %>%
  dplyr::select(ncbi_gene_id, log2FoldChange)

# Define inputs-----------------------------------------------------------------
upregulated_genes <- contrast_results %>%
  filter(log2FoldChange >= lfc_cutoff) %>%
  pull(ncbi_gene_id)

downregulated_genes <- contrast_results %>%
  filter(log2FoldChange <= -lfc_cutoff) %>%
  pull(ncbi_gene_id)

gene_sets <- list(
  upregulated = upregulated_genes,
  downregulated = downregulated_genes
)

#Lists the 3 aspects in GO.
onts <- c("BP", "MF", "CC")

# Stores results to be saved in Excel.
results_list_df <- list()


# ORA---------------------------------------------------------------------------
# KEGG
for (name in names(gene_sets)) {
  obj <- enrichKEGG(
    gene_sets[[name]],
    organism = "rno",
    universe = background_genes
  )
  key <- paste0(name, "_KEGG")
  results_list_df[[key]] <- as.data.frame(obj)
  saveRDS(obj, paste0("data/rds/", contrast_name, "_", key, ".rds"))
}

# GO
for (ont in onts) {
  for (name in names(gene_sets)) {
    obj <- enrichGO(
      gene_sets[[name]],
      OrgDb = "org.Rn.eg.db",
      ont = ont,
      universe = background_genes
    )
    key <- paste0(name, "_GO_", ont)
    results_list_df[[key]] <- as.data.frame(obj)
    saveRDS(obj, paste0("data/rds/", contrast_name, "_", key, ".rds"))
  }
}


# Save results to Excel with multiple sheets.
write_xlsx(
  results_list_df,
  paste0("results/tables/ora/", contrast_name, "_ora_results.xlsx")
)


# Save session info
writeLines(capture.output(session_info()), "results/sessioninfo/15_ora.txt")
