#Configuration------------------------------------------------------------------
library(sessioninfo)
library(clusterProfiler)
library(org.Rn.eg.db)
library(dplyr)
library(writexl)


contrast_name <- "drp1ca_vs_drp1ca_insulin"

#Extract DEGs from contrast results.
contrast_results <- readRDS(paste0("data/rds/", contrast_name, ".rds")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID") %>%
  dplyr::mutate(ncbi_gene_id = sub("GeneID:", "", GeneID)) %>%
  na.omit() %>%
  dplyr::select(ncbi_gene_id, log2FoldChange) %>%
  arrange(desc(log2FoldChange))


# Define inputs-----------------------------------------------------------------
gene_list <- contrast_results$log2FoldChange
names(gene_list) <- contrast_results$ncbi_gene_id

#Lists the 3 aspects in GO.
onts <- c("BP", "MF", "CC")

# Stores results to be saved in Excel.
results_list_df <- list()


# GSEA--------------------------------------------------------------------------
# KEGG
obj <- gseKEGG(geneList = gene_list, organism = "rno")
results_list_df[["KEGG"]] <- as.data.frame(obj)
saveRDS(obj, paste0("data/rds/", contrast_name, "_GSEA_KEGG.rds"))

# GO
for (ont in onts) {
  obj <- gseGO(geneList = gene_list, OrgDb = org.Rn.eg.db, ont = ont)
  key <- paste0("GO_", ont)
  results_list_df[[key]] <- as.data.frame(obj)
  saveRDS(obj, paste0("data/rds/", contrast_name, "_GSEA_", key, ".rds"))
}


# Save results to Excel with multiple sheets.
write_xlsx(
  results_list_df,
  paste0("results/tables/gsea/", contrast_name, "_gsea_results.xlsx")
)


# Save session info
writeLines(capture.output(session_info()), "results/sessioninfo/16_gsea.txt")
