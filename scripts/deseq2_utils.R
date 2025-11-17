# Load required libraries
library(dplyr) # for %>%, mutate, select, arrange.
library(org.Rn.eg.db) # for mapIds.


# Transform contrast results from dds to a labelled data frame.
# Dependencies: dplyr, org.Rn.eg.db, tibble.
annotate_results <- function(contrast_results) {
  contrast_results_df <- contrast_results %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GeneID") %>%
    mutate(
      ncbi_gene_id = sub("GeneID:", "", GeneID),
      gene_symbol = mapIds(
        org.Rn.eg.db,
        keys = ncbi_gene_id,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
    ) %>%
    dplyr::select(ncbi_gene_id, gene_symbol, everything(), -GeneID) %>%
    arrange(desc(log2FoldChange))
}

#Save transformed data frame into an excel file with gene symbols.
#Dependencies: dplyr, writexl.
save_excel <- function(
  contrast_results,
  file_name,
  alpha = 0.05,
  lfc_cutoff = 0.26
) {
  results_df <- annotate_results(contrast_results)

  full_path_all <- file.path("results/tables/deg/", paste0(file_name, ".xlsx"))
  full_path_signif <- file.path(
    "results/tables/deg/",
    paste0(file_name, "_significant.xlsx")
  )

  writexl::write_xlsx(results_df, path = full_path_all)
  writexl::write_xlsx(
    results_df %>%
      filter(!is.na(padj)) %>%
      filter(padj <= alpha & abs(log2FoldChange) >= lfc_cutoff),
    path = full_path_signif
  )
}
