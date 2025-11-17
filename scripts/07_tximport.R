library(sessioninfo)
library(tximport)

# Import sample metadata.
sample_table <- readxl::read_excel("config/sample_table.xlsx")

# Create named vector of quant files
quant_files <- setNames(sample_table$file, sample_table$sample_name)


# Import GTF file.
gtf_file <- rtracklayer::import(
  "data/raw/reference/GCF_036323735.1_GRCr8_genomic.gtf.gz"
)


# Extract transcript-to-gene mapping (using db_xref as gene ID).
tx2gene_df <- unique(
  data.frame(
    gtf_file[gtf_file$type == "exon"]
    )[, c("transcript_id", "db_xref")]
)


# Import Salmon quantifications.
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene_df)

# Save txi for DESeq2.
saveRDS(txi, "data/rds/txi.rds")

# Save session info
writeLines(
  capture.output(session_info()),
  "results/sessioninfo/07_tximport.txt"
)
