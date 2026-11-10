library(dplyr)
library(readr)

samples <- read_tsv('gdc_sample_sheet.2026-01-05.tsv')
full_manifest <- read_tsv('gdc_manifest.2026-01-05.173855.txt')
# Filter for tumor RNA-seq BAMs
tumor_rnaseq <- samples %>%
  filter(`Data Type` == "Aligned Reads",
         `Data Category` == "Sequencing Reads",
         `Tumor Descriptor` == "Primary",
         grepl(".rna_seq.genomic.gdc_realn.bam", `File Name`)) %>%
  group_by(`Case ID`) %>%          # group by TCGA case barcode
  slice(1) %>%                     # take first file per case
  ungroup() %>%
  select(`File ID`, `File Name`)

new_manifest <- filter(full_manifest, filename %in% tumor_rnaseq$`File Name`)

# Write manifest
write_tsv(new_manifest, "manifest_tumor_rnaseq.tsv")
