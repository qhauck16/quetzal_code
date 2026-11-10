library(dplyr)
library(readr)
library(jsonlite)
library(fs)

samples <- read_tsv('gdc_sample_sheet.2026-01-05.tsv')
full_manifest <- read_tsv('gdc_manifest.2026-01-05.173855.txt')

tumor_rnaseq <- samples %>%
  filter(`Data Type` == "Aligned Reads",
         `Data Category` == "Sequencing Reads",
         `Tumor Descriptor` == "Primary",
         grepl(".rna_seq.genomic.gdc_realn.bam", `File Name`)) %>%
  group_by(`Case ID`) %>%          # group by TCGA case barcode
  slice(1) %>%                     # take first file per case
  ungroup()


samples_full_info <- fromJSON("metadata.cohort.2026-01-06.json", flatten = TRUE)

temp <- samples_full_info %>%
  filter(data_format == 'BAM') %>%
  filter(file_id %in% tumor_rnaseq$`File ID`)

key <- bind_rows(temp$associated_entities)
key$file_id <- temp$file_id
  

full_info <- left_join(key, tumor_rnaseq, by = c('file_id' = 'File ID')) %>%
  left_join(select(tcga_metadata, c(rail_id, tcga_barcode, gdc_cases.project.project_id)), by = c('entity_submitter_id' = 'tcga_barcode'))



#actually moving the files 

base_dir <- "bams/"

full_info <- full_info %>%
  rename(
    project_id = gdc_cases.project.project_id
  )

for (i in seq_len(nrow(full_info))) {
  
  file_id   <- full_info$file_id[i]
  rail_id   <- full_info$rail_id[i]
  project   <- full_info$project_id[i]
  
  src_dir <- path(base_dir, file_id)
  dst_dir <- path(base_dir, project)
  
  # Create destination directory if needed
  dir_create(dst_dir)
  
  # Find BAM
  bam <- dir_ls(src_dir, regexp = "\\.bam$", type = "file")
  if (length(bam) != 1) {
    warning("Expected 1 BAM in ", src_dir)
    next
  }
  
  # Find BAI (.bai or .bam.bai)
  bai <- dir_ls(src_dir, regexp = "\\.bai$", type = "file")
  if (length(bai) != 1) {
    warning("Expected 1 BAI in ", src_dir)
    next
  }
  
  # Move + rename
  file_move(
    bam,
    path(dst_dir, paste0(rail_id, ".bam"))
  )

  file_move(
    bai,
    path(dst_dir, paste0(rail_id, ".bai"))
  )
  # 
  message(bam, " -> ", path(dst_dir, paste0(rail_id, ".bam")))
  message(bai, " -> ", path(dst_dir, paste0(rail_id, ".bai")))
}
