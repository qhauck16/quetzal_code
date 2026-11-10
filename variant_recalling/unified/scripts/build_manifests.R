
suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(data.table)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
gene_dir          <- args[1]   # e.g. nkap (relative to unified/)
existing_rna_bams <- args[2]   # e.g. ../../diff_splicing/nkap_bams ; "NA" if none

if (existing_rna_bams == "NA") existing_rna_bams <- NA_character_

samples_path  <- file.path(gene_dir, "samples.tsv")
manifest_path <- file.path(gene_dir, "gdc_manifest.txt")
metadata_path <- file.path(gene_dir, "metadata.cohort.json")

stopifnot(file.exists(samples_path))
stopifnot(file.exists(manifest_path))
stopifnot(file.exists(metadata_path))

samples <- fread(samples_path, data.table = FALSE) %>% as_tibble()

manifest <- fread(manifest_path, data.table = FALSE) %>% as_tibble()
metadata <- fromJSON(metadata_path, flatten = TRUE) %>% as_tibble()

metadata <- metadata %>%
  filter(data_format == "BAM",
         file_id %in% manifest$id)

# associated_entities: assume one aliquot per BAM (matches prior scripts)
key <- bind_rows(metadata$associated_entities)
key$file_id <- metadata$file_id
key$experimental_strategy <- metadata$experimental_strategy
if ("analysis.workflow_type" %in% colnames(metadata)) {
  key$workflow <- metadata$analysis.workflow_type
} else {
  key$workflow <- NA_character_
}

key <- key %>%
  mutate(
    patient_id  = substr(entity_submitter_id, 1, 12),
    sample_code = as.numeric(substr(entity_submitter_id, 14, 15)),
    sample_type = case_when(
      sample_code >= 1 & sample_code <= 9   ~ "Tumor",
      sample_code >= 10                     ~ "Normal",
      TRUE                                  ~ NA_character_
    )
  )

# Each rail_id is ONE aliquot. Until 2026-09-25 files were joined to
# rail_ids by 12-char patient ID, so a patient's other sample (another
# vial, a metastasis/recurrence) could stand in for the rail's own RNA:
# 2 SF3B1 + 2 NKAP rna_manifest rows were a different aliquot. Now:
#   RNA-Seq    -> the rail's exact aliquot barcode (tcga_barcode)
#   WXS tumour -> same SAMPLE as the rail (first 16 chars, e.g.
#                 TCGA-44-2665-01A; DNA aliquots differ from RNA ones)
#   WXS normal -> any normal of the same patient
rail_key <- samples %>%
  transmute(rail_id, role, tcga_barcode,
            rail_sample  = substr(tcga_barcode, 1, 16),
            rail_patient = substr(tcga_barcode, 1, 12))

key <- key %>% mutate(file_sample = substr(entity_submitter_id, 1, 16))

# ------- WXS: paired Tumor + Normal -------
wxs_t <- key %>%
  filter(experimental_strategy == "WXS", sample_type == "Tumor") %>%
  inner_join(rail_key, by = c("file_sample" = "rail_sample"))
wxs_n <- key %>%
  filter(experimental_strategy == "WXS", sample_type == "Normal") %>%
  inner_join(rail_key %>% dplyr::select(-rail_sample),
             by = c("patient_id" = "rail_patient"))
wxs <- bind_rows(wxs_t, wxs_n) %>%
  group_by(rail_id) %>%
  filter(all(c("Tumor", "Normal") %in% sample_type)) %>%
  group_by(rail_id, sample_type) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(rail_id, role, patient_id, sample_type,
                file_id, entity_submitter_id)

# ------- RNA-Seq: primary tumor STAR 2-Pass Genome -------
rna <- key %>%
  filter(experimental_strategy == "RNA-Seq",
         (is.na(workflow) | workflow == "STAR 2-Pass Genome")) %>%
  inner_join(rail_key %>% dplyr::select(rail_id, role, tcga_barcode),
             by = c("entity_submitter_id" = "tcga_barcode")) %>%
  group_by(rail_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(rail_id, role, patient_id,
                file_id, entity_submitter_id) %>%
  mutate(sample_type = "Tumor")

fwrite(wxs, file.path(gene_dir, "wxs_manifest.tsv"), sep = "\t")
fwrite(rna, file.path(gene_dir, "rna_manifest.tsv"), sep = "\t")

cat("WXS paired patients:", length(unique(wxs$patient_id)),
    " (", nrow(wxs), "BAMs)\n")
cat("RNA-Seq tumor BAMs available in metadata:", nrow(rna), "\n")

# ------- Stage existing RNA-seq BAMs (symlink) -------
bams_dir <- file.path(gene_dir, "bams")
dir.create(bams_dir, showWarnings = FALSE, recursive = TRUE)

candidate_dirs <- c(
  existing_rna_bams,
  file.path("..", basename(gene_dir), "rnaseq_calling/bams")  # legacy
)
candidate_dirs <- candidate_dirs[!is.na(candidate_dirs) & dir.exists(candidate_dirs)]

stage_rna <- function(rail_id) {
  target_bam <- file.path(bams_dir, paste0(rail_id, "_Tumor.rna.bam"))
  if (file.exists(target_bam) || file.exists(paste0(target_bam, ".bai"))) {
    return(target_bam)
  }
  for (d in candidate_dirs) {
    for (n in c(paste0(rail_id, "_Tumor.rna.bam"),
                paste0(rail_id, ".bam"))) {
      src <- file.path(d, n)
      if (file.exists(src)) {
        file.symlink(normalizePath(src), target_bam)
        # index: try adjacent .bai or .bam.bai
        for (idx_ext in c(".bai", ".bam.bai")) {
          src_idx <- if (idx_ext == ".bai") sub("\\.bam$", ".bai", src)
                     else                    paste0(src, ".bai")
          if (file.exists(src_idx)) {
            tgt_idx <- sub("\\.bam$", ".bai", target_bam)
            if (!file.exists(tgt_idx))
              file.symlink(normalizePath(src_idx), tgt_idx)
            break
          }
        }
        return(src)
      }
    }
  }
  return(NA_character_)
}

rna_staged <- vapply(rna$rail_id, stage_rna, character(1))
n_rna_symlinked <- sum(!is.na(rna_staged))
cat("RNA-Seq BAMs symlinked from:", paste(candidate_dirs, collapse = ", "),
    "->", n_rna_symlinked, "/", nrow(rna), "\n")

# Build to_download.tsv (only files not already staged)
to_dl_rna <- rna %>%
  mutate(symlinked_from = unname(rna_staged)) %>%
  filter(is.na(symlinked_from)) %>%
  mutate(
    target_name = paste0(rail_id, "_Tumor.rna.bam"),
    source      = "rna"
  ) %>%
  dplyr::select(file_id, source, target_name, rail_id, role)

to_dl_wxs <- wxs %>%
  mutate(
    target_name = paste0(rail_id, "_", sample_type, ".wxs.bam"),
    source      = "wxs"
  ) %>%
  dplyr::select(file_id, source, target_name, rail_id, role)

to_dl <- bind_rows(to_dl_wxs, to_dl_rna)
fwrite(to_dl, file.path(gene_dir, "to_download.tsv"), sep = "\t")
cat("Files left to download:", nrow(to_dl),
    "(WXS:", nrow(to_dl_wxs), ", RNA:", nrow(to_dl_rna), ")\n")

# Convenience: a sub-manifest for `gdc-client download -m`
gdc_manifest_full <- fread(manifest_path, data.table = FALSE) %>% as_tibble()
to_dl_manifest <- gdc_manifest_full %>%
  filter(id %in% to_dl$file_id)
fwrite(to_dl_manifest,
       file.path(gene_dir, "to_download.manifest.tsv"), sep = "\t")
cat("Wrote to_download.manifest.tsv (",
    nrow(to_dl_manifest), "files)\n")

writeLines("done", file.path(gene_dir, ".manifests_built"))
