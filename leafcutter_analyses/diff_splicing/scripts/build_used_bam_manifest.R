# GDC manifest for exactly the RNA-seq BAMs the leafcutter runs used, so
# they can be re-downloaded after the local BAMs are deleted.
#
# "Used" = every rail_id in any results/*/bams_used.tsv or
# control_results/*/bams_used.tsv. Those samples came from several
# download rounds (download/ plan of 2026-04-20 plus older per-gene
# rounds now in deprecated/, spop/, idh1/), so the file_id for each is
# looked up in every GDC metadata.cohort*.json saved next to a manifest.
# Each hit is checked against the BAM on disk by file size.
#
# 24 used samples had no saved JSON record or a size mismatch. For those,
# download/gdc_api_rnaseq_bams_2026-09-25.tsv holds a GDC API query (files
# endpoint, RNA-Seq BAMs) by BOTH the metadata barcode and the barcode in
# the BAM's own @RG SM tag (download/bam_sm_check_unresolved.tsv). 13 of
# them are a DIFFERENT ALIQUOT of the same patient than recount3 lists for
# that rail_id (e.g. metadata -11A normal, BAM -01A tumour): the manifest
# records the file that was actually on disk (`aliquot_on_disk`), since
# that is what the leafcutter results were computed from.
#
# Run from diff_splicing/:
#   /scratch/midway2/qhauck/conda_env/rstudio-server/bin/Rscript --vanilla scripts/build_used_bam_manifest.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(data.table)
})

OUT_MANIFEST <- "download/used_bams.gdc_manifest.txt"
OUT_MAP      <- "download/used_bams.rail_map.tsv"

used <- bind_rows(lapply(
  Sys.glob(c("results/*/bams_used.tsv", "control_results/*/bams_used.tsv")),
  function(f) read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    transmute(rail_id = as.integer(rail_id), run = basename(dirname(f)),
              set = basename(dirname(dirname(f)))))) %>%
  group_by(rail_id) %>%
  summarise(runs = paste(sort(unique(paste0(set, "/", run))), collapse = ","),
            .groups = "drop")
cat("used rail_ids:", nrow(used), "\n")

meta <- fread("../../../tcga_v2_samples.tsv", data.table = FALSE) %>%
  dplyr::select(rail_id, tcga_barcode)

jsons <- setdiff(list.files(".", pattern = "^metadata\\.cohort.*\\.json$",
                            recursive = TRUE, full.names = TRUE),
                 list.files(".snakemake", recursive = TRUE, full.names = TRUE))
files <- bind_rows(lapply(jsons, function(j) {
  d <- fromJSON(j, simplifyVector = FALSE)
  bind_rows(lapply(d, function(e) {
    if (!grepl("rna_seq\\.genomic\\.gdc_realn\\.bam$", e$file_name)) return(NULL)
    tibble(file_id = e$file_id, filename = e$file_name, md5 = e$md5sum,
           size = as.numeric(e$file_size), state = e$state,
           tcga_barcode = vapply(e$associated_entities,
                                 function(a) a$entity_submitter_id, ""),
           source_json = sub("^\\./", "", j))
  }))
}))
cat("RNA-seq BAM records across", length(jsons), "JSONs:", nrow(files), "\n")

sm  <- read_tsv("download/bam_sm_check_unresolved.tsv", show_col_types = FALSE)
api <- read_tsv("download/gdc_api_rnaseq_bams_2026-09-25.tsv", show_col_types = FALSE) %>%
  transmute(file_id, filename, md5, size = as.numeric(size), state,
            aliquot, source_json = "gdc_api_2026-09-25")

base <- used %>% left_join(meta, by = "rail_id") %>%
  left_join(sm %>% select(rail_id, bam_sm_barcode), by = "rail_id") %>%
  mutate(aliquot_on_disk = coalesce(bam_sm_barcode, tcga_barcode))
# Where the BAM's own SM tag is known to differ from the metadata barcode,
# use ONLY records for the aliquot actually on disk. Without this, once
# the BAMs are deleted (no size check to arbitrate) a JSON record for the
# metadata barcode could be picked for 55126 -- a file never used.
hits <- bind_rows(
  base %>% filter(is.na(bam_sm_barcode) | bam_sm_barcode == tcga_barcode) %>%
    left_join(files %>% rename(aliquot = tcga_barcode),
              by = c("tcga_barcode" = "aliquot"), relationship = "many-to-many"),
  base %>% inner_join(api, by = c("aliquot_on_disk" = "aliquot"),
                      relationship = "many-to-many"))

# On-disk BAM for the size check (inventory path relative to diff_splicing/).
inv <- read_tsv("download/bam_inventory.tsv", show_col_types = FALSE) %>%
  transmute(rail_id, home_bam, disk_size = file.size(home_bam))
# FROZEN after 2026-09-25: the used BAMs were deleted. An aliquot has
# several GDC RNA BAMs (genomic / transcriptome / chimeric) and only the
# byte-size match against the file on disk tells which one was used, so
# this script cannot be re-run without them. The saved, size-verified
# download/used_bams.* (read-only) are the record.
n_on_disk <- sum(file.exists(inv$home_bam[inv$rail_id %in% used$rail_id]))
if (n_on_disk < nrow(used))
  stop(sprintf(paste0("Only %d/%d used BAMs on disk -- refusing to rebuild. ",
                      "download/used_bams.gdc_manifest.txt is the frozen record."),
               n_on_disk, nrow(used)))
hits <- hits %>% left_join(inv, by = "rail_id") %>%
  mutate(size_match = !is.na(disk_size) & !is.na(size) & disk_size == size)

# One file per sample: prefer the size-verified record; several JSONs can
# list the same file_id, so collapse on it first.
res <- hits %>%
  group_by(rail_id) %>%
  arrange(desc(size_match), .by_group = TRUE) %>%
  filter(is.na(file_id) | file_id == first(file_id)) %>%
  summarise(across(c(tcga_barcode, aliquot_on_disk, runs, file_id, filename, md5, size, state,
                     home_bam, disk_size, size_match), first),
            source_json = paste(unique(source_json), collapse = ";"),
            .groups = "drop")

write_tsv(res, OUT_MAP)
res %>% filter(!is.na(file_id)) %>%
  transmute(id = file_id, filename, md5, size, state) %>%
  distinct() %>%   # 2 files back two rail_ids each (aliquot mix-up, see header)
  write_tsv(OUT_MANIFEST)

cat(sprintf("\nfound file_id: %d / %d\nsize-verified vs disk BAM: %d\nno file_id: %d\nsize mismatch or no disk BAM: %d\ndisk BAM is a different aliquot than metadata: %d\n",
            sum(!is.na(res$file_id)), nrow(res), sum(res$size_match),
            sum(is.na(res$file_id)), sum(!is.na(res$file_id) & !res$size_match),
            sum(res$aliquot_on_disk != res$tcga_barcode)))
cat("wrote", OUT_MANIFEST, "and", OUT_MAP, "\n")
