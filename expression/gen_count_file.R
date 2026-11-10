library(tidyverse)
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)

# --- Setup ------------------------------------------------------------------

gencode <- readRDS('../inputs/gencode_v46_granges.RDS')
tcga_metadata <- fread('../inputs/tcga_v2_samples.tsv', data.table = F)

# Ensembl ID -> gene_name mapping (version-stripped)
txdb          <- makeTxDbFromGRanges(gencode)
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
introns       <- unlist(introns_by_tx, use.names = FALSE)
tx2gene <- unique(mcols(gencode)[gencode$type == "transcript",
                                 c("transcript_id", "gene_id", "gene_name")])
introns$transcript_id <- names(introns_by_tx)[
  rep(seq_along(introns_by_tx), elementNROWS(introns_by_tx))
]
mcols(introns) <- cbind(mcols(introns),
                        tx2gene[match(introns$transcript_id, tx2gene$transcript_id), ])

ensembl_to_hgnc <- as.data.frame(introns) %>%
  dplyr::select(gene_id, gene_name) %>%
  mutate(gene_id_short = str_split_i(gene_id, "\\.", 1)) %>%
  dplyr::select(-gene_id) %>%
  unique()

# --- Per-cohort processing --------------------------------------------------

process_cohort <- function(file_path, cohort_name) {
  message(cohort_name)

  expr <- fread(file_path, data.table = FALSE)

  # Map external_id -> rail_id
  id_map    <- dplyr::select(tcga_metadata, external_id, rail_id) %>% distinct()
  new_names <- colnames(expr)
  new_names[-1] <- id_map$rail_id[match(new_names[-1], id_map$external_id)]
  colnames(expr) <- new_names

  # Keep all rail_id-mapped samples (drop only those whose external_id
  # has no rail_id mapping in tcga_metadata).
  keep_cols    <- !is.na(colnames(expr))
  keep_cols[1] <- TRUE                # gene_id column
  expr <- expr[, keep_cols, drop = FALSE]

  # Ensembl -> gene_name, deduplicate
  expr <- expr %>%
    mutate(gene_id_short = str_split_i(gene_id, "\\.", 1)) %>%
    inner_join(ensembl_to_hgnc, by = "gene_id_short") %>%
    dplyr::select(-gene_id, -gene_id_short) %>%
    group_by(gene_name) %>%
    slice_head(n = 1) %>%
    column_to_rownames("gene_name") %>%
    as.matrix()

  # Counts -> RPM
  rpm <- t(t(expr) / colSums(expr)) * 1e6

  write_tsv(
    rownames_to_column(as.data.frame(rpm), "gene_name"),
    paste0("tcga_counts/rpm_files/", cohort_name, "_rpm.tsv")
  )

  # Return both the unscaled RPM and the per-gene z-scored RPM so the
  # caller can build combined files for both.
  list(rpm    = rpm,
       scaled = t(scale(t(rpm))))
}

# --- Run --------------------------------------------------------------------

files   <- list.files("tcga_counts/",
                       pattern = "tcga\\.gene_sums\\.[A-Z]+\\.G[0-9]+\\.gz$",
                       full.names = TRUE)
cohorts <- str_extract(basename(files), "(?<=tcga\\.gene_sums\\.)[A-Z]+(?=\\.G[0-9]+\\.gz)")

all_results <- setNames(
  lapply(seq_along(files), function(i) process_cohort(files[i], cohorts[i])),
  cohorts
)

all_rpm    <- lapply(all_results, `[[`, "rpm")
all_scaled <- lapply(all_results, `[[`, "scaled")

# Align rows before cbind (common genes across all cohorts)
common_genes    <- Reduce(intersect, lapply(all_scaled, rownames))
combined_rpm    <- do.call(cbind, lapply(all_rpm,    function(m) m[common_genes, ]))
combined_scaled <- do.call(cbind, lapply(all_scaled, function(m) m[common_genes, ]))

write_tsv(
  rownames_to_column(as.data.frame(combined_rpm), "gene_name"),
  "combined_rpm.tsv"
)

write_tsv(
  rownames_to_column(as.data.frame(combined_scaled), "gene_name"),
  "combined_scaled_rpm.tsv"
)
