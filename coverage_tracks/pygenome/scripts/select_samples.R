
# Per-GSP sample selection straight from the current flashier RDS:
#   - condition  = top n_condition rail_ids by norm_L[, factor_key]
#   - controls   = n_control samples drawn from the rest of the cohort,
#                  per-project_id proportion-matched to the condition.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(flashier)
})

args <- commandArgs(trailingOnly = TRUE)
gsp           <- args[1]                 # e.g. "NKAP"
factor_key    <- args[2]                 # e.g. "k13"
rds_path      <- args[3]
tcga_metadata <- args[4]
n_condition   <- as.integer(args[5])
n_control     <- as.integer(args[6])
seed_val      <- as.integer(args[7])
out_path      <- args[8]                 # {gsp}/samples.tsv

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

cat("Loading RDS:", rds_path, "\n")
out <- readRDS(rds_path)

ldf <- flashier::ldf(out, type = "i")
L   <- ldf$L
colnames(L) <- paste0("k", seq_len(ncol(L)))
stopifnot(factor_key %in% colnames(L))

samples <- tibble(
  rail_id = as.integer(rownames(L)),
  loading = L[, factor_key]
) %>% arrange(desc(loading))

meta <- fread(tcga_metadata, data.table = FALSE) %>%
  dplyr::select(rail_id,
                tcga_barcode,
                project_id = gdc_cases.project.project_id) %>%
  mutate(project_short = str_split_i(project_id, "-", 2))

samples <- samples %>%
  inner_join(meta, by = "rail_id")

# --- condition ------------------------------------------------------
condition <- samples %>% slice_head(n = n_condition)

project_props <- condition %>%
  count(project_id, name = "count") %>%
  mutate(prop      = count / sum(count),
         n_control = round(prop * n_control))

cat("Condition (top", n_condition, "by", factor_key, "):\n")
print(condition %>% dplyr::select(rail_id, project_short, loading))
cat("\nControls per project_id (target", n_control, "):\n")
print(project_props)

# --- controls -------------------------------------------------------
control_pool <- samples %>%
  anti_join(condition, by = "rail_id")

set.seed(seed_val)
controls <- control_pool %>%
  inner_join(project_props %>% dplyr::select(project_id, n_control),
             by = "project_id") %>%
  group_by(project_id) %>%
  group_modify(~ {
    n <- dplyr::first(.x$n_control)
    dplyr::slice_sample(.x, n = min(n, nrow(.x)))
  }) %>%
  ungroup() %>%
  dplyr::select(-n_control)

# --- write samples.tsv ---------------------------------------------
samples_df <- bind_rows(
  condition %>% mutate(role = "condition"),
  controls  %>% mutate(role = "control")
) %>%
  distinct(rail_id, .keep_all = TRUE) %>%
  dplyr::select(rail_id, role, loading,
                tcga_barcode, project_id, project_short) %>%
  arrange(role, desc(loading))

fwrite(samples_df, out_path, sep = "\t")
cat("\nWrote", out_path, "(",
    sum(samples_df$role == "condition"), "condition /",
    sum(samples_df$role == "control"),   "controls)\n")
