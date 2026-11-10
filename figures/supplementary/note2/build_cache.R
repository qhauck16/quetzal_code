# Cache for Supplementary Note 2, Figure 1 (was Supplementary Fig. S8): TCGA tumours carrying the NKAP
# neurodevelopmental-syndrome substitutions (p.R330C, p.R330H) vs their cohorts.
# Ported from tcga/genome_wide/flashier/diff_splicing/nkap_rare_disease.Rmd
# (same event set). Extended: PSI for EVERY sample of the two cohorts, not
# just the cohort median, so the figure can show the carrier's rank.
#
# Events: NKAP 20v60 LeafCutter, dPSI > 0.05, FC > 2.5, alt 3'ss whose cryptic
# acceptor is -8..+20 nt from the canonical 3'ss (df_final convention).
# Cohort membership comes from tcga_v2_samples (not the flashier sample set):
# the R330H tumour (857257, LAML) was removed by the pre-flashier QC, so it has
# no GSP #13 score but does have PSI.
# Any cwd. Reads one psi_table.tsv per event gene -> run on SLURM.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(tidyverse)
})
OUT <- file.path(QROOT, "figures/supplementary/note2/cache")

CARRIERS <- tibble(rail_id = c("886325", "857257"),
                   change  = c("R330C", "R330H"),
                   cohort  = c("TCGA-READ", "TCGA-LAML"))

df <- fread(file.path(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results"), "nkap/20v60/dpsi0.05_fc2.5/df_final.tsv"),
            data.table = FALSE)
events <- df %>% filter(type == "alt 3'ss", !is.na(distance_to_cryptic),
                        distance_to_cryptic >= -8, distance_to_cryptic <= 20)
cat(sprintf("[s19] df_final %d junctions; alt 3'ss in [-8, 20]: %d\n", nrow(df), nrow(events)))

meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
  transmute(rail_id = as.character(rail_id), cohort = gdc_cases.project.project_id) %>%
  filter(cohort %in% CARRIERS$cohort)

# Same column matching as the original: PSI headers are cluster:start-end;
# try start and start-1 and keep whichever matches more of the gene's events.
match_cols <- function(grp, psi_cols) {
  se <- sub("^[^:]+:", "", psi_cols)
  a <- paste0(grp$start, "-", grp$end); b <- paste0(grp$start - 1, "-", grp$end)
  if (sum(a %in% se) == 0 && sum(b %in% se) == 0) return(rep(NA_integer_, nrow(grp)))
  match(if (sum(a %in% se) >= sum(b %in% se)) a else b, se)
}

long <- events %>% split(list(.$chrom, .$gene_name), drop = TRUE) %>%
  lapply(function(grp) {
    f <- file.path(file.path(QROOT, "leafcutter_analyses/psi_files"), "tcga", grp$chrom[1], grp$gene_name[1], "psi_table.tsv")
    if (!file.exists(f)) return(NULL)
    psi <- fread(f, data.table = FALSE); colnames(psi)[1] <- "rail_id"
    psi$rail_id <- as.character(psi$rail_id)
    m <- match_cols(grp, setdiff(colnames(psi), "rail_id")); keep <- !is.na(m)
    if (!any(keep)) return(NULL)
    cols <- setdiff(colnames(psi), "rail_id")[m[keep]]
    psi <- psi[psi$rail_id %in% meta$rail_id, c("rail_id", cols), drop = FALSE]
    psi %>% pivot_longer(-rail_id, names_to = "event", values_to = "psi") %>%
      mutate(event = paste(grp$gene_name[1], event, sep = ":"))
  }) %>% bind_rows() %>%
  left_join(meta, by = "rail_id") %>% filter(is.finite(psi))

cat(sprintf("[s19] events with PSI: %d | samples: %s\n", n_distinct(long$event),
            paste(names(table(meta$cohort)), table(meta$cohort), collapse = ", ")))
for (i in seq_len(nrow(CARRIERS)))
  cat(sprintf("[s19] %s (%s) events with PSI: %d\n", CARRIERS$change[i], CARRIERS$rail_id[i],
              sum(long$rail_id == CARRIERS$rail_id[i])))
fwrite(long, file.path(OUT, "psi_long.tsv"), sep = "\t")
fwrite(CARRIERS, file.path(OUT, "carriers.tsv"), sep = "\t")
cat("[s19] wrote psi_long.tsv\n")
