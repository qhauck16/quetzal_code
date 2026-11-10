# Within-cohort NULL comparisons for the IDH1 / SPOP note: wild-type vs
# wild-type LeafCutter runs in TCGA-LGG (IDH1 wild-type) and TCGA-PRAD (SPOP
# wild-type), 20 v 60, so the mutant runs' differential-splicing counts can be
# compared with what the same cohort gives by chance. Same universe, aliquot
# dedupe and matching as select_samples.R: a random (seeded) 20 "test" from the
# wild-type pool, then MatchIt 1:3 nearest on mean read length + FFPE from the
# remaining wild-type tumours. Roles are labelled test/control so the same
# Snakefile runs them; neither group carries a mutation in the gene.
# Writes samples/NULL_<cohort>.tsv and samples/download_list_null.tsv.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(yaml); library(data.table); library(tidyverse); library(MatchIt) })
ROOT <- QROOT
OUT  <- file.path(ROOT, "mutation_ds", "samples")
SEED <- 2026; RATIO <- 3; N_TEST <- 20

meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
  transmute(rail_id, tcga_barcode, project = gdc_cases.project.project_id, avg_len,
            ffpe = ifelse(is.na(gdc_cases.samples.is_ffpe), FALSE, gdc_cases.samples.is_ffpe)) %>%
  mutate(sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
    x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
  filter(!is.na(avg_len))
in_flash <- as.integer(fread(file.path(ROOT, "figures/figure2/cache/SF3B1_full_info.tsv"),
                             select = "rail_id")[[1]])
meta <- meta %>% mutate(vial = substr(tcga_barcode, 16, 16)) %>%
  arrange(sampleId, vial != "A", !rail_id %in% in_flash, rail_id) %>%
  distinct(sampleId, .keep_all = TRUE)
am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), data.table = FALSE)
called <- unique(am$sampleId)

NULLS <- tribble(~comparison, ~gene, ~cohort, ~gsp,
  "NULL_LGG",  "IDH1", "TCGA-LGG",  1L,
  "NULL_PRAD", "SPOP", "TCGA-PRAD", 12L)
set.seed(SEED)
sheets <- lapply(seq_len(nrow(NULLS)), function(i) {
  g  <- NULLS[i, ]
  wt <- meta %>% filter(project == g$cohort, sampleId %in% called,
                        !sampleId %in% am$sampleId[am$hugoGeneSymbol == g$gene])
  fake <- wt %>% slice_sample(n = N_TEST)
  d <- bind_rows(fake %>% mutate(group = "test"),
                 wt %>% filter(!rail_id %in% fake$rail_id) %>% mutate(group = "control")) %>%
    mutate(group = factor(group, levels = c("control", "test")))
  md <- match.data(matchit(group ~ avg_len + ffpe, data = d, method = "nearest", ratio = RATIO)) %>%
    as_tibble() %>%
    transmute(comparison = g$comparison, gene = g$gene, gsp = g$gsp, role = as.character(group),
              rail_id, tcga_barcode, project, avg_len, ffpe)
  cat(sprintf("[null] %s: %d vs %d (pool %d %s-wild-type)\n", g$comparison,
              sum(md$role == "test"), sum(md$role == "control"), nrow(wt), g$gene))
  fwrite(md, file.path(OUT, paste0(g$comparison, ".tsv")), sep = "\t")
  md
})
all_s <- bind_rows(sheets)
stopifnot(!any(duplicated(all_s[, c("comparison", "rail_id")])))
dl <- all_s %>% distinct(rail_id, tcga_barcode, project)
have <- function(r) file.exists(file.path(dirname(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results")), "junc_files", paste0(r, ".junc"))) |
                    file.exists(file.path(ROOT, "mutation_ds", "junc_files", paste0(r, ".junc")))
dl$have_junc <- have(dl$rail_id)
fwrite(dl, file.path(OUT, "download_list_null.tsv"), sep = "\t")
cat(sprintf("[null] unique samples: %d; .junc on disk: %d; to download: %d\n",
            nrow(dl), sum(dl$have_junc), sum(!dl$have_junc)))
