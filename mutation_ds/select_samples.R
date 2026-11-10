# Supplementary Note 1 (IDH1 / SPOP): mutation-defined differential-splicing
# comparisons within a single TCGA cohort, to test whether each gene's own
# splicing effect matches the GSP it associates with (IDH1 -> GSP #1 in LGG,
# SPOP -> GSP #12 in PRAD).
#
# One comparison per mutation group, mutants vs 3x matched wild-type controls
# (same 1:3 design and matching as the GSP comparisons in the main text):
#   IDH1 R132 (LGG)  : 20 mutants drawn at random (seeded) from all carriers
#   SPOP F133 / W131 / F102 / Y87 (PRAD): ALL carriers of that residue
# Mutation calls: the AlphaMissense-merged PanCanAtlas + GDC table (any call at
# the residue). Controls: same cohort, have mutation calls, NO call in the gene;
# MatchIt 1:3 nearest neighbour on mean read length + FFPE (exact on cohort).
# Writes mutation_ds/samples/<comparison>.tsv and samples/download_list.tsv.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(tidyverse); library(MatchIt)
})
ROOT <- QROOT
OUT  <- file.path(ROOT, "mutation_ds", "samples"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED <- 123; RATIO <- 3; N_IDH1 <- 20

meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
  transmute(rail_id, tcga_barcode, project = gdc_cases.project.project_id, avg_len,
            ffpe = ifelse(is.na(gdc_cases.samples.is_ffpe), FALSE, gdc_cases.samples.is_ffpe)) %>%
  mutate(sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
    x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
  filter(!is.na(avg_len))
# One RNA aliquot per tumour sample: TCGA has samples with both 01A and 01B
# vials, and some identical barcodes under two rail_ids. Prefer vial A, then a
# rail_id in the flashier fit, then the lowest rail_id.
in_flash <- as.integer(fread(file.path(ROOT, "figures/figure2/cache/SF3B1_full_info.tsv"),
                             select = "rail_id")[[1]])
meta <- meta %>%
  mutate(vial = substr(tcga_barcode, 16, 16), patient = substr(tcga_barcode, 1, 12)) %>%
  arrange(sampleId, vial != "A", !rail_id %in% in_flash, rail_id) %>%
  distinct(sampleId, .keep_all = TRUE)
# Different tumour samples of one patient (e.g. primary 01 + metastasis 06) are
# kept as separate samples, as in the GSP comparisons of the main text.
am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), data.table = FALSE) %>%
  mutate(res = as.integer(str_extract(proteinChange, "\\d+")))
called <- unique(am$sampleId)

GROUPS <- tribble(
  ~comparison,  ~gene,  ~cohort,     ~residue, ~gsp,
  "IDH1_R132",  "IDH1", "TCGA-LGG",  132L,     1L,
  "SPOP_F133",  "SPOP", "TCGA-PRAD", 133L,     12L,
  "SPOP_W131",  "SPOP", "TCGA-PRAD", 131L,     12L,
  "SPOP_F102",  "SPOP", "TCGA-PRAD", 102L,     12L,
  "SPOP_Y87",   "SPOP", "TCGA-PRAD", 87L,      12L)

set.seed(SEED)
sheets <- lapply(seq_len(nrow(GROUPS)), function(i) {
  g <- GROUPS[i, ]
  coh <- meta %>% filter(project == g$cohort, sampleId %in% called)
  gene_mut <- unique(am$sampleId[am$hugoGeneSymbol == g$gene])
  carriers <- unique(am$sampleId[am$hugoGeneSymbol == g$gene & am$res %in% g$residue])
  mut <- coh %>% filter(sampleId %in% carriers)
  if (g$gene == "IDH1") mut <- mut %>% slice_sample(n = min(N_IDH1, nrow(mut)))
  wt  <- coh %>% filter(!sampleId %in% gene_mut)
  d <- bind_rows(mut %>% mutate(group = "test"), wt %>% mutate(group = "control")) %>%
    mutate(group = factor(group, levels = c("control", "test")))
  m <- matchit(group ~ avg_len + ffpe, data = d, method = "nearest", ratio = RATIO)
  md <- match.data(m) %>% as_tibble() %>%
    transmute(comparison = g$comparison, gene = g$gene, gsp = g$gsp, role = as.character(group),
              rail_id, tcga_barcode, project, avg_len, ffpe)
  cat(sprintf("[select] %-10s %s: %d mutants (of %d carriers) vs %d controls (pool %d wild-type)\n",
              g$comparison, g$cohort, sum(md$role == "test"), length(intersect(carriers, coh$sampleId)),
              sum(md$role == "control"), nrow(wt)))
  fwrite(md, file.path(OUT, paste0(g$comparison, ".tsv")), sep = "\t")
  md
})
all_s <- bind_rows(sheets)

# A sample may be a control in one SPOP comparison and never a case in another
# (residue groups are disjoint; controls carry no SPOP call at all).
stopifnot(!any(duplicated(all_s[, c("comparison", "rail_id")])),
          !any(duplicated(paste(all_s$comparison, all_s$tcga_barcode))))
dl <- all_s %>% distinct(rail_id, tcga_barcode, project)
# All TCGA BAMs were deleted on 2026-09-25, but the diff-splicing .junc files
# (same regtools settings, relabelled by the aliquot fix) were kept: a sample
# with a .junc there needs no download.
JUNC <- file.path(dirname(file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results")), "junc_files")
dl$have_junc <- file.exists(file.path(JUNC, paste0(dl$rail_id, ".junc")))
fwrite(dl, file.path(OUT, "download_list.tsv"), sep = "\t")
cat(sprintf("[select] unique samples: %d (LGG %d, PRAD %d); .junc on disk: %d; to download: %d\n",
            nrow(dl), sum(dl$project == "TCGA-LGG"), sum(dl$project == "TCGA-PRAD"),
            sum(dl$have_junc), sum(!dl$have_junc)))
