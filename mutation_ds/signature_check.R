# Does the within-cohort mutant splicing signature carry over to mutants in
# OTHER cohorts? Signature = the comparison's LeafCutter junctions (cluster
# p.adjust < 0.05, |dPSI| > 0.05, mutant - control): score per tumour =
# mean PSI over the up junctions - mean PSI over the down junctions (PSI from
# the TCGA LeafCutter PSI tables, every flashier tumour).
#   Positive control: HELD-OUT LGG (R132 not among the 20 DS tests vs IDH1-
#     wild-type not among the 60 DS controls) -- does the signature generalise
#     within the cohort it was learnt in?
#   Test: each non-brain R132 tumour's percentile among the IDH1-wild-type
#     tumours of its OWN cohort; under no effect the percentiles are uniform
#     (mean 0.5): one-sample Wilcoxon of percentile - 0.5, one-sided > 0.
#   SPOP: the same with the F133 run's signature; the 3 hotspot carriers
#     outside PRAD are reported individually (too few to test).
# Out: out/signature_scores.tsv, out/signature_tests.tsv, out/signature_psi_<comp>.rds
# (junction PSI matrix of the home cohort, for the PCA). Any cwd -> SLURM.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(yaml); library(data.table); library(tidyverse); library(GenomicRanges) })
MD <- file.path(QROOT, "mutation_ds")
FIG <- file.path(QROOT, "figures")

# ---- samples ---------------------------------------------------------------
L <- fread(file.path(FIG, "shared_cache/gsp_loadings.tsv"), select = "rail_id")
meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "tcga_barcode", "gdc_cases.project.project_id"),
              data.table = FALSE) %>%
  transmute(rail_id = as.character(rail_id), cohort = gdc_cases.project.project_id,
            sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
              x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
  filter(rail_id %in% as.character(L$rail_id))
am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), select = c("sampleId", "hugoGeneSymbol", "proteinChange"), data.table = FALSE)
meta <- meta %>% filter(sampleId %in% am$sampleId) %>%
  mutate(idh1_any = sampleId %in% am$sampleId[am$hugoGeneSymbol == "IDH1"],
         r132     = sampleId %in% am$sampleId[am$hugoGeneSymbol == "IDH1" & grepl("R132", am$proteinChange)],
         spop_any = sampleId %in% am$sampleId[am$hugoGeneSymbol == "SPOP"],
         spop_hot = sampleId %in% am$sampleId[am$hugoGeneSymbol == "SPOP" &
                                                grepl("(F133|W131|F102|Y87)", am$proteinChange)])

# ---- signature junctions ---------------------------------------------------
load_sig <- function(comp) {
  base <- file.path(MD, "results", comp)
  ds <- fread(file.path(base, "leafcutter_ds_effect_sizes.txt"), data.table = FALSE)
  if (fread(file.path(base, "leaf_pheno_data.txt"), header = FALSE, nrows = 1)$V2 == "test")
    ds$deltapsi <- -ds$deltapsi
  sig <- fread(file.path(base, "leafcutter_ds_cluster_significance.txt"), data.table = FALSE) %>%
    mutate(cluster_id = sub(".*:", "", cluster))
  ds %>% separate(intron, c("chr", "s", "e", "cluster_id"), sep = ":", convert = TRUE) %>%
    mutate(strand = if_else(grepl("_\\+$", cluster_id), "+", "-"),
           start = as.integer(s) + 1L, end = as.integer(e) - 1L) %>%
    left_join(sig %>% select(cluster_id, p.adjust), by = "cluster_id") %>%
    filter(!is.na(p.adjust), p.adjust < 0.05, abs(deltapsi) > 0.05) %>%
    transmute(chr, strand, start, end, dir = ifelse(deltapsi > 0, "up", "down"), deltapsi)
}

# junction -> gene(s) whose span contains it on the same strand; PSI columns
# are cluster:start-end, matched on start or start-1 (build_s19_cache.R).
gr <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
genes <- gr[gr$type == "gene"]
psi_matrix <- function(j) {
  jg <- GRanges(j$chr, IRanges(j$start, j$end), strand = j$strand)
  ov <- findOverlaps(jg, genes, type = "within")
  cand <- tibble(i = queryHits(ov), gene = genes$gene_name[subjectHits(ov)], chr = j$chr[queryHits(ov)]) %>%
    distinct()
  out <- list(); found <- rep(FALSE, nrow(j))
  for (key in unique(paste(cand$chr, cand$gene))) {
    cg <- cand %>% filter(paste(chr, gene) == key, !found[i])
    if (!nrow(cg)) next
    f <- file.path(file.path(QROOT, "leafcutter_analyses/psi_files"), "tcga", cg$chr[1], cg$gene[1], "psi_table.tsv")
    if (!file.exists(f) || file.size(f) < 10) next
    hdr <- colnames(fread(f, nrows = 0)); se <- sub("^[^:]+:", "", hdr)
    a <- paste0(j$start[cg$i], "-", j$end[cg$i]); b <- paste0(j$start[cg$i] - 1, "-", j$end[cg$i])
    m <- ifelse(a %in% se, match(a, se), match(b, se))
    keep <- !is.na(m)
    if (!any(keep)) next
    x <- fread(f, select = c(1, m[keep]), data.table = FALSE)
    rows <- match(meta$rail_id, as.character(x[[1]]))
    mm <- as.matrix(x[rows, -1, drop = FALSE]); colnames(mm) <- as.character(cg$i[keep])
    out[[key]] <- mm; found[cg$i[keep]] <- TRUE
  }
  M <- do.call(cbind, out)
  list(M = M, idx = as.integer(colnames(M)), found = sum(found))
}

# keep: cohorts whose PSI matrix is saved for the PCA (note1/signature_plots.R)
score <- function(comp, keep) {
  j <- load_sig(comp); pm <- psi_matrix(j)
  dir <- j$dir[pm$idx]
  ks <- meta$cohort %in% keep
  M <- pm$M[ks, , drop = FALSE]; rownames(M) <- meta$rail_id[ks]
  saveRDS(list(M = M, junctions = j[pm$idx, ], cohort = meta$cohort[ks]),
          file.path(MD, "out", sprintf("signature_psi_%s.rds", comp)))
  cat(sprintf("[sig] %s: %d junctions (%d up, %d down); %d with PSI\n", comp, nrow(j),
              sum(j$dir == "up"), sum(j$dir == "down"), pm$found))
  up <- rowMeans(pm$M[, dir == "up", drop = FALSE], na.rm = TRUE)
  dn <- rowMeans(pm$M[, dir == "down", drop = FALSE], na.rm = TRUE)
  up - ifelse(is.finite(dn), dn, 0)
}

ds_ids <- function(comp) fread(file.path(MD, "samples", paste0(comp, ".tsv")), data.table = FALSE) %>%
  mutate(rail_id = as.character(rail_id))
meta$idh1_score <- score("IDH1_R132", c("TCGA-LGG", "TCGA-GBM"))
meta$spop_score <- score("SPOP_F133", "TCGA-PRAD")
fwrite(meta, file.path(MD, "out", "signature_scores.tsv"), sep = "\t")

pct_in_cohort <- function(d, score_col, mut_col, wt_filter) {
  d %>% group_by(cohort) %>% group_modify(function(x, k) {
    ref <- x[[score_col]][wt_filter(x)]
    x %>% filter(.data[[mut_col]]) %>%
      mutate(pct = vapply(.data[[score_col]], function(v) mean(ref < v, na.rm = TRUE), numeric(1)),
             n_ref = sum(is.finite(ref)))
  }) %>% ungroup()
}

idh_ds <- ds_ids("IDH1_R132")
lgg <- meta %>% filter(cohort == "TCGA-LGG", !rail_id %in% idh_ds$rail_id)
pc <- wilcox.test(lgg$idh1_score[lgg$r132], lgg$idh1_score[!lgg$idh1_any], alternative = "greater")
nb <- pct_in_cohort(meta %>% filter(!cohort %in% c("TCGA-LGG", "TCGA-GBM")), "idh1_score", "r132",
                    function(x) !x$idh1_any) %>% filter(is.finite(pct))
nbt <- wilcox.test(nb$pct - 0.5, alternative = "greater")

spop_ds <- ds_ids("SPOP_F133")
prad <- meta %>% filter(cohort == "TCGA-PRAD", !rail_id %in% spop_ds$rail_id)
pcs <- wilcox.test(prad$spop_score[prad$spop_any], prad$spop_score[!prad$spop_any], alternative = "greater")
ns <- pct_in_cohort(meta %>% filter(cohort != "TCGA-PRAD"), "spop_score", "spop_hot",
                    function(x) !x$spop_any)

tests <- tibble(
  test = c("IDH1 signature, held-out LGG: R132 vs IDH1-wt",
           "IDH1 signature, non-brain R132: own-cohort percentile > 0.5",
           "SPOP signature, held-out PRAD: SPOP-mut vs SPOP-wt"),
  n_mut = c(sum(lgg$r132), nrow(nb), sum(prad$spop_any)),
  n_ref = c(sum(!lgg$idh1_any), NA, sum(!prad$spop_any)),
  median_pct = c(NA, median(nb$pct), NA),
  p = c(pc$p.value, nbt$p.value, pcs$p.value))
print(as.data.frame(tests), digits = 3)
cat("\n[sig] non-brain R132 percentiles by cohort:\n")
print(as.data.frame(nb %>% group_by(cohort) %>% summarise(n = n(), median_pct = median(pct), .groups = "drop") %>%
                      arrange(desc(n))), digits = 3)
cat("\n[sig] SPOP hotspot carriers outside PRAD:\n")
print(as.data.frame(ns %>% select(cohort, rail_id, spop_score, pct, n_ref)), digits = 3)
fwrite(tests, file.path(MD, "out", "signature_tests.tsv"), sep = "\t")
fwrite(bind_rows(nb %>% mutate(gene = "IDH1"), ns %>% mutate(gene = "SPOP")) %>%
         select(gene, cohort, rail_id, pct, n_ref), file.path(MD, "out", "signature_percentiles.tsv"), sep = "\t")
