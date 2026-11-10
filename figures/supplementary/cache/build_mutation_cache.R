# Cache for supplementary figures S10 and S12: the AlphaMissense (>0.9 or
# NA) mutation matrix and cohort-stratified, one-sided CMH tests, written
# as the small tables the figures need. Loadings come from the
# *_full_plus_variants.tsv tables (k10 / k13 columns = norm_L), so the
# 3.7 GB flashier RDS is never loaded.
#   sbatch --mem=24G --wrap "Rscript --vanilla cache/build_mutation_cache.R"
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
source(file.path(QROOT, "figures/scripts/vt_paths.R"))  # vt_path()
OUT <- file.path(QROOT, "figures/supplementary/cache")

tcga_metadata <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
  mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
    x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }))
am_mut_data <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"),
                     data.table = FALSE)
mut_data_am09 <- am_mut_data %>% filter(is.na(am_pathogenicity) | am_pathogenicity > 0.9)

m <- mut_data_am09 %>% dplyr::select(hugoGeneSymbol, sampleId) %>% mutate(value = 1) %>%
  unique() %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = value, values_fill = 0) %>%
  right_join(dplyr::select(tcga_metadata, c(rail_id, sampleId, gdc_cases.project.project_id)),
             by = "sampleId") %>%
  dplyr::select(-sampleId) %>% column_to_rownames("rail_id") %>% na.omit()
meta_col <- "gdc_cases.project.project_id"
gcols <- setdiff(colnames(m), meta_col)
m <- m[, c(gcols[colSums(m[, gcols, drop = FALSE]) > 0], meta_col), drop = FALSE]
mut_mat <- m
cat("mut_mat:", dim(mut_mat), "\n")

loadings <- function(file, col) {
  d <- fread(vt_path(file), data.table = FALSE)
  setNames(d[[col]], as.character(d$rail_id))
}
L10 <- loadings("SF3B1_full_plus_variants", "k10")
L13 <- loadings("NKAP_full_plus_variants",  "k13")

cmh <- function(top) {
  is_top <- as.integer(rownames(mut_mat) %in% top)
  gene_cols  <- setdiff(colnames(mut_mat), meta_col)
  mut_to_use <- as.matrix(mut_mat[, gene_cols])
  mut_to_use <- mut_to_use[, colSums(mut_to_use[is_top == 1, , drop = FALSE]) > 2, drop = FALSE]
  keep <- colnames(mut_to_use); n <- length(keep)
  strata <- split(seq_along(is_top), mut_mat[, meta_col])
  strata <- strata[sapply(strata, length) >= 2]
  sa <- sE <- sV <- sad <- sbc <- numeric(n); nv <- integer(n)
  for (idx in strata) {
    ns <- length(idx); n1 <- sum(is_top[idx]); n0 <- ns - n1
    if (n1 == 0 || n0 == 0) next
    m1 <- colSums(mut_to_use[idx, , drop = FALSE]); m0 <- ns - m1
    ok <- m1 > 0 & m0 > 0
    a <- as.numeric(crossprod(is_top[idx], mut_to_use[idx, , drop = FALSE]))
    b <- n1 - a; c <- m1 - a; d <- n0 - c
    sa <- sa + a; sE <- sE + ifelse(ok, n1 * m1 / ns, 0)
    sV <- sV + ifelse(ok, n1 * n0 * m1 * m0 / (ns^2 * (ns - 1)), 0)
    sad <- sad + ifelse(ok, a * d / ns, 0); sbc <- sbc + ifelse(ok, b * c / ns, 0)
    nv <- nv + as.integer(ok)
  }
  z <- ifelse(nv >= 1 & sV > 0, (sa - sE) / sqrt(sV), NA_real_)
  r <- data.frame(gene = keep, n_top_mutated = as.integer(sa),
                  odds_ratio = round(sad / pmax(sbc, 1e-10), 2),
                  p_value = pnorm(z, lower.tail = FALSE))
  r$padj <- pmin(r$p_value * length(gene_cols), 1)
  list(res = r[order(r$p_value), ], mut_to_use = mut_to_use)
}

## ---- S9: SUGP1 along GSP #10 (top 100, SF3B1-mutated excluded) -------
top10 <- names(sort(L10, decreasing = TRUE))[1:100]
t10 <- top10[top10 %in% rownames(mut_mat)]
t10 <- t10[mut_mat[t10, "SF3B1"] == 0]
r10 <- cmh(t10)$res
fwrite(r10, file.path(OUT, "s9_cmh_sf3b1_excluded.tsv"), sep = "\t")

carriers <- function(g) mut_data_am09 %>% filter(hugoGeneSymbol == g) %>%
  inner_join(dplyr::select(tcga_metadata, rail_id, sampleId), by = "sampleId") %>%
  mutate(rail_id = as.character(rail_id)) %>% group_by(rail_id) %>%
  summarise(label = paste(unique(proteinChange), collapse = ","), .groups = "drop")
hot <- c('R625H','R625C','K700E','K666T','G740E','K666E','K666N','K741E','K741N','K741Q','G740V')
cls <- function(l) { if (is.na(l)) return(NA_character_)
  h <- strsplit(l, ",")[[1]]; p <- suppressWarnings(as.numeric(regmatches(h, regexpr("\\d+", h))))
  if (any(h %in% hot)) return("SF3B1 hotspot mutation")
  if (any(p > 620 & p < 800, na.rm = TRUE)) return("Mutation Near SF3B1 Hotspot")
  "Other SF3B1 Mutation" }
s9 <- tibble(rail_id = names(L10), k10 = unname(L10)) %>% arrange(desc(k10)) %>%
  slice_head(n = 100) %>% mutate(rank = row_number()) %>%
  left_join(carriers("SF3B1") %>% rename(sf3b1_label = label), by = "rail_id") %>%
  left_join(carriers("SUGP1") %>% rename(sugp1_label = label), by = "rail_id") %>%
  mutate(sf3b1_class = vapply(sf3b1_label, cls, character(1)),
         category = case_when(!is.na(sf3b1_class) ~ sf3b1_class,
                              !is.na(sugp1_label) ~ "SUGP1 carrier",
                              TRUE ~ "No SF3B1/SUGP1 Mutation"))
fwrite(s9, file.path(OUT, "s9_top100_gsp10.tsv"), sep = "\t")

## ---- S12: co-occurrence in top 50 on GSP #13 without NKAP mutation ----
top13 <- names(sort(L13, decreasing = TRUE))[1:50]
t13 <- top13[top13 %in% rownames(mut_mat)]
t13 <- t13[mut_mat[t13, "NKAP"] == 0]
cr <- cmh(t13); r13 <- cr$res
fwrite(r13, file.path(OUT, "s12_cmh_nkap_excluded.tsv"), sep = "\t")
sig <- r13$gene[r13$padj < 0.05]
sig <- sig[order(r13$padj[match(sig, r13$gene)])]
M <- vapply(sig, function(g) cr$mut_to_use[t13, g], numeric(length(t13)))
co <- crossprod(M)
fwrite(data.frame(row = rep(rownames(co), ncol(co)), col = rep(colnames(co), each = nrow(co)),
                  count = as.vector(co)), file.path(OUT, "s12_cooccurrence.tsv"), sep = "\t")
cat("S9 SUGP1:", format(r10$padj[r10$gene == "SUGP1"]), " S12 sig:", paste(sig, collapse = ","), "\n")
