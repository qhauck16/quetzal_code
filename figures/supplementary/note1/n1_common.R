# Shared inputs for the Supplementary Note 1 figures (sourced after common.R).
#   n1_tumours(): flashier tumours WITH mutation calls (AlphaMissense-merged
#     table), with cohort, GSP #1 / #12 log-loadings and IDH1 / SPOP status.
#   n1_sf3b1(), n1_nkap(): the Fig. 2d / 3a tables (variant tables) with cohort.
#   fmt_p(), fmt_fold(): label formatting.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source(file.path(FIG_ROOT, "scripts", "vt_paths.R"))
source(file.path(FIG_ROOT, "figure3", "scripts", "nkap_variants.R"))   # classify_variants
N1_MDS <- file.path(QROOT, "mutation_ds")
N1_OUT <- file.path(OUT_DIR, "note1"); dir.create(N1_OUT, showWarnings = FALSE)
N1_HOT <- "(F133|W131|F102|Y87)"

n1_cohorts <- function() fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "gdc_cases.project.project_id"),
                               data.table = FALSE) %>%
  transmute(rail_id = as.character(rail_id), cohort = gdc_cases.project.project_id)

n1_tumours <- function() {
  L <- fread(file.path(FIG_ROOT, "shared_cache", "gsp_loadings.tsv"), select = c("rail_id", "k1", "k12"),
             data.table = FALSE) %>% mutate(rail_id = as.character(rail_id))
  am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), select = c("sampleId", "hugoGeneSymbol", "proteinChange"),
              data.table = FALSE)
  g <- function(gene, re = ".") unique(am$sampleId[am$hugoGeneSymbol == gene & grepl(re, am$proteinChange)])
  hot <- am %>% filter(hugoGeneSymbol == "SPOP", grepl(N1_HOT, proteinChange)) %>%
    group_by(sampleId) %>% summarise(lab = paste(unique(proteinChange), collapse = ", "), .groups = "drop")
  hot_lab <- setNames(hot$lab, hot$sampleId)
  fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "tcga_barcode", "gdc_cases.project.project_id"),
        data.table = FALSE) %>%
    transmute(rail_id = as.character(rail_id), cohort = gdc_cases.project.project_id,
              sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
                x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
    inner_join(L, by = "rail_id") %>%
    filter(sampleId %in% am$sampleId) %>%
    mutate(logk1 = log(k1), logk12 = log(k12),
           idh1_any = sampleId %in% g("IDH1"), r132 = sampleId %in% g("IDH1", "R132"),
           spop_any = sampleId %in% g("SPOP"), spop_hot = sampleId %in% g("SPOP", N1_HOT),
           spop_change = unname(hot_lab[sampleId]))
}
n1_sf3b1 <- function() fread(vt_path("SF3B1_full_info"), data.table = FALSE) %>%
  mutate(rail_id = as.character(rail_id)) %>% left_join(n1_cohorts(), by = "rail_id")
n1_nkap <- function() classify_variants(fread(vt_path("NKAP_full_plus_variants"), data.table = FALSE)) %>%
  mutate(rail_id = as.character(rail_id)) %>% left_join(n1_cohorts(), by = "rail_id")

fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
fmt_fold <- function(f) ifelse(f >= 100, paste0(formatC(f, format = "d", big.mark = ","), "-fold"),
                               sprintf("%.1f-fold", f))
