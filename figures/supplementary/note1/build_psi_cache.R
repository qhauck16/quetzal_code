.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
# Cache for the Note 1 PSI CDFs: per-junction mean PSI of the GSP #1 and
# GSP #12 driving junctions (mutation_ds/gsp_top/fact_{1,12}_top_junctions.tsv,
# all splice types) in four sample groups each, to separate a tissue
# (covariate) effect from a mutation (splicing) effect:
#   GSP #1 : brain (LGG + GBM) vs non-brain;  IDH1 R132 vs other brain
#   GSP #12: PRAD vs non-PRAD;                SPOP mutant vs other PRAD
# Samples: tumours in the flashier fit (shared_cache/gsp_loadings.tsv). The
# mutation contrasts use only samples with mutation calls (AlphaMissense-
# merged table); "other" = no call at the residue (IDH1) / in the gene (SPOP).
# PSI: LeafCutter PSI tables, {psi_root}/tcga/{chr}/{gene}/psi_table.tsv;
# headers cluster:start-end, matched on start or start-1 (as build_s19_cache.R).
# Out: note1/cache/psi_group_means.tsv. Any cwd; reads ~500 PSI tables -> SLURM.
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(yaml); library(data.table); library(tidyverse) })
FIG <- file.path(QROOT, "figures")
MD  <- file.path(QROOT, "mutation_ds")
OUT <- file.path(FIG, "supplementary/note1/cache"); dir.create(OUT, showWarnings = FALSE)

L <- fread(file.path(FIG, "shared_cache/gsp_loadings.tsv"), select = "rail_id")
meta <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), select = c("rail_id", "tcga_barcode", "gdc_cases.project.project_id"),
              data.table = FALSE) %>%
  transmute(rail_id = as.character(rail_id), cohort = gdc_cases.project.project_id,
            sampleId = vapply(strsplit(tcga_barcode, "-"), function(x) {
              x[4] <- sub("[A-Za-z]+$", "", x[4]); paste(x[1:4], collapse = "-") }, character(1))) %>%
  filter(rail_id %in% as.character(L$rail_id))
am <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), select = c("sampleId", "hugoGeneSymbol", "proteinChange"), data.table = FALSE)
called <- meta$sampleId %in% am$sampleId
r132 <- meta$sampleId %in% am$sampleId[am$hugoGeneSymbol == "IDH1" & grepl("R132", am$proteinChange)]
spop <- meta$sampleId %in% am$sampleId[am$hugoGeneSymbol == "SPOP"]
brain <- meta$cohort %in% c("TCGA-LGG", "TCGA-GBM")
prad  <- meta$cohort == "TCGA-PRAD"

GROUPS <- list(
  `1`  = list(`Brain (LGG + GBM)` = brain, `Non-brain` = !brain,
              `IDH1 R132 (brain)` = brain & r132, `Other brain` = brain & called & !r132),
  `12` = list(`PRAD` = prad, `Non-PRAD` = !prad,
              `SPOP mutant (PRAD)` = prad & spop, `Other PRAD` = prad & called & !spop))
for (g in names(GROUPS)) for (k in names(GROUPS[[g]]))
  cat(sprintf("[psi] GSP #%s %-20s n = %d\n", g, k, sum(GROUPS[[g]][[k]])))

match_cols <- function(grp, psi_cols) {
  se <- sub("^[^:]+:", "", psi_cols)
  a <- paste0(grp$start, "-", grp$end); b <- paste0(grp$start - 1, "-", grp$end)
  if (sum(a %in% se) == 0 && sum(b %in% se) == 0) return(rep(NA_integer_, nrow(grp)))
  match(if (sum(a %in% se) >= sum(b %in% se)) a else b, se)
}

res <- bind_rows(lapply(names(GROUPS), function(g) {
  tj <- fread(file.path(MD, "gsp_top", sprintf("fact_%s_top_junctions.tsv", g)), data.table = FALSE) %>%
    transmute(gene_name, chrom = chromosome, start = as.integer(start), end = as.integer(end), type)
  bind_rows(lapply(split(tj, list(tj$chrom, tj$gene_name), drop = TRUE), function(grp) {
    f <- file.path(file.path(QROOT, "leafcutter_analyses/psi_files"), "tcga", grp$chrom[1], grp$gene_name[1], "psi_table.tsv")
    if (!file.exists(f)) return(NULL)
    psi <- fread(f, data.table = FALSE); colnames(psi)[1] <- "rail_id"
    m <- match_cols(grp, colnames(psi)[-1]); keep <- !is.na(m)
    if (!any(keep)) return(NULL)
    rows <- match(meta$rail_id, as.character(psi$rail_id))
    bind_rows(lapply(which(keep), function(j) {
      v <- psi[[colnames(psi)[-1][m[j]]]][rows]
      tibble(gsp = as.integer(g), gene_name = grp$gene_name[j], start = grp$start[j], end = grp$end[j],
             type = grp$type[j], group = names(GROUPS[[g]]),
             mean_psi = vapply(GROUPS[[g]], function(s) mean(v[s], na.rm = TRUE), numeric(1)),
             n = vapply(GROUPS[[g]], function(s) sum(s & !is.na(v)), integer(1)))
    }))
  })) %>% { cat(sprintf("[psi] GSP #%s: %d / %d junctions with PSI\n", g,
                        n_distinct(paste(.$gene_name, .$start, .$end)), nrow(tj))); . }
}))
fwrite(res, file.path(OUT, "psi_group_means.tsv"), sep = "\t")
cat("[psi] wrote", file.path(OUT, "psi_group_means.tsv"), "\n")
