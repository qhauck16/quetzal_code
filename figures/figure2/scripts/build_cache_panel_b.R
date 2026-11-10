.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

# Cache builder for figure 2, panel b (neuronal GSP ECDF).
#
# Faithful extraction of the upstream chain in
#   tcga/genome_wide/flashier/neuronal_gsp.Rmd
# chunks `loading re-used files`, `finding top junctions and nearby
# annotated junctions`, and `GTEx and TCGA PSI`.
#
# Produces two caches so the ECDF can be replotted without re-running
# any of it:
#   panel_b_top_junc.tsv    <- top_junc_df   (first loop; needs the 3.5 G RDS)
#   panel_b_running_psi.tsv <- running_psi_df (second loop; 2 big PSI tables/gene)
#
# The second loop is the slow one and is RESUMABLE: it appends one row
# per gene and skips genes already present on a re-run.
#
# NOTE: this reproduces the Rmd's behaviour exactly, including quirks --
# do not "fix" them here or the cache stops matching the published panel.
# See the report accompanying this script for the known quirks.
#
#   Rscript --vanilla build_cache_panel_b.R

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(flashier)
  library(fastTopics)
  library(GenomicRanges)
})

# S4Vectors (via GenomicRanges) masks dplyr::rename and friends -- a bare
# rename() here fails with "object 'V1' not found" on every gene. Bind
# the dplyr versions explicitly.
rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter

# ------------------------------------------------------------------
# Absolute paths (runnable from any cwd)
# ------------------------------------------------------------------
FLASH_RDS  <- file.path(QROOT, "quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS")
TCGA_META  <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
GTEX_META  <- file.path(QROOT, "inputs/gtexv2_metadata.tsv")
GENCODE    <- file.path(QROOT, "inputs/gencode_v46_granges.RDS")
OBJ_DIR    <- file.path(QROOT, "quetzal_runs/genome_wide")
ALL_GENES  <- file.path(QROOT, "inputs/snaptron/tcga")
COLS_TMPL  <- file.path(QROOT, "inputs/snaptron/tcga/chr1/snaptron_output/ABCA4_snaptron.tsv")
PSI_DIR    <- file.path(QROOT, "leafcutter_analyses/psi_files")
CACHE_DIR  <- file.path(QROOT, "figures/figure2/cache")

TOP_JUNC_TSV <- file.path(CACHE_DIR, 'panel_b_top_junc.tsv')
RUN_PSI_TSV  <- file.path(CACHE_DIR, 'panel_b_running_psi.tsv')

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

CUTOFF <- 1   # norm_F[, fact] > cutoff
FACT   <- 1   # factor / GSP index (GSP #1)

# ------------------------------------------------------------------
# Metadata (cheap, needed by loop 2 regardless)
# ------------------------------------------------------------------
cat('[panel_b] reading metadata\n')
tcga_metadata <- fread(TCGA_META, data.table = FALSE)
gtex_metadata <- fread(GTEX_META, data.table = FALSE)

# ------------------------------------------------------------------
# LOOP 1 -- top_junc_df
# Skipped entirely (including the 3.5 G RDS load) when the cache is
# present and newer than the RDS it was derived from.
# ------------------------------------------------------------------
top_junc_fresh <- file.exists(TOP_JUNC_TSV) &&
  file.info(TOP_JUNC_TSV)$mtime > file.info(FLASH_RDS)$mtime

if (top_junc_fresh) {
  cat('[panel_b] top_junc cache is newer than the RDS -- skipping loop 1\n')
  top_junc_df <- fread(TOP_JUNC_TSV, data.table = FALSE)
  cat(sprintf('[panel_b] loaded %d cached top junctions\n', nrow(top_junc_df)))
} else {
  cat('[panel_b] building top_junc_df (loads the 3.5 G RDS; slow)\n')
  gencode <- readRDS(GENCODE)
  out <- readRDS(FLASH_RDS)

  ldf <- flashier::ldf(out, type = 'i')
  # absorb the D into the F matrix
  norm_L <- ldf$L
  norm_F <- t(diag(ldf$D) %*% t(ldf$F))
  rm(out, ldf); gc()

  facts_to_analyze <- sum(norm_F[, FACT] > CUTOFF)
  top_poisson_facts <- rownames(norm_F)[order(norm_F[, FACT], decreasing = TRUE)][1:facts_to_analyze]
  cat(sprintf('[panel_b] %d junction groups pass norm_F[, %d] > %g\n',
              facts_to_analyze, FACT, CUTOFF))

  cols <- c('gene_name', colnames(fread(COLS_TMPL)))
  top_junc_df <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)

  # `strand` persists across iterations in the Rmd; seed it so the
  # length() guard below cannot error on the first gene.
  strand <- character(0)
  n_tot <- length(top_poisson_facts)
  n_fail <- 0L

  for (idx in seq_along(top_poisson_facts)) {
    pnmf_factor <- top_poisson_facts[idx]
    gene_name <- str_split_i(pnmf_factor, '\\.', 1)
    cat(sprintf('[panel_b][loop1][%d/%d] %s\n', idx, n_tot, gene_name))

    ok <- tryCatch({
      pnmf_fact_num <- str_split_i(pnmf_factor, '\\.', 2)

      for (i in c(1:22, 'X')) {
        if (file.exists(paste0(OBJ_DIR, '/chr', i, '/FastTopics_output/', gene_name, '/res.RDS'))) {
          res <- readRDS(paste0(OBJ_DIR, '/chr', i, '/FastTopics_output/', gene_name, '/res.RDS'))
          dat <- fread(paste0(OBJ_DIR, '/chr', i, '/FastTopics_output/', gene_name, '/dat.tsv'), data.table = FALSE)
          chr <- i
          gene_table <- fread(paste0(ALL_GENES, '/chr', chr, '/snaptron_output/', gene_name, '_snaptron.tsv'))

          gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c('exon')]
          strand <- unique(as.character(strand(gene_data)))
          next
        }
      }
      if (length(strand) == 0) {
        TRUE   # Rmd `next`s; nothing appended
      } else {
        res <- poisson2multinom(res)

        unn_juncs <- gene_table %>%
          filter(annotated == 0) %>%
          transmute(junc = paste0(start, '-', end)) %>%
          unlist() %>%
          as.vector()

        # unannotated junc with highest weighting in factor
        top_junc <- rownames(res$F[rownames(res$F) %in% unn_juncs, pnmf_fact_num, drop = FALSE])[
          which.max(res$F[rownames(res$F) %in% unn_juncs, pnmf_fact_num])]

        junc_start <- as.numeric(str_split_i(top_junc, '-', 1))
        junc_end   <- as.numeric(str_split_i(top_junc, '-', 2))

        top_junc_info <- as.vector(unlist(filter(gene_table, start == junc_start & end == junc_end)))
        to_keep <- c(gene_name, top_junc_info)
        names(to_keep) <- cols
        to_keep <- as.data.frame(as.list(to_keep))

        top_junc_df <- rbind(top_junc_df, to_keep)
        TRUE
      }
    }, error = function(e) {
      cat(sprintf('[panel_b][loop1][FAIL] %s: %s\n', gene_name, conditionMessage(e)))
      FALSE
    })
    if (!isTRUE(ok)) n_fail <- n_fail + 1L
  }

  cat(sprintf('[panel_b][loop1] done: %d rows, %d failures\n', nrow(top_junc_df), n_fail))
  write.table(top_junc_df, TOP_JUNC_TSV, sep = '\t', quote = FALSE, row.names = FALSE)

  # Free the big objects before loop 2 -- other cache scripts may be
  # running alongside this one.
  rm(gencode, norm_L, norm_F, top_poisson_facts); gc()
}

# ------------------------------------------------------------------
# LOOP 2 -- running_psi_df (resumable, appends per gene)
# ------------------------------------------------------------------
psi_cols <- c('gene_name', 'TCGA_brain_psi', 'TCGA_other_psi',
              'GTEx_brain_psi', 'GTEx_other_psi')

done_genes <- character(0)
if (file.exists(RUN_PSI_TSV)) {
  prev <- fread(RUN_PSI_TSV, data.table = FALSE)
  done_genes <- as.character(prev$gene_name)
  cat(sprintf('[panel_b] resuming: %d genes already in %s\n',
              length(done_genes), basename(RUN_PSI_TSV)))
} else {
  write.table(setNames(data.frame(matrix(ncol = length(psi_cols), nrow = 0)), psi_cols),
              RUN_PSI_TSV, sep = '\t', quote = FALSE, row.names = FALSE)
}

n_tot  <- nrow(top_junc_df)
n_fail <- 0L
n_new  <- 0L

for (i in seq_len(n_tot)) {
  gene_name <- top_junc_df$gene_name[i]
  if (gene_name %in% done_genes) {
    cat(sprintf('[panel_b][loop2][%d/%d] %s -- cached, skipping\n', i, n_tot, gene_name))
    next
  }
  cat(sprintf('[panel_b][loop2][%d/%d] %s\n', i, n_tot, gene_name))

  ok <- tryCatch({
    chrom      <- top_junc_df$chromosome[i]
    junc_start <- top_junc_df$start[i]
    junc_end   <- top_junc_df$end[i]

    gtex_path <- paste0(PSI_DIR, '/gtex/', chrom, '/', gene_name, '/psi_table.tsv')
    tcga_path <- paste0(PSI_DIR, '/tcga/', chrom, '/', gene_name, '/psi_table.tsv')

    if (!file.exists(gtex_path) || !file.exists(tcga_path) ||
        file.size(gtex_path) < 100 || file.size(tcga_path) < 100) {
      TRUE   # Rmd `next`
    } else {
      suppressWarnings({
        gtex_psi_table <- fread(gtex_path, data.table = FALSE) %>% rename(rail_id = V1)
        tcga_psi_table <- fread(tcga_path, data.table = FALSE) %>% rename(rail_id = V1)
      })

      junc_colnames <- intersect(
        colnames(tcga_psi_table)[grepl(paste0(junc_start, '-', junc_end), colnames(tcga_psi_table))],
        colnames(gtex_psi_table)[grepl(paste0(junc_start, '-', junc_end), colnames(gtex_psi_table))]
      )[1]

      if (is.na(junc_colnames)) {
        TRUE   # Rmd `next`
      } else {
        tcga_to_plot <- tcga_psi_table[, c('rail_id', junc_colnames)] %>%
          left_join(dplyr::select(tcga_metadata, c(rail_id, gdc_cases.project.project_id)), by = 'rail_id')
        colnames(tcga_to_plot) <- c('rail_id', 'junc_psi', 'gdc_cases.project.project_id')

        gtex_to_plot <- gtex_psi_table[, c('rail_id', junc_colnames)] %>%
          left_join(dplyr::select(gtex_metadata, c(rail_id, study)), by = 'rail_id')
        colnames(gtex_to_plot) <- c('rail_id', 'junc_psi', 'gtex_study')

        tcga_stats <- tcga_to_plot %>%
          mutate(group_id = ifelse(gdc_cases.project.project_id %in% c('TCGA-GBM', 'TCGA-LGG'),
                                   'TCGA: Brain Cohort', 'TCGA: Other')) %>%
          group_by(group_id) %>%
          summarize(mean_psi = mean(junc_psi, na.rm = TRUE))

        gtex_stats <- gtex_to_plot %>%
          mutate(group_id = ifelse(gtex_study %in% c('BRAIN'),
                                   'GTEx: Brain Cohort', 'GTEx: Other')) %>%
          group_by(group_id) %>%
          summarize(mean_psi = mean(junc_psi, na.rm = TRUE))

        new_row <- tibble::tibble(
          gene_name      = gene_name,
          TCGA_brain_psi = tcga_stats$mean_psi[tcga_stats$group_id == 'TCGA: Brain Cohort'],
          TCGA_other_psi = tcga_stats$mean_psi[tcga_stats$group_id == 'TCGA: Other'],
          GTEx_brain_psi = gtex_stats$mean_psi[gtex_stats$group_id == 'GTEx: Brain Cohort'],
          GTEx_other_psi = gtex_stats$mean_psi[gtex_stats$group_id == 'GTEx: Other']
        )

        # Append immediately so a killed run loses at most one gene.
        write.table(new_row[, psi_cols], RUN_PSI_TSV, sep = '\t', quote = FALSE,
                    row.names = FALSE, col.names = FALSE, append = TRUE)
        n_new <- n_new + 1L

        rm(gtex_psi_table, tcga_psi_table, tcga_to_plot, gtex_to_plot)
        TRUE
      }
    }
  }, error = function(e) {
    cat(sprintf('[panel_b][loop2][FAIL] %s: %s\n', gene_name, conditionMessage(e)))
    FALSE
  })
  if (!isTRUE(ok)) n_fail <- n_fail + 1L
  if (i %% 25 == 0) gc()
}

cat(sprintf('[panel_b][loop2] done: %d new rows, %d failures\n', n_new, n_fail))
cat(sprintf('[panel_b] wrote %s\n', TOP_JUNC_TSV))
cat(sprintf('[panel_b] wrote %s\n', RUN_PSI_TSV))
