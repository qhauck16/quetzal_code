# GSP #30 (YLPM1) coverage-track inputs for Fig. 4b (KDM4C, TRIM38): top-13 GSP #30 samples
# vs project-matched controls (set.seed(123)). Extracted from poisson_NMF_junction_plots.Rmd
# (chunk "prep_pygenometrack" and its KDM4C / TRIM38 calls), code unchanged apart from paths.
# Writes fig4b_tracks/<GENE>/{<GENE>_sample_bigwigs.tsv, <GENE>_groups_file.txt,
# junctions.bed}; then run gsp_tracks.sh, which averages the recount3 bigwigs per group
# (GenometracksByGenotype AggregateBigwigsForPlotting.py) into the High_k30 / Matched_Controls
# .bw and MeanPSI .links files the figure reads.
#   Rscript gsp_tracks.R [OUT]           (run from coverage_tracks/; OUT default "fig4b_tracks")
suppressPackageStartupMessages({ library(tidyverse); library(readr); library(data.table) })
args <- commandArgs(TRUE); OUT <- if (length(args)) args[1] else "fig4b_tracks"
psi_dir <- '../leafcutter_analyses/psi_files'

prep_pygenometrack <- function(gene_name, chr, strand, region,
                               n_top = 13, n_control_mult = 10,
                               out_dir = OUT) {

  psi_path <- paste0(psi_dir, '/tcga/', chr, '/', gene_name, '/psi_table.tsv')
  pygenome_outdir <- file.path(out_dir, gene_name)
  dir.create(pygenome_outdir, showWarnings = FALSE, recursive = TRUE)

  # Top n_top YLPM1 samples by k30
  set.seed(123)
  YLPM1_full_info <- read_tsv('../figures/figure4/cache/YLPM1_full_plus_variants.tsv',
         # mutation cols are ~all NA: readr would guess logical and
         # silently NA every call (2026-10-05)
         col_types = cols(.default = col_guess(), gene_mutation = col_character(), wxs_mutation = col_character(), rnaseq_mutation = col_character()), show_col_types = FALSE) %>%
    arrange(desc(k30))

  condition <- YLPM1_full_info %>%
    slice_head(n = n_top) %>%
    pull(rail_id)

  # Project-matched controls
  project_props <- YLPM1_full_info %>%
    slice_head(n = n_top) %>%
    group_by(project_id) %>%
    summarize(count = n()) %>%
    mutate(prop = count / sum(count),
           n_control = round(prop * n_top * n_control_mult))

  control <- YLPM1_full_info %>%
    filter(!rail_id %in% condition) %>%
    inner_join(project_props, by = "project_id") %>%
    group_by(project_id) %>%
    group_modify(~ {
      n <- dplyr::first(.x$n_control)
      dplyr::slice_sample(.x, n = min(n, nrow(.x)))
    }) %>%
    ungroup() %>%
    pull(rail_id)

  sample_df_input <- YLPM1_full_info %>%
    filter(rail_id %in% c(condition, control))

  # Sample bigwigs
  tcga_bw_dir <- "../inputs/bigwigs/tcga"
  all_samples <- c(condition, control)
  group_labels <- c(rep("High_k30", length(condition)),
                    rep("Matched_Controls", length(control)))

  sample_bigwigs <- tibble(
    sample = all_samples,
    bigwigPath = file.path(
      tcga_bw_dir,
      str_split_i(sample_df_input$project_id[match(all_samples, sample_df_input$rail_id)], '-', 2),
      paste0(all_samples, ".bw")
    ),
    Group_label = group_labels,
    strand = "."
  )

  write.table(sample_bigwigs,
              file = file.path(pygenome_outdir, paste0(gene_name, "_sample_bigwigs.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Groups file
  groups_file <- tibble(
    Group_label = c("High_k30", "Matched_Controls"),
    color = c("#7fc97f", "#969696"),
    Bedgz = "", Supergroup = ""
  )

  write.table(groups_file,
              file = file.path(pygenome_outdir, paste0(gene_name, "_groups_file.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # PSI BED
  make_bed_df <- function(mat, chr, strand) {
    df_vals <- as.data.frame(t(mat * 100))
    df_vals$junc_id <- rownames(df_vals)

    df_meta <- df_vals %>%
      separate(junc_id, into = c("clu", "range"), sep = ":", remove = FALSE) %>%
      separate(range, into = c("start", "end"), sep = "-", convert = TRUE)

    df_meta$start <- df_meta$start - 1
    df_meta$`#Chrom` <- chr
    df_meta$strand <- strand

    df_meta$junc <- paste0(df_meta$`#Chrom`, ":", df_meta$start, ":", df_meta$end,
                           ":clu_", df_meta$clu, "_", df_meta$strand)
    df_meta$gid  <- paste0(df_meta$`#Chrom`, "_clu_", df_meta$clu, "_", df_meta$strand)

    df_meta %>%
      relocate(`#Chrom`, start, end, junc, gid, strand) %>%
      dplyr::select(-clu, -junc_id, -junc)
  }

  suppressWarnings({
    tcga_norm_dat <- fread(psi_path, data.table = F) %>%
      filter(V1 %in% sample_bigwigs$sample) %>%
      column_to_rownames('V1')
  })
  tcga_norm_dat[is.na(tcga_norm_dat)] <- 0

  tcga_bed <- make_bed_df(tcga_norm_dat, chr, strand) %>%
    group_by(start, end) %>%
    slice_head(n = 1) %>%
    ungroup()

  num_cols <- sapply(tcga_bed, is.numeric)
  tcga_bed[, num_cols] <- lapply(tcga_bed[, num_cols, drop = FALSE],
                                 function(x) round(as.double(x), 4))

  write.table(tcga_bed,
              file = file.path(pygenome_outdir, "junctions.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Print terminal commands
  message("Files written to ", pygenome_outdir, "/")
  message("Run these commands in the terminal, then call make_sashimi_ini():\n")

  cat(paste0(
    "(head -n 1 ", pygenome_outdir, "/junctions.bed && ",
    "tail -n +2 ", pygenome_outdir, "/junctions.bed | sort -k1,1 -k2,2n) ",
    "| bgzip -c > ", pygenome_outdir, "/junctions.bed.gz\n"
  ))

  cat(paste0("tabix -p bed ", pygenome_outdir, "/junctions.bed.gz\n"))

  cat(paste0(
    "python ~/software/GenometracksByGenotype/AggregateBigwigsForPlotting.py ",
    "--BigwigList ", file.path(pygenome_outdir, paste0(gene_name, "_sample_bigwigs.tsv")),
    " --BigwigListType KeyFile ",
    " --OutputPrefix ", pygenome_outdir, "/",
    " --Region ", region,
    " --GroupSettingsFile ", file.path(pygenome_outdir, paste0(gene_name, "_groups_file.txt")),
    " --Bed12GenesToIni ~/software/GenometracksByGenotype/PremadeTracks/gencode.v26.FromGTEx.genes.bed12.gz",
    " --Normalization WholeGenome",
    " --BedfileForSashimiLinks ", pygenome_outdir, "/junctions.bed.gz\n"
  ))

  cat(paste0(
    "pyGenomeTracks --region ", region,
    " --tracks ", pygenome_outdir, "/tracks.ini",
    " -o ", pygenome_outdir, "/simple_coverage.pdf\n"
  ))

  invisible(list(pygenome_outdir = pygenome_outdir, region = region))
}

prep_pygenometrack('TRIM38', 'chr6', '+', 'chr6:25,973,077-25,984,559')
prep_pygenometrack('KDM4C', 'chr9', '+', 'chr9:7,011,426-7,014,212')
