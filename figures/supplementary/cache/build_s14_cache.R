
# Cache for Supplementary Figure 15 (NKAP extended sashimi plots).
# Sections 1-5 of tcga/genome_wide/flashier/NKAP_detailed/custom_plot.R,
# looped over the top-15 events; the plotting half is redrawn in
# supplementary/figS15/build.R. One RDS per gene in cache/s14/, plus
# cache/s14/events.tsv (event order read by figS15/build.R).
#
# Self-contained since 2026-10-06 -- the two inputs that used to be static
# files in NKAP_detailed/ are rebuilt here:
#   events   <- NKAP_detailed/pick_events.R (config: alt 3'ss, cryptic
#               acceptor -12..-3 nt, rank by p.adjust then |dPSI|, gene must
#               have a PSI table; MANE-Select CDS/UTR class of the cryptic site)
#   samples  <- pygenome/scripts/select_samples.R (top 15 by k13 + 150
#               project-proportion-matched controls, set.seed(123)); verified
#               to give the same 15 + 150 rail_ids as the old NKAP/samples.tsv
#               when k13 comes from the quetzal_code NKAP table.
# Any cwd (absolute paths). Heavy (bigwig reads) -> SLURM.

if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
})
OUT_DIR  <- file.path(QROOT, "figures/supplementary/cache/s14")
TCGA_BW  <- file.path(QROOT, "inputs/bigwigs/tcga")
PSI_DIR  <- file.path(QROOT, "leafcutter_analyses/psi_files/tcga")
GENCODE  <- file.path(QROOT, "inputs/gencode_v46_granges.RDS")
MANE     <- file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt")
DF_FINAL <- file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results/nkap/20v60/dpsi0.05_fc2.5/df_final.tsv")
NKAP_TAB <- file.path(QROOT, "figures/figure3/cache/NKAP_full_plus_variants.tsv")
TCGA_META <- file.path(QROOT, "inputs/tcga_v2_samples.tsv")
N_EVENTS <- 15; DIST_MIN <- -12; DIST_MAX <- -3
N_COND <- 15; N_CTRL <- 150; CTRL_SEED <- 123
EXON_FLANK <- 50L; INTRON_FLANK <- 100L
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df_final <- fread(DF_FINAL, data.table = FALSE)
gencode  <- readRDS(GENCODE)
mane_tbl <- read_tsv(MANE, show_col_types = FALSE, comment = "")

# ---- samples (select_samples.R) -------------------------------------
meta <- fread(TCGA_META, data.table = FALSE) %>%
  dplyr::select(rail_id, tcga_barcode, project_id = gdc_cases.project.project_id) %>%
  mutate(project_short = str_split_i(project_id, "-", 2))
ranked <- fread(NKAP_TAB, data.table = FALSE) %>%
  transmute(rail_id = as.integer(rail_id), loading = k13) %>%
  arrange(desc(loading)) %>% inner_join(meta, by = "rail_id")
cond <- ranked %>% slice_head(n = N_COND)
props <- cond %>% count(project_id, name = "count") %>%
  mutate(n_control = round(count / sum(count) * N_CTRL))
set.seed(CTRL_SEED)
ctrl <- ranked %>% anti_join(cond, by = "rail_id") %>%
  inner_join(props %>% dplyr::select(project_id, n_control), by = "project_id") %>%
  group_by(project_id) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(dplyr::first(.x$n_control), nrow(.x)))) %>%
  ungroup() %>% dplyr::select(-n_control)
fwrite(bind_rows(cond %>% mutate(role = "condition"), ctrl %>% mutate(role = "control")),
       file.path(OUT_DIR, "samples.tsv"), sep = "\t")

# ---- events (pick_events.R) -----------------------------------------
mane_map <- setNames(mane_tbl$Ensembl_nuc, mane_tbl$symbol)
classify_region <- function(chr, pos, strand, gene) {
  tx <- mane_map[gene]
  if (is.na(tx)) return("no_mane_transcript")
  ents <- gencode[gencode$transcript_id == tx & gencode$type %in% c("CDS", "UTR")]
  if (length(ents) == 0) {
    tx_ids  <- unique(gencode$transcript_id)
    matched <- tx_ids[sub("\\.\\d+$", "", tx_ids) == sub("\\.\\d+$", "", tx)]
    if (length(matched) > 0)
      ents <- gencode[gencode$transcript_id %in% matched & gencode$type %in% c("CDS", "UTR")]
  }
  if (length(ents) == 0) return("no_cds_or_utr")
  ov <- findOverlaps(GRanges(chr, IRanges(pos, pos), strand = strand), ents, ignore.strand = TRUE)
  if (length(ov) == 0) return("intronic")
  if ("CDS" %in% as.character(ents$type[subjectHits(ov)])) "CDS" else "UTR"
}
cand <- df_final %>%
  filter(type == "alt 3'ss", !is.na(distance_to_cryptic),
         distance_to_cryptic >= DIST_MIN, distance_to_cryptic <= DIST_MAX) %>%
  arrange(p.adjust, desc(abs(deltapsi)))
cand <- cand[file.exists(file.path(PSI_DIR, cand$chrom, cand$gene_name, "psi_table.tsv")), ] %>%
  slice_head(n = N_EVENTS)
events <- cand %>%
  mutate(cryptic_pos   = if_else(strand == "+", as.integer(end), as.integer(start)),
         annotated_pos = if_else(strand == "+", as.integer(mane_intron_end), as.integer(mane_intron_start)),
         gene_region   = mapply(classify_region, chrom, cryptic_pos, strand, gene_name,
                                SIMPLIFY = TRUE, USE.NAMES = FALSE),
         rank = row_number()) %>%
  transmute(rank, gene = gene_name, chr = chrom, strand, distance_to_cryptic, gene_region,
            deltapsi, p_adjust = p.adjust, annotated_pos, cryptic_pos)
fwrite(events, file.path(OUT_DIR, "events.tsv"), sep = "\t")
cat(sprintf("[s14] %d events: %s\n", nrow(events), paste(events$gene, collapse = " ")))
samples <- bind_rows(cond %>% mutate(role = "condition"), ctrl %>% mutate(role = "control"))

for (i in seq_len(nrow(events))) {
  ev <- events[i, ]; GENE <- ev$gene
  chr <- ev$chr; strand <- ev$strand
  annot_acc <- ev$annotated_pos; cryp_acc <- ev$cryptic_pos
  shift_nt  <- ev$distance_to_cryptic
  if (strand == "-") { win_lo <- annot_acc - EXON_FLANK;   win_hi <- annot_acc + INTRON_FLANK
  } else             { win_lo <- annot_acc - INTRON_FLANK; win_hi <- annot_acc + EXON_FLANK }

  row_df <- df_final %>% filter(gene_name == GENE, type == "alt 3'ss",
                                distance_to_cryptic == shift_nt) %>% slice(1)
  donor_pos   <- if (strand == "-") as.integer(row_df$end) else as.integer(row_df$start)
  cryp_lo_hi  <- sort(c(as.integer(cryp_acc), donor_pos))
  annot_lo_hi <- sort(c(as.integer(annot_acc), donor_pos))

  region_gr <- GRanges(chr, IRanges(win_lo, win_hi)); n_bp <- width(region_gr)
  extract_cov <- function(rid, proj_short) {
    p <- file.path(TCGA_BW, proj_short, paste0(rid, ".bw"))
    if (!file.exists(p)) return(rep(NA_real_, n_bp))
    v <- tryCatch(as.numeric(import(p, which = region_gr, as = "NumericList")[[1]]),
                  error = function(e) NULL)
    if (is.null(v) || length(v) != n_bp) return(rep(NA_real_, n_bp))
    v
  }
  cond_mat <- do.call(rbind, mapply(extract_cov, cond$rail_id, cond$project_short, SIMPLIFY = FALSE))
  ctrl_mat <- do.call(rbind, mapply(extract_cov, ctrl$rail_id, ctrl$project_short, SIMPLIFY = FALSE))
  pos_vec <- seq(win_lo, win_hi)
  cov_df <- rbind(data.frame(pos = pos_vec, cov = colMeans(cond_mat, na.rm = TRUE), group = "cond"),
                  data.frame(pos = pos_vec, cov = colMeans(ctrl_mat, na.rm = TRUE), group = "ctrl"))

  psi <- fread(file.path(PSI_DIR, chr, GENE, "psi_table.tsv"), data.table = FALSE)
  setnames(psi, 1, "rail_id"); psi$rail_id <- as.integer(psi$rail_id)
  find_junc_col <- function(lo_hi) {
    hits <- grep(sprintf(":%d-%d$", lo_hi[1], lo_hi[2]), colnames(psi), value = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }
  cryp_col <- find_junc_col(cryp_lo_hi); annot_col <- find_junc_col(annot_lo_hi)
  group_psi <- function(rail_ids, col) {
    if (is.na(col)) return(NA_real_)
    mean(as.numeric(psi[[col]][match(rail_ids, psi$rail_id)]), na.rm = TRUE) * 100
  }
  psi_tbl <- data.frame(junction = c("cryptic", "annotated"),
                        acceptor = c(cryp_acc, annot_acc),
                        cond = c(group_psi(cond$rail_id, cryp_col), group_psi(cond$rail_id, annot_col)),
                        ctrl = c(group_psi(ctrl$rail_id, cryp_col), group_psi(ctrl$rail_id, annot_col)))

  mane_tx <- setNames(mane_tbl$Ensembl_nuc, mane_tbl$symbol)[GENE]
  tx_ents <- gencode[which(gencode$transcript_id == mane_tx)]
  if (length(tx_ents) == 0) {
    tx_ids  <- unique(gencode$transcript_id)
    matched <- tx_ids[sub("\\.\\d+$", "", tx_ids) == sub("\\.\\d+$", "", mane_tx)]
    tx_ents <- gencode[gencode$transcript_id %in% matched]
  }
  exons  <- as.data.frame(tx_ents[tx_ents$type == "exon"])[, c("start", "end")]
  codons <- as.data.frame(tx_ents[tx_ents$type %in% c("start_codon", "stop_codon")])
  codons <- data.frame(type = as.character(codons$type), pos = (codons$start + codons$end) / 2)

  saveRDS(list(gene = GENE, chr = chr, strand = strand, rank = ev$rank,
               shift = shift_nt, region = ev$gene_region,
               annot_acc = annot_acc, cryp_acc = cryp_acc, donor = donor_pos,
               win = c(win_lo, win_hi), cov = cov_df, psi = psi_tbl,
               exons = exons, codons = codons,
               n_valid = c(cond = sum(rowSums(!is.na(cond_mat)) > 0),
                           ctrl = sum(rowSums(!is.na(ctrl_mat)) > 0))),
          file.path(OUT_DIR, paste0(GENE, ".rds")))
  cat(sprintf("%s done: psi cond %s | ctrl %s\n", GENE,
              paste(round(psi_tbl$cond, 1), collapse = "/"),
              paste(round(psi_tbl$ctrl, 1), collapse = "/")))
}
