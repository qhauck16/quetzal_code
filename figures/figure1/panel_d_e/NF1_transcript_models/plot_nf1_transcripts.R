# One-off: plot the MANE NF1 transcript model and the NF1-216
# (ENST00000495910.6, an NMD isoform) directly from the GENCODE v46
# GRanges that powers the plotting_genes pipeline. Run interactively in
# RStudio with cwd = plotting_genes/NF1_transcript_models/.


if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggbio)
  library(GenomicRanges)
  library(rtracklayer)
})

# These two paths were relative in the original file; updated to
# absolute here so the script still runs from quetzal_code/. They must
# match `gencode_granges` / `mane_summary` in quetzal_code/paths.yaml.
GENCODE_RDS <- file.path(QROOT, "inputs/gencode_v46_granges.RDS")
MANE_TSV    <- file.path(QROOT, "inputs/MANE.GRCh38.v1.4.summary.txt")
OUT_DIR     <- "."

# ---- Load GENCODE + MANE map -----------------------------------------
gencode <- readRDS(GENCODE_RDS)
mane    <- fread(MANE_TSV, data.table = FALSE)

nf1_exons <- gencode[gencode$gene_name == "NF1" & gencode$type == "exon"]
cat(sprintf("NF1 exon rows in GENCODE: %d (across %d transcripts)\n",
            length(nf1_exons), length(unique(nf1_exons$transcript_id))))

# ---- MANE NF1 transcript ---------------------------------------------
mane_nf1_tx     <- mane$Ensembl_nuc[mane$symbol == "NF1"]
mane_exons      <- nf1_exons[nf1_exons$transcript_id %in% mane_nf1_tx]
mane_tx_models  <- rtracklayer::split(mane_exons, mane_exons$transcript_id)
cat(sprintf("MANE NF1 transcript_id: %s (%d exons)\n",
            paste(unique(mane_exons$transcript_id), collapse = ", "),
            length(mane_exons)))

# ---- NF1-216 (ENST00000495910.6) -------------------------------------
nf1_216_id      <- "ENST00000495910.6"
nf1_216_exons   <- nf1_exons[nf1_exons$transcript_id == nf1_216_id]
stopifnot("ENST00000495910.6 not found in GENCODE NF1 exons" =
            length(nf1_216_exons) > 0)
nf1_216_models  <- rtracklayer::split(nf1_216_exons, nf1_216_exons$transcript_id)
cat(sprintf("NF1-216 (%s) exons: %d\n",
            nf1_216_id, length(nf1_216_exons)))

# ---- Compute shared x-axis range ------------------------------------
# Union of both transcripts' coordinates. Every plot below is forced to
# this exact range with `expand = FALSE`, so a 1 bp shift in genomic
# coords maps to the same number of inches in every PDF -- you can stack
# the saved PDFs in Illustrator/Inkscape and the exons and gap arrows
# line up to the bp.
all_exons <- c(mane_exons, nf1_216_exons)
X_MIN <- min(start(all_exons))
X_MAX <- max(end(all_exons))
cat(sprintf("Shared x-range: %d - %d (%d kb)\n",
            X_MIN, X_MAX, round((X_MAX - X_MIN) / 1000)))

# ---- Plot helper -----------------------------------------------------
# `exon.rect.h` fixed across plots so exon thickness is identical; same
# `coord_cartesian` range so genomic coordinates align between PDFs.
EXON_H   <- 0.05
PDF_W_IN <- 10                            # match Illustrator paper width
PDF_H_IN <- 1                             # per-transcript height
plot_tx <- function(tx_models, title) {
  ggplot() +
    ggbio::geom_alignment(tx_models, gap.geom = getOption("figure1.gap_geom", "arrow"),
                          label = FALSE, exon.rect.h = EXON_H) +
    scale_x_continuous(
      limits = c(X_MIN, X_MAX),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(X_MIN, X_MAX), expand = FALSE, clip = "off") +
    theme_minimal() +
    theme(plot.margin = margin(5, 5, 5, 5)) +
    ggtitle(title)
}

mane_p   <- plot_tx(mane_tx_models, "NF1 (MANE transcript)")
nf1216_p <- plot_tx(nf1_216_models, "NF1-216 (ENST00000495910.6, NMD)")

# ---- Combined plot: both transcripts stacked, same coord system -----
combined_models <- rtracklayer::split(
  all_exons,
  c(mane_exons$transcript_id, nf1_216_exons$transcript_id)
)
combined_p <- ggplot() +
  ggbio::geom_alignment(combined_models, gap.geom = getOption("figure1.gap_geom", "arrow"),
                        label = TRUE, exon.rect.h = EXON_H) +
  scale_x_continuous(
    limits = c(X_MIN, X_MAX),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(X_MIN, X_MAX), expand = FALSE, clip = "off") +
  theme_minimal() +
  theme(plot.margin = margin(5, 5, 5, 5)) +
  ggtitle("NF1: MANE vs NF1-216 (ENST00000495910.6, NMD)")

# ---- Save -----------------------------------------------------------
# All single-transcript PDFs share PDF_W_IN width and PDF_H_IN height so
# bp -> inches is identical; combined plot is 2x tall.
ggplot2::ggsave(file.path(OUT_DIR, "NF1_MANE.pdf"),
                mane_p,   width = PDF_W_IN, height = PDF_H_IN, device = cairo_pdf)
ggplot2::ggsave(file.path(OUT_DIR, "NF1-216.pdf"),
                nf1216_p, width = PDF_W_IN, height = PDF_H_IN, device = cairo_pdf)
ggplot2::ggsave(file.path(OUT_DIR, "NF1_MANE_and_216.pdf"),
                combined_p, width = PDF_W_IN, height = 2 * PDF_H_IN,
                device = cairo_pdf)

cat(sprintf("Wrote %s, %s, %s\n",
            file.path(OUT_DIR, "NF1_MANE.pdf"),
            file.path(OUT_DIR, "NF1-216.pdf"),
            file.path(OUT_DIR, "NF1_MANE_and_216.pdf")))
