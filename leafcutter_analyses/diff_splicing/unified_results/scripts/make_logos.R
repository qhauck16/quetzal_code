
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(ggseqlogo)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [df_final, genome_fa, gene, tag, disrupted_out, mane_out]
df_path        <- args[1]
genome_path    <- args[2]
gene           <- toupper(args[3])
tag            <- args[4]
disrupted_out  <- args[5]
mane_out       <- args[6]

df <- read_tsv(df_path, show_col_types = FALSE)

fasta <- FaFile(genome_path)
open(fasta)

get_3ss <- function(chr, s, e, strand) {
  if (strand == "+") {
    r <- GRanges(chr, IRanges(e - 19, e + 3))
    as.character(getSeq(fasta, r)[[1]])
  } else {
    r <- GRanges(chr, IRanges(s - 3, s + 19))
    as.character(reverseComplement(getSeq(fasta, r)[[1]]))
  }
}

get_5ss <- function(chr, s, e, strand) {
  if (strand == "+") {
    r <- GRanges(chr, IRanges(s - 3, s + 5))
    as.character(getSeq(fasta, r)[[1]])
  } else {
    r <- GRanges(chr, IRanges(e - 5, e + 3))
    as.character(reverseComplement(getSeq(fasta, r)[[1]]))
  }
}

chr_col <- if ("chromosome" %in% names(df)) "chromosome" else "chrom"

dir.create(dirname(disrupted_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(mane_out),      recursive = TRUE, showWarnings = FALSE)

# ---- disrupted splice sites (the junction itself) ----
disrupted_plot <- function() {
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste0(gene, " ", tag, " - disrupted ss (no data)")))
  }
  seq3 <- mapply(get_3ss, df[[chr_col]], df$start, df$end, df$strand)
  seq5 <- mapply(get_5ss, df[[chr_col]], df$start, df$end, df$strand)

  keep3 <- nchar(seq3) == 23 & !grepl("[^ACGT]", seq3)
  keep5 <- nchar(seq5) == 9  & !grepl("[^ACGT]", seq5)

  l5 <- if (sum(keep5) > 0)
          ggseqlogo(seq5[keep5], method = "bits") +
            ggtitle(paste0(gene, " ", tag, " 5'ss disrupted (n = ", sum(keep5), ")"))
        else ggplot() + theme_void() + labs(title = "5'ss disrupted (no data)")
  l3 <- if (sum(keep3) > 0)
          ggseqlogo(seq3[keep3], method = "bits") +
            ggtitle(paste0(gene, " ", tag, " 3'ss disrupted (n = ", sum(keep3), ")"))
        else ggplot() + theme_void() + labs(title = "3'ss disrupted (no data)")
  l5 / l3
}

ggplot2::ggsave(disrupted_out, disrupted_plot(), width = 8, height = 6, dpi = 300)

# ---- MANE splice sites for the introns being disrupted ----
mane_plot <- function() {
  mane_df <- df %>%
    filter(type %in% c("alt 3'ss", "alt 5'ss"),
           !is.na(mane_intron_start), !is.na(mane_intron_end))

  if (nrow(mane_df) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste0(gene, " ", tag, " - MANE ss (no data)")))
  }

  seq3 <- mapply(get_3ss, mane_df[[chr_col]],
                 mane_df$mane_intron_start, mane_df$mane_intron_end,
                 mane_df$strand)
  seq5 <- mapply(get_5ss, mane_df[[chr_col]],
                 mane_df$mane_intron_start, mane_df$mane_intron_end,
                 mane_df$strand)

  keep3 <- nchar(seq3) == 23 & !grepl("[^ACGT]", seq3)
  keep5 <- nchar(seq5) == 9  & !grepl("[^ACGT]", seq5)

  l5 <- if (sum(keep5) > 0)
          ggseqlogo(seq5[keep5], method = "bits") +
            ggtitle(paste0(gene, " ", tag, " 5'ss MANE (n = ", sum(keep5), ")"))
        else ggplot() + theme_void() + labs(title = "5'ss MANE (no data)")
  l3 <- if (sum(keep3) > 0)
          ggseqlogo(seq3[keep3], method = "bits") +
            ggtitle(paste0(gene, " ", tag, " 3'ss MANE (n = ", sum(keep3), ")"))
        else ggplot() + theme_void() + labs(title = "3'ss MANE (no data)")
  l5 / l3
}

ggplot2::ggsave(mane_out, mane_plot(), width = 8, height = 6, dpi = 300)
