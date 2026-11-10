
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [matched, genome_fa, features_out, seq3_out, seq5_out]
matched_path  <- args[1]
genome_path   <- args[2]
features_out  <- args[3]
seq3_out      <- args[4]
seq5_out      <- args[5]

matched <- read_tsv(matched_path, show_col_types = FALSE)

dir.create(dirname(features_out), recursive = TRUE, showWarnings = FALSE)

if (nrow(matched) == 0) {
  write_tsv(tibble(), features_out)
  writeLines(character(0), seq3_out)
  writeLines(character(0), seq5_out)
  quit(save = "no", status = 0)
}

fasta <- FaFile(genome_path)
open(fasta)

# ---- intron sequence / GC / length ----
gr <- GRanges(seqnames = matched$chromosome,
              ranges   = IRanges(matched$start, matched$end),
              strand   = matched$strand)

seqs <- getSeq(fasta, gr)
gc <- letterFrequency(seqs, letters = c("G", "C"), as.prob = TRUE)
matched$gc_content <- rowSums(gc)
matched$length     <- matched$end - matched$start + 1

# ---- 3'ss and 5'ss sequences (MaxEntScan windows) ----
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

matched$seq_3ss <- mapply(get_3ss, matched$chromosome, matched$start,
                          matched$end, matched$strand)
matched$seq_5ss <- mapply(get_5ss, matched$chromosome, matched$start,
                          matched$end, matched$strand)

# ---- AG distance helpers (upstream / downstream of 3'ss) ----
get_downstream_AG <- function(chr, three_ss, strand, window = 100) {
  if (strand == "+") {
    r <- GRanges(chr, IRanges(three_ss + 1, three_ss + window))
    s <- as.character(getSeq(fasta, r)[[1]])
  } else {
    r <- GRanges(chr, IRanges(three_ss - window, three_ss - 1))
    s <- as.character(reverseComplement(getSeq(fasta, r)[[1]]))
  }
  m <- regexpr("AG", s, fixed = TRUE)
  if (m[1] == -1) NA_integer_ else as.integer(m[1] + 1)
}

get_upstream_AG <- function(chr, three_ss, strand, window = 100) {
  if (strand == "+") {
    r <- GRanges(chr, IRanges(three_ss - window, three_ss - 1))
    s <- as.character(getSeq(fasta, r)[[1]])
  } else {
    r <- GRanges(chr, IRanges(three_ss + 1, three_ss + window))
    s <- as.character(reverseComplement(getSeq(fasta, r)[[1]]))
  }
  matches <- gregexpr("AG", s, fixed = TRUE)[[1]]
  if (matches[1] == -1) NA_integer_
  else as.integer(window - matches[length(matches)])
}

matched <- matched %>%
  mutate(three_ss = if_else(strand == "+", end, start))

matched$downstream_AG_distance <- mapply(get_downstream_AG,
                                         matched$chromosome,
                                         matched$three_ss,
                                         matched$strand)
matched$upstream_AG_distance   <- mapply(get_upstream_AG,
                                         matched$chromosome,
                                         matched$three_ss,
                                         matched$strand)

write_tsv(matched, features_out)

writeLines(matched$seq_3ss, seq3_out)
writeLines(matched$seq_5ss, seq5_out)
