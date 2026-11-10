# Variant parsing for figure 3 panels a and b.
#
# NKAP_full_plus_variants.tsv carries mutations as free text in three
# columns (gene_mutation, wxs_mutation, rnaseq_mutation), each a
# ";"-separated list of "GENE:change1,change2" blocks. Nothing is
# pre-parsed, so both panels have to dig positions and labels out of
# those strings.
#
# Only ONE source wins per row in gene_mutation, so a WXS PRPF8 hit can
# be hidden behind a FAM50A AlphaMissense call. Searching all three
# columns prevents that masking.

# Every amino-acid position for `gene` in one mutation string.
extract_all_positions <- function(mutation_string, gene) {
  if (is.na(mutation_string) || mutation_string == "" ||
      mutation_string == "NA") return(numeric(0))
  # Whole gene block up to the next ";", then split on "," so every
  # mutation in the block contributes a position -- not just the first.
  blocks <- regmatches(mutation_string,
                       gregexpr(paste0(gene, ":[^;]+"), mutation_string))[[1]]
  if (length(blocks) == 0) return(numeric(0))
  pos <- unlist(lapply(blocks, function(b) {
    muts <- strsplit(sub(paste0("^", gene, ":"), "", b), ",")[[1]]
    suppressWarnings(as.numeric(regmatches(muts, regexpr("\\d+", muts))))
  }))
  pos[!is.na(pos)]
}

# ... across all three mutation sources on a row.
gene_positions <- function(row, gene) {
  c(extract_all_positions(row[["gene_mutation"]],   gene),
    extract_all_positions(row[["wxs_mutation"]],    gene),
    extract_all_positions(row[["rnaseq_mutation"]], gene))
}

# Every mutation STRING (not position) for `gene` on a row.
gene_changes <- function(row, gene) {
  out <- character(0)
  for (col in c("gene_mutation", "wxs_mutation", "rnaseq_mutation")) {
    s <- row[[col]]
    if (is.na(s) || s == "" || s == "NA") next
    for (b in regmatches(s, gregexpr(paste0(gene, ":[^;]+"), s))[[1]])
      out <- c(out, strsplit(sub(paste0("^", gene, ":"), "", b), ",")[[1]])
  }
  unique(trimws(out[out != ""]))
}

# Priority windows. NKAP's functional region is AA 300-413 and PRPF8's
# is 1570-1900; FAM50A and SLU7 count anywhere. Priority order matters:
# a sample with both an NKAP and a PRPF8 hit is called NKAP.
VARIANT_SPEC <- list(
  list(label = "NKAP variant",   gene = "NKAP",   lo = 300,  hi = 413),
  list(label = "PRPF8 variant",  gene = "PRPF8",  lo = 1570, hi = 1900),
  list(label = "FAM50A variant", gene = "FAM50A", lo = -Inf, hi = Inf),
  list(label = "SLU7 variant",   gene = "SLU7",   lo = -Inf, hi = Inf)
)
VARIANT_OTHER <- "No NKAP/PRPF8/FAM50A/SLU7 variant"   # same gene order as the title and the other legend entries

# Returns a data frame with `variant_class` and `proteinChange`.
classify_variants <- function(df) {
  cls <- character(nrow(df))
  lab <- character(nrow(df))
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    hit <- NULL
    for (spec in VARIANT_SPEC) {
      p <- gene_positions(row, spec$gene)
      if (any(p >= spec$lo & p <= spec$hi)) { hit <- spec; break }
    }
    if (is.null(hit)) { cls[i] <- VARIANT_OTHER; lab[i] <- ""; next }
    cls[i] <- hit$label
    cand <- gene_changes(row, hit$gene)
    if (length(cand) == 0) { lab[i] <- ""; next }
    pos <- suppressWarnings(as.numeric(regmatches(cand, regexpr("\\d+", cand))))
    keep <- which(!is.na(pos) & pos >= hit$lo & pos <= hit$hi)
    # Fall back to every candidate when none is in range, so a row is
    # never labelled blank just because its position parse failed.
    lab[i] <- paste(cand[if (length(keep)) keep else seq_along(cand)],
                    collapse = ", ")
  }
  df$variant_class  <- cls
  df$proteinChange  <- lab
  df
}

# Legend order = priority order, "Other" last.
VARIANT_LEVELS <- c(vapply(VARIANT_SPEC, `[[`, character(1), "label"),
                    VARIANT_OTHER)
# Variant-class colours.
VARIANT_COL <- setNames(
  c("#8E44AD", "#E74C3C", "#3498DB", "#27AE60", "grey65"),
  VARIANT_LEVELS)
