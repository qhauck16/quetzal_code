#!/usr/bin/env Rscript
# Lift GRCh37 mutation rows to GRCh38 and merge pan_can_atlas + gdc into one file.
# Output: merged_grch38.tsv (canonical) + stats logged to stdout.

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
})

mut_dir   <- "mutations"
work_dir  <- "."   # run from mutation_processing/
chain_gz  <- "../inputs/alphamissense/hg19ToHg38.over.chain.gz"
chain_un  <- sub("\\.gz$", "", chain_gz)
liftover  <- Sys.getenv("LIFTOVER", Sys.which("liftOver"))   # UCSC liftOver binary
out_path  <- file.path(work_dir, "merged_grch38.tsv")

stopifnot(file.exists(chain_gz), file.exists(liftover))
if (!file.exists(chain_un)) system2("gunzip", c("-k", chain_gz))

# ---------- read inputs ----------
cat("Reading inputs ...\n")
pc  <- fread(file.path(mut_dir, "pan_can_atlas/compiled.tsv"))
gdc <- fread(file.path(mut_dir, "gdc/compiled.tsv"))
pc[,  source := "pan_can_atlas"]
gdc[, source := "gdc"]
if (!"normalRefCount" %in% names(gdc)) gdc[, normalRefCount := NA_integer_]

keep_cols <- c("hugoGeneSymbol","entrezGeneId","sampleId","patientId","studyId",
  "center","mutationStatus","validationStatus","tumorAltCount","tumorRefCount",
  "normalRefCount","startPosition","endPosition","referenceAllele","variantAllele",
  "proteinChange","mutationType","variantType","keyword","ncbiBuild","chr",
  "refseqMrnaId","proteinPosStart","proteinPosEnd","source")
pc  <- pc[,  ..keep_cols]
gdc <- gdc[, ..keep_cols]

mut <- rbindlist(list(pc, gdc), use.names = TRUE)
mut[, .row_id := .I]
cat(sprintf("Loaded %d rows  (pan_can_atlas=%d, gdc=%d)\n",
            nrow(mut), nrow(pc), nrow(gdc)))
cat("\nrows by source x ncbiBuild:\n"); print(mut[, .N, by = .(source, ncbiBuild)])

# ---------- normalize chromosome ----------
mut[, chr := as.character(chr)]
mut[chr == "23", chr := "X"]
mut[chr == "24", chr := "Y"]
mut[, chr_ucsc := paste0("chr", chr)]

g37 <- mut[ncbiBuild == "GRCh37"]
g38 <- mut[ncbiBuild == "GRCh38"]
cat(sprintf("\nGRCh37 to lift: %d   |   already GRCh38: %d\n", nrow(g37), nrow(g38)))

# ---------- liftover ----------
bed_in  <- file.path(work_dir, "tmp_in.bed")
bed_out <- file.path(work_dir, "tmp_out.bed")
bed_unm <- file.path(work_dir, "tmp_unmapped.bed")

fwrite(
  data.table(
    chrom  = g37$chr_ucsc,
    start  = g37$startPosition - 1L,
    end    = g37$endPosition,
    name   = g37$.row_id,
    score  = 0L,
    strand = "+"),
  bed_in, sep = "\t", col.names = FALSE)

cat("\nRunning liftOver ...\n")
rc_run <- system2(liftover,
  args = c("-bedPlus=6", "-tab", bed_in, chain_un, bed_out, bed_unm))
stopifnot(rc_run == 0)

mapped <- fread(bed_out, header = FALSE,
  col.names = c("chrom","start","end","name","score","strand"))
n_unm <- length(readLines(bed_unm))   # blocks are 2-line entries; just rough
cat(sprintf("Mapped: %d / %d (%.3f%%)\n",
            nrow(mapped), nrow(g37), 100 * nrow(mapped)/nrow(g37)))

cat("\nstrand distribution post-lift:\n"); print(table(mapped$strand))
canon <- paste0("chr", c(1:22, "X", "Y"))
n_nc  <- sum(!mapped$chrom %in% canon)
cat(sprintf("non-canonical contigs post-lift: %d (these will be dropped)\n", n_nc))
mapped <- mapped[chrom %in% canon]

# join lifted coords back
g37_lift <- merge(g37, mapped, by.x = ".row_id", by.y = "name", all = FALSE)
flip_idx <- g37_lift$strand == "-"
cat(sprintf("strand-flipped rows: %d (alleles will be reverse-complemented)\n", sum(flip_idx)))

rc_alleles <- function(s) {
  out <- s
  ok <- !is.na(s) & s != "-" & s != "" & !grepl("[^ACGTNacgtn]", s)
  if (any(ok)) {
    out[ok] <- as.character(reverseComplement(DNAStringSet(s[ok])))
  }
  out
}
if (any(flip_idx)) {
  g37_lift[flip_idx, referenceAllele := rc_alleles(referenceAllele)]
  g37_lift[flip_idx, variantAllele   := rc_alleles(variantAllele)]
}

g37_lift[, chr           := sub("^chr", "", chrom)]
g37_lift[, startPosition := start + 1L]
g37_lift[, endPosition   := end]
g37_lift[, ncbiBuild     := "GRCh38"]
g37_lift[, lifted        := TRUE]
g37_lift[, c("chrom","start","end","score","strand","chr_ucsc",".row_id") := NULL]

g38[, lifted   := FALSE]
g38[, c("chr_ucsc",".row_id") := NULL]

# ---------- merge & dedupe ----------
# put gdc / native-GRCh38 first so it survives unique() on conflict
ord <- order(factor(g38$source, levels = c("gdc","pan_can_atlas")))
g38 <- g38[ord]

cols <- intersect(names(g37_lift), names(g38))
merged <- rbindlist(list(g38[, ..cols], g37_lift[, ..cols]), use.names = TRUE)

key_cols <- c("sampleId","chr","startPosition","endPosition","referenceAllele","variantAllele")
merged[, .dk := do.call(paste, c(.SD, sep = "|")), .SDcols = key_cols]

src_by_key <- merged[, .(source = paste(sort(unique(source)), collapse = ",")), by = .dk]
merged <- merged[!duplicated(.dk)]
merged[, source := src_by_key$source[match(.dk, src_by_key$.dk)]]
merged[, .dk := NULL]

cat(sprintf("\npost-dedupe rows: %d\n", nrow(merged)))
cat("source breakdown:\n"); print(merged[, .N, by = source])

# ---------- write ----------
fwrite(merged, out_path, sep = "\t")
cat(sprintf("\nWrote %s\n", out_path))

# tidy
file.remove(bed_in); file.remove(bed_out); file.remove(bed_unm)
