# How gencode_v46_introns.tsv was built.
#
# Source recipe was originally written and then commented out at the top
# of deprecated/nmd_prevalence.R; this script reproduces it as a
# documentation/regeneration artifact.
#
# Inputs:
#   ../../hg38_granges.RDS         # GENCODE v46 GTF read in as a GRanges
#                                    (see flashier/just_for_fun/... for the
#                                    builder; importable via rtracklayer)
# Outputs:
#   gencode_v46_introns.tsv        # 1-based inclusive intron coordinates,
#                                    one row per (transcript, intron),
#                                    annotated with transcript_id, gene_id,
#                                    gene_name.


suppressPackageStartupMessages({
  library(GenomicFeatures)   # pulls in GenomicRanges + S4Vectors
  library(readr)             # write_tsv only
})

GRANGES_PATH <- "../../inputs/gencode_v46_granges.RDS"
OUT_PATH     <- "gencode_v46_introns.tsv"

# 1. Load the cached GENCODE v46 GRanges.
gencode <- readRDS(GRANGES_PATH)

# 2. Build a TxDb, extract per-transcript introns, flatten to one GRanges.
txdb          <- makeTxDbFromGRanges(gencode)
introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
introns       <- unlist(introns_by_tx, use.names = FALSE)

# 3. Re-attach the parent transcript_id, then join transcript->gene metadata
#    from the original GRanges so downstream code can filter by gene_name.
introns$transcript_id <- names(introns_by_tx)[
  rep(seq_along(introns_by_tx), elementNROWS(introns_by_tx))
]
tx2gene <- unique(
  mcols(gencode)[gencode$type == "transcript",
                 c("transcript_id", "gene_id", "gene_name")]
)
mcols(introns) <- cbind(
  mcols(introns),
  tx2gene[match(introns$transcript_id, tx2gene$transcript_id), ]
)

# 4. Dump as TSV (1-based inclusive coords, the GRanges convention).
intron_df <- as.data.frame(introns)
write_tsv(intron_df, OUT_PATH)

cat(sprintf("Wrote %s: %d rows x %d cols\n",
            OUT_PATH, nrow(intron_df), ncol(intron_df)))
