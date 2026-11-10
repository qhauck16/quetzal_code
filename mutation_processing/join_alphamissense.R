#!/usr/bin/env Rscript
# Aggregate AlphaMissense scores per (chr,pos,ref,alt) (max pathogenicity across
# transcripts) and left-join onto merged_grch38.tsv. Adds AM columns; rows that
# are not missense SNVs or have no AM coverage get NA.

suppressPackageStartupMessages({ library(data.table) })

work <- "."   # run from mutation_processing/
setwd(work)

cat("Reading inputs ...\n")
am <- fread("am_hits.tsv", header = FALSE,
  col.names = c("chr_chr","pos","ref","alt","genome","uniprot_id",
                "transcript_id","protein_variant","am_pathogenicity","am_class"))
am[, chr := sub("^chr","", chr_chr)]
am[, c("chr_chr","genome") := NULL]

cat(sprintf("AM rows matching keys: %d (over %d unique loci)\n",
            nrow(am), uniqueN(am[, .(chr, pos, ref, alt)])))

# Aggregate: keep transcript with max pathogenicity per locus.  Also count and
# carry n_transcripts.
n_tx <- am[, .(am_n_transcripts = .N), by = .(chr, pos, ref, alt)]
setorder(am, chr, pos, ref, alt, -am_pathogenicity)
am_agg <- am[, .SD[1], by = .(chr, pos, ref, alt)]
am_agg <- merge(am_agg, n_tx, by = c("chr","pos","ref","alt"))

cat("\nAM class breakdown (aggregated, all loci):\n")
print(am_agg[, .N, by = am_class])

cat("\nLoading merged file ...\n")
m <- fread("merged_grch38.tsv")

cat("Joining ...\n")
# rename for clarity in output
setnames(am_agg,
  c("uniprot_id","transcript_id","protein_variant","am_pathogenicity","am_class","am_n_transcripts"),
  c("am_uniprot_id","am_transcript_id","am_protein_variant","am_pathogenicity","am_class","am_n_transcripts"))

m_out <- merge(m, am_agg,
  by.x = c("chr","startPosition","referenceAllele","variantAllele"),
  by.y = c("chr","pos","ref","alt"),
  all.x = TRUE, sort = FALSE)

# Coverage stats among missense SNVs
miss_mask <- m_out$mutationType == "Missense_Mutation" &
             nchar(m_out$referenceAllele) == 1 & nchar(m_out$variantAllele) == 1 &
             m_out$referenceAllele %in% c("A","C","G","T") &
             m_out$variantAllele   %in% c("A","C","G","T")
miss_mask[is.na(miss_mask)] <- FALSE
n_miss   <- sum(miss_mask)
n_scored <- sum(miss_mask & !is.na(m_out$am_pathogenicity))
cat(sprintf("\nMissense SNV rows: %d\n", n_miss))
cat(sprintf("Scored by AM:      %d  (%.2f%%)\n", n_scored, 100*n_scored/n_miss))

cat("\nClass breakdown among scored missense rows:\n")
print(m_out[miss_mask & !is.na(am_pathogenicity), .N, by = am_class])

cat("\nPathogenicity distribution among scored missense rows:\n")
print(summary(m_out[miss_mask & !is.na(am_pathogenicity), am_pathogenicity]))

cat("\nWriting merged_grch38_with_am.tsv ...\n")
fwrite(m_out, "merged_grch38_with_am.tsv", sep = "\t")
cat(sprintf("rows: %d   cols: %d\n", nrow(m_out), ncol(m_out)))
cat("Done.\n")
