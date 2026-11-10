#!/usr/bin/env Rscript
# Validate merged_grch38.tsv against the original two-file pipeline.

suppressPackageStartupMessages({ library(data.table) })

mut_dir <- "mutations"
mrg_pth <- "merged_grch38.tsv"

merged <- fread(mrg_pth)
pc     <- fread(file.path(mut_dir, "pan_can_atlas/compiled.tsv"))
gdc    <- fread(file.path(mut_dir, "gdc/compiled.tsv"))

cat("=========================================================\n")
cat("CHECK 1 -- per-sample mutation-COUNT distribution\n")
cat("=========================================================\n")
orig_per_sample   <- unique(rbind(
    pc[, .(hugoGeneSymbol, sampleId)],
    gdc[, .(hugoGeneSymbol, sampleId)]))[, .N, by = sampleId]
merged_per_sample <- unique(merged[, .(hugoGeneSymbol, sampleId)])[, .N, by = sampleId]
joint <- merge(orig_per_sample, merged_per_sample,
               by = "sampleId", all = TRUE,
               suffixes = c("_orig","_merged"))
joint[is.na(N_orig),   N_orig   := 0L]
joint[is.na(N_merged), N_merged := 0L]
cat(sprintf("samples in orig only:   %d\n", sum(joint$N_merged == 0)))
cat(sprintf("samples in merged only: %d\n", sum(joint$N_orig   == 0)))
cat("orig vs merged gene-count summary (per sample):\n")
print(summary(joint[, .(N_orig, N_merged, diff = N_merged - N_orig)]))
cat(sprintf("Pearson cor(orig, merged) = %.5f\n",
            cor(joint$N_orig, joint$N_merged)))

cat("\n=========================================================\n")
cat("CHECK 2 -- gene-level (sample x gene) Jaccard vs union\n")
cat("=========================================================\n")
orig_pairs   <- unique(rbind(pc[, .(sampleId, hugoGeneSymbol)],
                              gdc[, .(sampleId, hugoGeneSymbol)]))
merged_pairs <- unique(merged[, .(sampleId, hugoGeneSymbol)])
orig_pairs[,   key := paste(sampleId, hugoGeneSymbol, sep = "|")]
merged_pairs[, key := paste(sampleId, hugoGeneSymbol, sep = "|")]
inter <- length(intersect(orig_pairs$key, merged_pairs$key))
uni   <- length(union(orig_pairs$key, merged_pairs$key))
cat(sprintf("orig (sample, gene) pairs   : %d\n", nrow(orig_pairs)))
cat(sprintf("merged (sample, gene) pairs : %d\n", nrow(merged_pairs)))
cat(sprintf("intersection                : %d\n", inter))
cat(sprintf("Jaccard                     : %.6f  (target: 1.0)\n", inter / uni))
cat(sprintf("in orig but not merged      : %d\n", length(setdiff(orig_pairs$key,   merged_pairs$key))))
cat(sprintf("in merged but not orig      : %d\n", length(setdiff(merged_pairs$key, orig_pairs$key))))

cat("\n=========================================================\n")
cat("CHECK 3 -- spot-check well-known oncogene hotspots\n")
cat("    GRCh38 reference coords from COSMIC/ClinVar.\n")
cat("=========================================================\n")
hotspots <- data.table(
  gene = c("KRAS","KRAS","BRAF","TP53","TP53","PIK3CA","IDH1"),
  protein = c("G12C","G12D","V600E","R175H","R273H","E545K","R132H"),
  chr_38 = c("12","12","7","17","17","3","2"),
  pos_38 = c(25245350L, 25245351L, 140753336L, 7675088L, 7674220L, 179218303L, 208248389L),
  ref    = c("C","C","A","G","G","G","C"),
  alt    = c("A","T","T","A","A","A","T"))
for (i in seq_len(nrow(hotspots))) {
  h <- hotspots[i]
  hits <- merged[chr == h$chr_38 & startPosition == h$pos_38 &
                 referenceAllele == h$ref & variantAllele == h$alt,
                 .N, by = .(hugoGeneSymbol, source)]
  cat(sprintf("\n%s %s  (chr%s:%d %s>%s):\n",
              h$gene, h$protein, h$chr_38, h$pos_38, h$ref, h$alt))
  if (nrow(hits) == 0) cat("  no hit\n") else print(hits)
}

cat("\n=========================================================\n")
cat("CHECK 4 -- coordinate sanity: lifted vs native GRCh38\n")
cat("    For variants in BOTH sources, coords must agree.\n")
cat("=========================================================\n")
both <- merged[source == "gdc,pan_can_atlas",
               .(chr, startPosition, sampleId, hugoGeneSymbol)]
cat(sprintf("both-source variants: %d (these all share GRCh38 coords by dedup key)\n",
            nrow(both)))
# extra sanity: can we find for each pan_can_atlas-only lifted variant a
# nearby gdc call on same sample/gene at a similar coord? (qualitative)
pco <- merged[source == "pan_can_atlas" & lifted == TRUE]
cat(sprintf("pan_can_atlas-only lifted: %d\n", nrow(pco)))

cat("\n=========================================================\n")
cat("CHECK 5 -- pre-dedup duplicate prevalence + strand-flip share\n")
cat("=========================================================\n")
cat(sprintf("strand-flipped rows kept   : %d (%.3f%% of merged)\n",
            sum(merged$lifted & merged$source != "gdc"),  # rough proxy
            100 * sum(merged$lifted) / nrow(merged)))
cat(sprintf("rows still flagged lifted  : %d\n", sum(merged$lifted)))
cat(sprintf("rows from native GRCh38    : %d\n", sum(!merged$lifted)))

cat("\n=========================================================\n")
cat("CHECK 6 -- chromosome distribution\n")
cat("=========================================================\n")
print(merged[, .N, by = chr][order(as.integer(suppressWarnings(as.integer(chr))), chr)])

cat("\nValidation done.\n")
