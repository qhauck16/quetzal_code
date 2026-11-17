# Supplementary Table 3: COSMIC Tier 1 genes with GSP-associated alternative
# splice sites. For each GSP's 20v60 LeafCutter comparison, every junction that
# is upregulated (dPSI > 1%, PSI fold change > 2.5, cluster BH-adjusted P <= 0.05)
# and is an alt 3'ss (GSP #10 / SF3B1, GSP #13 / NKAP) or an alt 3'ss / alt 5'ss
# (GSP #30 / YLPM1), in a COSMIC Tier 1 gene. Splice types and distance_to_cryptic
# are relative to MANE introns (unified_results df_final, dpsi0.01_fc2.5).
# `in_text_count` marks the GSP #13 rows behind "80 COSMIC Tier 1 genes carried
# GSP #13-associated 3'ss displacements": alt 3'ss 1-27 nt downstream of the
# MANE 3'ss (figure3/paper_stats.Rmd).
# Output: figures/supplementary/out/tables/STable3_cosmic_alt_ss.tsv
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

UNI <- file.path(QROOT, "leafcutter_analyses/diff_splicing/unified_results")
SETS <- list(list(gsp = 10, gene = "sf3b1", label = "GSP #10 (SF3B1)", types = "alt 3'ss"),
             list(gsp = 13, gene = "nkap",  label = "GSP #13 (NKAP)",  types = "alt 3'ss"),
             list(gsp = 30, gene = "ylpm1", label = "GSP #30 (YLPM1)", types = c("alt 3'ss", "alt 5'ss")))
cos <- fread(file.path(QROOT, "inputs/cosmic_tier_1.tsv"), data.table = FALSE)
role <- setNames(cos$`Role in Cancer`, cos$`Gene Symbol`)

tab <- bind_rows(lapply(SETS, function(s) {
  fread(file.path(UNI, s$gene, "20v60/dpsi0.01_fc2.5/df_final.tsv"), data.table = FALSE) %>%
    filter(deltapsi > 0.01, test / control > 2.5, type %in% s$types, gene_name %in% names(role)) %>%
    transmute(gene = gene_name, cosmic_role = role[gene_name], gsp = s$label,
              comparison = "20v60", junction = paste0(chrom, ":", start, "-", end), strand,
              splice_type = type,
              mane_intron = paste0(chrom, ":", mane_intron_start, "-", mane_intron_end),
              distance_to_cryptic_nt = distance_to_cryptic,
              psi_control = control, psi_test = test, delta_psi = deltapsi,
              psi_fold_change = test / control, log_effect_size = logef,
              cluster_id, cluster_loglr = loglr, cluster_df = df,
              cluster_p = p, cluster_p_adjusted = p.adjust,
              in_text_count = s$gsp == 13 & type == "alt 3'ss" & !is.na(distance_to_cryptic) &
                distance_to_cryptic >= -27 & distance_to_cryptic <= -1)
})) %>% arrange(gsp, gene, desc(delta_psi))

dir.create(file.path(QROOT, "figures/supplementary/out/tables"), showWarnings = FALSE, recursive = TRUE)
out <- file.path(QROOT, "figures/supplementary/out/tables/STable3_cosmic_alt_ss.tsv")
fwrite(tab, out, sep = "\t", quote = FALSE)
print(tab %>% group_by(gsp) %>% summarise(junctions = n(), genes = n_distinct(gene),
                                          genes_in_text_count = n_distinct(gene[in_text_count])))
cat("wrote", out, "\n")
