# Allele-specific copy number at YLPM1 from PanCanAtlas ABSOLUTE (tcga/pancanatlas/, downloaded 2026-10-30
# from GDC PanCanAtlas: TCGA_mastercalls.abs_segtabs.fixed.txt + abs_tables purity/ploidy). Segments are
# hg19: YLPM1 = chr14:75,230,069-75,326,186 (hg19), midpoint used for the overlap.
# Groups: strong double mutants (461535, 866724), likely single mutants (other 11 of top 13 on GSP #30),
# LoF mutants with GSP #30 < 1e-4 (lof_low_gsp30/samples.tsv). Sample = tumour barcode's first 15 chars.
# Out: ylpm1_absolute_cn.tsv
suppressPackageStartupMessages(library(data.table))
PA <- "../../inputs/pancanatlas"   # run from variant_recalling/ylpm1_biallelic/
seg <- fread(file.path(PA, "TCGA_mastercalls.abs_segtabs.fixed.txt"))[Chromosome == 14 & Start <= 75278000 & End >= 75278000]
pur <- fread(file.path(PA, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"))
top <- fread("top13_recall/samples.tsv")[, grp := fifelse(rail_id %in% c(461535, 866724), "double", "likely_single")]
low <- fread("lof_low_gsp30/samples.tsv")[, grp := "low_GSP30_LoF"]
s <- rbind(top[, .(rail_id, project_id, rna_barcode, proteinChange, grp)], low[, .(rail_id, project_id, rna_barcode, proteinChange, grp)])
s[, array := substr(rna_barcode, 1, 15)]
x <- merge(s, seg[, .(array = Sample, HSCN1 = Modal_HSCN_1, HSCN2 = Modal_HSCN_2, totalCN = Modal_Total_CN, LOH, homdel = Homozygous_deletion)],
           by = "array", all.x = TRUE)
x <- merge(x, pur[, .(array, purity, ploidy)], by = "array", all.x = TRUE)
x[, minor := pmin(HSCN1, HSCN2)]
fwrite(x[order(grp, -LOH)], "ylpm1_absolute_cn.tsv", sep = "\t")
print(x[, .(n = .N, with_CN = sum(!is.na(LOH)), LOH = sum(LOH == 1, na.rm = TRUE), homdel = sum(homdel == 1, na.rm = TRUE),
            pct_LOH = round(100 * mean(LOH == 1, na.rm = TRUE)), median_purity = median(purity, na.rm = TRUE)), by = grp])
options(width = 200); print(x[grp != "low_GSP30_LoF", .(grp, rail_id, project_id, proteinChange, HSCN1, HSCN2, totalCN, LOH, purity, ploidy)][order(grp)])
