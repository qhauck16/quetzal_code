# Supplementary Table 1: cancer-type-specific junction groups (gene-level Quetzal).
# One row per elevated junction group x TCGA cohort: a gene-level Poisson NMF
# junction group (factor) whose loading is elevated in that cohort (linear model,
# Bonferroni) and that has >= 1 unique junction with Bonferroni-significant higher
# PSI in the cohort than in matched GTEx tissue (beta-binomial GLM, beta > 0).
# The text's count of junction groups (e.g. 4,279) is the number of distinct
# gene x factor pairs here; a group elevated in two cohorts has two rows.
# Inputs: gene_level_analysis/junction_level_sig_upregulated.tsv (one row per
# significant junction x hit; gene_cancer_type_hits.Rmd), inputs/cosmic_tier_1.tsv.
# Output: figures/supplementary/out/tables/STable1_cancer_type_splicing.tsv
#   Rscript --vanilla build_table_1_cancer_type_splicing.R   (from supplementary/)
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(data.table); library(dplyr) })

j <- fread(file.path(QROOT, "gene_level_analysis/junction_level_sig_upregulated.tsv"), data.table = FALSE)
cos <- fread(file.path(QROOT, "inputs/cosmic_tier_1.tsv"), data.table = FALSE)
role <- setNames(cos$`Role in Cancer`, cos$`Gene Symbol`)

fmt <- function(x, d = 3) trimws(formatC(x, digits = d, format = "g"))
tab <- j %>%
  arrange(gene_name, factor, cancer_type, bb_pvalue) %>%
  mutate(junction = paste0(chr, ":", junction_id)) %>%
  group_by(gene = gene_name, chromosome = chr, junction_group = factor,
           tcga_cohort = cancer_type, gtex_tissue = gtex_study) %>%
  summarise(
    loading_coefficient       = dplyr::first(factor_cancer_beta),
    loading_p                 = dplyr::first(factor_cancer_pvalue),
    n_significant_junctions   = n(),
    junctions                 = paste(junction, collapse = ","),
    junction_bb_beta          = paste(fmt(bb_beta), collapse = ","),
    junction_bb_p             = paste(fmt(bb_pvalue), collapse = ","),
    junction_bb_p_bonferroni  = paste(fmt(bb_pvalue_bonf), collapse = ","),
    junction_median_psi_tcga  = paste(fmt(median_psi_tcga), collapse = ","),
    junction_median_psi_gtex  = paste(fmt(median_psi_gtex), collapse = ","),
    .groups = "drop") %>%
  mutate(notes = ifelse(gene %in% names(role),
                        paste0("COSMIC Tier 1 (", role[gene], ")"), "")) %>%
  arrange(gene, junction_group, tcga_cohort)

dir.create(file.path(QROOT, "figures/supplementary/out/tables"), showWarnings = FALSE, recursive = TRUE)
out <- file.path(QROOT, "figures/supplementary/out/tables/STable1_cancer_type_splicing.tsv")
fwrite(tab, out, sep = "\t", quote = FALSE)
cat(sprintf("Table 1: %d rows | %d junction groups (gene x factor) | %d genes | %d COSMIC Tier 1 genes\n",
            nrow(tab), n_distinct(paste(tab$gene, tab$junction_group)), n_distinct(tab$gene),
            n_distinct(tab$gene[tab$notes != ""])))
cat("wrote", out, "\n")
