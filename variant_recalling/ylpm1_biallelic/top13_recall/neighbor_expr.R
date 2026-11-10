# Possible YLPM1-locus deletions in the top-13 GSP #30 samples: per-cohort z-scored RPM
# (tcga_gene_expression/combined_scaled_rpm.tsv, the Fig. 4g matrix: z within each TCGA project) of the
# genes around YLPM1 + YLPM1 itself (see GENES).
# One panel, all 13 samples; the 4 with putative two-LoF evidence are coloured and labelled.
# Out: neighbor_expr.pdf
source("../../../figures/supplementary/common.R")
# The N_NEAR protein-coding genes nearest YLPM1 by strand-IGNORANT distance (distance() is strand-aware and
# NAs opposite-strand genes -- that bug first dropped LTBP2, AREL1, PROX2, RPS6KL1), + YLPM1, chr14 order.
N_NEAR <- 10
suppressPackageStartupMessages(library(GenomicRanges))
gg <- readRDS("../../../inputs/gencode_v46_granges.RDS")
gg <- gg[gg$type == "gene" & gg$gene_type == "protein_coding" & seqnames(gg) == "chr14"]
gg <- gg[!duplicated(gg$gene_name)]
yl <- gg[gg$gene_name == "YLPM1"]
dd <- distance(gg, yl, ignore.strand = TRUE)
near <- gg[order(dd)][1:(N_NEAR + 1)]
GENES <- near$gene_name[order(start(near))]
cat("[neighbors] genes:", paste(sprintf("%s(%dkb)", GENES, round(dd[match(GENES, gg$gene_name)] / 1e3)), collapse = " "), "\n")
# Three groups (Quinn, 2026-10-30):
#   strong double mutants (2): 461535 (Y1506*/P1507fs + R1719*), 866724 (A252fs + Q437*);
#   likely single mutants: the other 11 of the top 13 on GSP #30;
#   LoF mutants with GSP #30 < 1e-4: variant table, LoF = protein change matching fs|*|splice (Fig. 4a rule).
source("../../../figures/scripts/vt_paths.R")
DOUBLE <- c("461535", "866724")
s <- fread("samples.tsv", colClasses = list(character = "rail_id"))
vt <- fread(vt_path("YLPM1_full_plus_variants"))[!is.na(k30)]
low <- vt[k30 < 1e-4 & grepl("fs|\\*|splice", proteinChange)]
low <- low[, .(rail_id = as.character(rail_id), rank = NA_integer_, project_id)]
ALL <- rbind(s[, .(rail_id, rank, project_id)], low)
cat(sprintf("[neighbors] %d top-13 + %d LoF mutants with GSP #30 < 1e-4\n", nrow(s), nrow(low)))
E <- "../../../expression/combined_scaled_rpm.tsv"
cmd <- sprintf("awk -F'\\t' -v G='%s' -v R='%s' 'BEGIN{n=split(G,g,\",\"); for(i=1;i<=n;i++) w[g[i]]=1; m=split(R,r,\",\"); for(i=1;i<=m;i++) want[r[i]]=1}
  NR==1{for(i=2;i<=NF;i++) if($i in want) c[i]=$i; next} ($1 in w){for(i in c) print $1\"\\t\"c[i]\"\\t\"$i}' %s",
  paste(GENES, collapse = ","), paste(ALL$rail_id, collapse = ","), E)
d <- fread(cmd = cmd, header = FALSE, col.names = c("gene", "rail_id", "z"), colClasses = c("character", "character", "numeric"))
d <- merge(d, ALL, by = "rail_id")
G3 <- c("Strong double mutants (2)", sprintf("Likely single mutants (%d)", nrow(s) - length(DOUBLE)),
        sprintf("LoF mutants, GSP #30 < 1e-4 (%d)", uniqueN(d[is.na(rank)]$rail_id)))
d[, grp := factor(fifelse(rail_id %in% DOUBLE, G3[1], fifelse(!is.na(rank), G3[2], G3[3])), G3)]
d[, gene := factor(gene, GENES)]
print(d[, .(median_z = round(median(z), 2), n = .N), by = .(gene, grp)][order(gene, grp)])
p <- ggplot(d, aes(gene, z, fill = grp)) + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.75, position = position_dodge(0.85)) +
  # points on ALL groups (dodge must see all 3 levels to align), hidden for the big GSP #30 < 1e-4 group
  geom_point(aes(alpha = grp), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.85, seed = 1), size = 0.9) +
  scale_alpha_manual(values = setNames(c(1, 1, 0), G3), guide = "none") +
  coord_cartesian(ylim = c(-2.6, 2.6)) +
  scale_fill_manual(values = setNames(c("#D62728", "#FF9896", "#C7C7C7"), G3), name = NULL) +
  labs(x = "Genes in chromosome order (chr14, around YLPM1)", y = "Expression z-score within cancer type") +
  theme_publication() + angled_x_labels() + theme(legend.position = "inside", legend.position.inside = c(0.78, 0.92))
ggsave_pub("neighbor_expr.pdf", p, width = 10, height = 4.2)
