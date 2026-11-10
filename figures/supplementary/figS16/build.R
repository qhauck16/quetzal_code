# YLPM1 splice-acceptor variant in TCGA-CJ-4638 (KIRC; rank 2 on GSP #30),
# chr14:74812626 A>G = the A of MANE exon 11's AG acceptor. One figure in the
# style of kinnex/fig_zoom/MTIF2/fig_MTIF2_slide.pdf (Quinn, 2026-10-26): no
# titles or panel letters, track labels on the left, grey tracks, gene model
# with chevron introns, red only for the variant, and a callout zoom from the
# exon-11 acceptor to the reference / variant sequence. Top to bottom:
#   RNA-seq coverage with junction arcs (all the same height; width = reads)
#   RNA-seq reads, squished (variant base = red tick), genomic axis
#   YLPM1 MANE exons 9-12 + callout to the intron | exon 11 sequence
# Reads: GDC slices from fetch.sh. The five-panel version with the tumour WGS
# reads is build_full.R. Run from supplementary/:
#   Rscript --vanilla figS16/build.R
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
source("common.R")
suppressPackageStartupMessages({ library(GenomicRanges); library(GenomicAlignments); library(Rsamtools) })
# Bioconductor masks these dplyr verbs (S4Vectors / BiocGenerics)
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; count <- dplyr::count
first <- dplyr::first; desc <- dplyr::desc
D    <- "figS16"
SITE <- 74812626; CHR <- "chr14"
WIN  <- c(74810150, 74816420)                 # MANE exons 9-12
TX   <- "ENST00000325680.12"                  # YLPM1 MANE Select
FA   <- FaFile(file.path(QROOT, "inputs/GRCh38.primary_assembly.genome.fa"))
INK  <- "grey25"; TRACK <- "grey45"; LIGHT <- "grey65"; VAR <- "#D62728"; REF <- "#2CA02C"
TXT  <- SIZE_BASE * 0.9 / .pt

gr  <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
exd <- as.data.frame(sort(gr[gr$type == "exon" & gr$transcript_id == TX]))[, c("start", "end", "exon_number")] %>%
  mutate(exon_number = as.integer(exon_number)) %>% filter(end >= WIN[1], start <= WIN[2]) %>%
  mutate(start = pmax(start, WIN[1]), end = pmin(end, WIN[2]))

pack <- function(s, e, gap = 1) {             # greedy row packing (IGV "squished")
  o <- order(s); row_end <- numeric(0); r <- integer(length(s))
  for (i in o) {
    k <- which(row_end < s[i] - gap)[1]
    if (is.na(k)) { row_end <- c(row_end, e[i]); k <- length(row_end) } else row_end[k] <- e[i]
    r[i] <- k
  }
  r
}

# ---- reads -------------------------------------------------------------------
ra <- readGAlignments(file.path(D, "bams/rna_YLPM1.bam"),
                      param = ScanBamParam(which = GRanges(CHR, IRanges(WIN[1], WIN[2])), what = "seq",
                                           mapqFilter = 20))
ra  <- ra[start(ra) <= WIN[2] & end(ra) >= WIN[1]]
cov <- as.numeric(coverage(ra)[[CHR]][WIN[1]:WIN[2]])
covd <- tibble(pos = WIN[1]:WIN[2], depth = cov)
jd <- as.data.frame(unlist(junctions(ra))) %>% count(start, end, name = "n") %>%
  filter(start >= WIN[1], end <= WIN[2])
blk <- grglist(ra); rid <- pack(start(ra), end(ra))
rb <- as.data.frame(blk) %>% transmute(read = group, start, end) %>% mutate(row = rid[read])
rl <- tibble(start = start(ra), end = end(ra), row = rid)
site_gr <- GRanges(CHR, IRanges(SITE, SITE))
on_site <- which(overlapsAny(blk, site_gr))
rsite <- tibble(read = on_site, base = as.character(stackStringsFromGAlignments(ra[on_site], site_gr))) %>%
  mutate(row = rid[read])
idep <- function(a, b) mean(cov[(a:b) - WIN[1] + 1])
cat(sprintf("[ylpm1] %d reads, %d rows; site bases: %s; intron 10 depth %.1f vs intron 11 %.1f\n",
            length(ra), max(rid), paste(names(table(rsite$base)), table(rsite$base), collapse = " "),
            idep(74811739, 74812627), idep(74812783, 74816202)))
print(jd)

lab_y <- function() theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = SIZE_BASE,
                                                      colour = INK, family = FONT_FAMILY))
no_x  <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
               axis.ticks.x = element_blank(), axis.line.x = element_blank())

# ---- 1: coverage + arcs (every arc the same height; width ∝ reads) ----------
ymax <- max(cov); arc_base <- ymax * 1.05; arc_top <- ymax * 1.45
arcs <- jd %>% rowwise() %>%
  reframe(start, end, n, x = seq(start, end, length.out = 80), t = seq(0, 1, length.out = 80)) %>%
  mutate(y = arc_base + (arc_top - arc_base) * 4 * t * (1 - t), id = paste(start, end))
p1 <- ggplot() +
  geom_col(data = covd, aes(pos, depth), width = 1, fill = TRACK, colour = NA) +
  geom_path(data = arcs, aes(x, y, group = id, linewidth = n), colour = INK, lineend = "round") +
  geom_text(data = jd, aes((start + end) / 2, arc_top + ymax * 0.04, label = n), vjust = 0,
            size = TXT, family = FONT_FAMILY, colour = INK) +
  scale_linewidth(range = c(0.25, 1), guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(breaks = c(0, round(ymax))) +
  coord_cartesian(xlim = WIN, ylim = c(0, arc_top + ymax * 0.22), expand = FALSE) +
  labs(y = "RNA-seq\ncoverage") +
  theme_publication() + no_x + lab_y()

# ---- 2: squished reads + genomic axis ----------------------------------------
p2 <- ggplot() +
  geom_segment(data = rl, aes(x = start, xend = end, y = row, yend = row), colour = LIGHT, linewidth = 0.12) +
  geom_rect(data = rb, aes(xmin = start, xmax = end + 1, ymin = row - 0.47, ymax = row + 0.47),
            fill = TRACK, colour = NA) +
  geom_rect(data = filter(rsite, base != "A"), aes(xmin = SITE - 8, xmax = SITE + 9,
                                                   ymin = row - 0.48, ymax = row + 0.48), fill = VAR) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_reverse() +
  coord_cartesian(xlim = WIN, expand = FALSE) +
  labs(x = sprintf("%s  (+ strand)", CHR), y = "RNA-seq\nreads") +
  theme_publication() + lab_y() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.title.x = element_text(hjust = 0, size = SIZE_BASE, colour = INK))

# ---- 3: gene model + callout to the sequence ----------------------------------
# chevrons every ~250 bp on the intron line, skipped inside exons
chev <- tibble(x = seq(WIN[1] + 120, WIN[2] - 120, by = 250)) %>%
  filter(!vapply(x, function(v) any(v >= exd$start - 40 & v <= exd$end + 40), logical(1)))
seq_ref <- strsplit(as.character(getSeq(FA, GRanges(CHR, IRanges(SITE - 4, SITE + 5)))), "")[[1]]
pos_seq <- (SITE - 4):(SITE + 5)
is_exon <- pos_seq >= 74812628
show <- function(b, ex) ifelse(ex, toupper(b), tolower(b))
seq_var <- seq_ref; seq_var[pos_seq == SITE] <- "G"
# the letters are spread over the window in data units, so the zoom lives in
# the same panel as the gene model and the callout lines are plain segments
mid <- mean(WIN); step <- diff(WIN) / 30; gap <- step * 0.5
xs <- mid + (seq_along(pos_seq) - 5.5) * step + ifelse(is_exon, gap / 2, -gap / 2)
bar_x <- (max(xs[!is_exon]) + min(xs[is_exon])) / 2
Y_GENE <- 0; Y_HEAD <- -1.55; Y_REF <- -2.25; Y_VAR <- -2.95
sq <- bind_rows(tibble(x = xs, y = Y_REF, b = show(seq_ref, is_exon), col = "ref"),
                tibble(x = xs, y = Y_VAR, b = show(seq_var, is_exon),
                       col = ifelse(pos_seq == SITE, "var", "ref")))
box <- c(SITE - 70, SITE + 70)
p3 <- ggplot() +
  annotate("segment", x = WIN[1], xend = WIN[2], y = Y_GENE, yend = Y_GENE, colour = LIGHT, linewidth = 0.5) +
  geom_text(data = chev, aes(x, Y_GENE, label = ">"), colour = LIGHT, size = TXT * 1.1, family = FONT_FAMILY) +
  geom_rect(data = exd, aes(xmin = start, xmax = end, ymin = Y_GENE - 0.32, ymax = Y_GENE + 0.32), fill = INK) +
  annotate("point", x = SITE, y = Y_GENE + 0.55, shape = 25, size = 1.6, fill = VAR, colour = VAR) +
  annotate("text", x = WIN[2], y = Y_GENE + 0.95, hjust = 1, vjust = 0, colour = VAR, size = TXT,
           family = FONT_FAMILY, label = sprintf("▼ variant (%s:%s A>G)", CHR, scales::comma(SITE))) +
  annotate("rect", xmin = box[1], xmax = box[2], ymin = Y_GENE - 0.45, ymax = Y_GENE + 0.45,
           fill = NA, colour = "black", linewidth = 0.3) +
  annotate("segment", x = box, xend = c(min(xs) - step, max(xs) + step), y = Y_GENE - 0.45, yend = Y_HEAD + 0.35,
           linewidth = 0.3) +
  annotate("text", x = mean(xs[!is_exon]), y = Y_HEAD, label = "intron 10", size = TXT, family = FONT_FAMILY,
           colour = "grey40") +
  annotate("text", x = mean(xs[is_exon]), y = Y_HEAD, label = "exon 11", size = TXT, family = FONT_FAMILY,
           colour = "grey40") +
  annotate("segment", x = bar_x, xend = bar_x, y = Y_REF + 0.35, yend = Y_VAR - 0.35, linewidth = 0.3) +
  geom_text(data = sq, aes(x, y, label = b, colour = col), size = TXT * 1.5, family = FONT_FAMILY) +
  annotate("text", x = min(xs) - step * 1.2, y = c(Y_REF, Y_VAR), hjust = 1, size = TXT, family = FONT_FAMILY,
           colour = c(REF, VAR), label = c("reference", "variant")) +
  scale_colour_manual(values = c(ref = "black", var = VAR), guide = "none") +
  scale_x_continuous(labels = scales::comma) +
  coord_cartesian(xlim = WIN, ylim = c(Y_VAR - 0.45, Y_GENE + 1.3), expand = FALSE, clip = "off") +
  labs(y = "YLPM1") +
  theme_publication() + no_x +
  theme(axis.title.y = element_text(face = "italic", angle = 0, vjust = 0.82, hjust = 1, size = SIZE_BASE,
                                    family = FONT_FAMILY, colour = "black"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

fig <- p1 / p2 / p3 + plot_layout(heights = c(1, 0.85, 1.05))   # reads squished tighter

ggsave_pub(file.path(OUT_DIR, "figS16.pdf"), fig, width = 8, height = 5.2)
cat("[ylpm1] wrote", file.path(OUT_DIR, "figS16.pdf"), "\n")
