# YLPM1 splice-acceptor variant in TCGA-CJ-4638 (KIRC; rank 2 on GSP #30):
# chr14:74812626 A>G, the A of MANE exon 11's AG acceptor (intron 10, -2).
# Reads come from GDC BAM slices (fetch.sh). Top to bottom:
#   a  RNA-seq coverage over exons 9-12 with junction arcs (read counts)
#   b  RNA-seq reads, squished (blocks = aligned bases, thin lines = splices)
#   c  MANE exons 9-12 in the same window, variant marked
#   d  tumour WGS reads +-40 bp around the variant, base at the site coloured,
#      reference sequence above
#   e  whole YLPM1 MANE gene, window boxed
# Run from supplementary/:  Rscript --vanilla figS16/build_full.R
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
WIN  <- c(74810150, 74816420)            # exons 9-12 of the MANE transcript
WGSW <- c(SITE - 40, SITE + 40)
FA   <- FaFile(file.path(QROOT, "inputs/GRCh38.primary_assembly.genome.fa"))
TX   <- "ENST00000325680.12"             # YLPM1 MANE Select
BASE_COL <- c(A = "#2E9E44", C = "#1F5FBF", G = "#E08A00", T = "#C0392B")
REF_GREY <- "grey70"

gr <- readRDS(file.path(QROOT, "inputs/gencode_v46_granges.RDS"))
ex <- sort(gr[gr$type == "exon" & gr$transcript_id == TX])
exd <- as.data.frame(ex)[, c("start", "end", "exon_number")] %>% mutate(exon_number = as.integer(exon_number))

# ---- greedy row packing (IGV "squished") -------------------------------------
pack <- function(s, e, gap = 3) {
  o <- order(s); row_end <- numeric(0); r <- integer(length(s))
  for (i in o) {
    k <- which(row_end < s[i] - gap)[1]
    if (is.na(k)) { row_end <- c(row_end, e[i]); k <- length(row_end) } else row_end[k] <- e[i]
    r[i] <- k
  }
  r
}

# ---- RNA ---------------------------------------------------------------------
rw <- GRanges(CHR, IRanges(WIN[1], WIN[2]))
ra <- readGAlignments(file.path(D, "bams/rna_YLPM1.bam"),
                      param = ScanBamParam(which = rw, what = "seq", mapqFilter = 20))
ra <- ra[start(ra) <= WIN[2] & end(ra) >= WIN[1]]
cov <- as.numeric(coverage(ra)[[CHR]][WIN[1]:WIN[2]])
covd <- tibble(pos = WIN[1]:WIN[2], depth = cov)
jn <- unlist(junctions(ra))
jd <- as.data.frame(jn) %>% count(start, end, name = "n") %>%
  filter(start >= WIN[1] - 2000, end <= WIN[2] + 2000) %>%
  mutate(canonical_10_11 = start == 74811739 & end == 74812627)
print(jd)
blk <- grglist(ra)                         # aligned blocks per read (N removed)
rid <- pack(start(ra), end(ra))
rb <- as.data.frame(blk) %>% transmute(read = group, start, end) %>% mutate(row = rid[read])
rl <- tibble(read = seq_along(ra), start = start(ra), end = end(ra), row = rid)
# base at the variant in RNA reads that cover it (not spliced over it)
# stackStrings returns one string per read OVERLAPPING the site, so take the
# reads whose aligned blocks (not splice gaps) cover it first
site_gr <- GRanges(CHR, IRanges(SITE, SITE))
on_site <- which(overlapsAny(blk, site_gr))
rbase <- as.character(stackStringsFromGAlignments(ra[on_site], site_gr))
rsite <- tibble(read = on_site, base = rbase) %>% filter(base %in% names(BASE_COL)) %>%
  mutate(row = rid[read])
cat(sprintf("[ylpm1] RNA reads in window: %d; covering the site: %s\n", length(ra),
            paste(names(table(rsite$base)), table(rsite$base), collapse = " ")))

# intron depths for the caption
idep <- function(a, b) mean(cov[(a:b) - WIN[1] + 1])
dep <- c(exon10 = idep(74811620, 74811738), intron10 = idep(74811739, 74812627),
         exon11 = idep(74812628, 74812782), intron11 = idep(74812783, 74816202))
print(round(dep, 1))

# A VIEW, not scale limits: limits would delete reads / arcs that cross the
# window edges before drawing them
x_win <- list(scale_x_continuous(labels = function(x) scales::comma(x)),
              coord_cartesian(xlim = WIN, expand = FALSE))
theme_track_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(), axis.line.x = element_blank())

# a: coverage + arcs (arc height and width scale with read count)
arc_y0 <- max(cov) * 1.05
jd_plot <- jd %>% filter(start >= WIN[1], end <= WIN[2]) %>%
  mutate(mid = (start + end) / 2, h = arc_y0 + max(cov) * (0.25 + 0.5 * n / max(n)))
arc_pts <- jd_plot %>% rowwise() %>% reframe(
  start, end, n, canonical_10_11,
  x = seq(start, end, length.out = 60), t = seq(0, 1, length.out = 60), h = h) %>%
  mutate(y = arc_y0 + (h - arc_y0) * 4 * t * (1 - t), id = paste(start, end))
p_a <- ggplot() +
  annotate("rect", xmin = 74811739, xmax = 74812627, ymin = -Inf, ymax = Inf, fill = "#FDEBD0", alpha = 0.7) +
  geom_area(data = covd, aes(pos, depth), fill = "grey55", colour = NA) +
  geom_path(data = arc_pts, aes(x, y, group = id, linewidth = n,
                                colour = canonical_10_11)) +
  geom_text(data = jd_plot, aes(mid, h + max(cov) * 0.06, label = n,
                                colour = canonical_10_11), size = SIZE_BASE * 0.85 / .pt,
            family = FONT_FAMILY, vjust = 0) +
  scale_linewidth(range = c(0.3, 1.1), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "#D62728", `FALSE` = "grey25"), guide = "none") +
  scale_x_continuous(labels = function(x) scales::comma(x)) +
  # headroom for the arc counts, which a bare coord_cartesian() would clip
  coord_cartesian(xlim = WIN, ylim = c(0, max(jd_plot$h) + max(cov) * 0.2), expand = FALSE) +
  labs(y = "RNA-seq\ncoverage", title = "TCGA-CJ-4638 (TCGA-KIRC): YLPM1 exon 11 acceptor variant",
       subtitle = sprintf(paste0("Exon 10→11 junction: %d reads (red); intron 10 (shaded) retained: ",
                                 "mean depth %.1f vs %.1f in intron 11; all %d RNA reads at the site carry G"),
                          jd$n[jd$canonical_10_11], dep["intron10"], dep["intron11"], nrow(rsite))) +
  theme_publication() + theme_track_x +
  theme(plot.title = md_title(), plot.subtitle = element_text(size = SIZE_BASE, hjust = 0.5, family = FONT_FAMILY))
p_a$labels$title <- gene_md(p_a$labels$title)

# b: squished RNA reads
p_b <- ggplot() +
  annotate("rect", xmin = 74811739, xmax = 74812627, ymin = -Inf, ymax = Inf, fill = "#FDEBD0", alpha = 0.7) +
  geom_segment(data = rl, aes(x = start, xend = end, y = row, yend = row), colour = "grey75", linewidth = 0.15) +
  geom_rect(data = rb, aes(xmin = start, xmax = end + 1, ymin = row - 0.4, ymax = row + 0.4),
            fill = "grey35", colour = NA) +
  geom_rect(data = rsite, aes(xmin = SITE - 6, xmax = SITE + 7, ymin = row - 0.45, ymax = row + 0.45, fill = base),
            colour = NA) +
  scale_fill_manual(values = BASE_COL, guide = "none") +
  x_win + scale_y_reverse(expand = c(0.01, 0.01)) +
  labs(y = "RNA-seq reads") +
  theme_publication() + theme_track_x +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

# c: MANE exons 9-12 in the window
exw <- exd %>% filter(end >= WIN[1], start <= WIN[2])
p_c <- ggplot() +
  annotate("segment", x = WIN[1], xend = WIN[2], y = 0, yend = 0, linewidth = 0.4) +
  geom_rect(data = exw, aes(xmin = pmax(start, WIN[1]), xmax = pmin(end, WIN[2]), ymin = -0.35, ymax = 0.35),
            fill = "#1F4E79") +
  geom_text(data = exw, aes(x = (pmax(start, WIN[1]) + pmin(end, WIN[2])) / 2, y = -0.75,
                            label = paste("exon", exon_number)), size = SIZE_BASE * 0.85 / .pt, family = FONT_FAMILY) +
  annotate("segment", x = SITE, xend = SITE, y = 0.4, yend = 1.0, colour = "#D62728", linewidth = 0.5) +
  annotate("point", x = SITE, y = 1.05, colour = "#D62728", size = 1.8) +
  annotate("text", x = SITE, y = 1.35, label = "chr14:74,812,626 A>G (acceptor -2)", colour = "#D62728",
           size = SIZE_BASE * 0.85 / .pt, family = FONT_FAMILY, vjust = 0) +
  x_win + scale_y_continuous(limits = c(-1, 1.8)) +
  labs(x = "chr14 position", y = NULL) +
  theme_publication() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

# d: tumour WGS reads +-40 bp, reference above
wa <- readGAlignments(file.path(D, "bams/wgs_tumour_site.bam"),
                      param = ScanBamParam(which = GRanges(CHR, IRanges(WGSW[1], WGSW[2])), what = "seq",
                                           mapqFilter = 20))
wa <- wa[overlapsAny(grglist(wa), GRanges(CHR, IRanges(SITE, SITE)))]   # aligned base at the site
wb <- as.character(stackStringsFromGAlignments(wa, GRanges(CHR, IRanges(SITE, SITE))))
wrow <- pack(pmax(start(wa), WGSW[1]), pmin(end(wa), WGSW[2]), gap = 1)
wd <- tibble(start = pmax(start(wa), WGSW[1]), end = pmin(end(wa), WGSW[2]), row = wrow, base = wb)
tw <- table(factor(wb, levels = c("A", "C", "G", "T")))
cat(sprintf("[ylpm1] WGS reads at site: %s\n", paste(names(tw), tw, collapse = " ")))
ref <- strsplit(as.character(getSeq(FA, GRanges(CHR, IRanges(WGSW[1], WGSW[2]))))[[1]], "")[[1]]
refd <- tibble(pos = WGSW[1]:WGSW[2], base = ref)
p_d <- ggplot() +
  annotate("rect", xmin = SITE + 1.5, xmax = WGSW[2] + 0.5, ymin = -Inf, ymax = 0.2, fill = "#DCE6F2", alpha = 0.6) +
  annotate("rect", xmin = SITE - 0.5, xmax = SITE + 1.5, ymin = -2.1, ymax = -1.0, fill = NA, colour = "black",
           linewidth = 0.3) +                                  # the AG acceptor
  geom_text(data = refd, aes(pos, -1.55, label = base, colour = base), size = SIZE_BASE * 0.8 / .pt,
            family = FONT_FAMILY, fontface = "bold") +
  geom_rect(data = wd, aes(xmin = start - 0.5, xmax = end + 0.5, ymin = row - 0.4, ymax = row + 0.4),
            fill = REF_GREY, colour = NA) +
  # IGV style: only bases that differ from the reference are coloured
  geom_tile(data = filter(wd, base != ref[SITE - WGSW[1] + 1]), aes(x = SITE, y = row, fill = base),
            width = 1, height = 0.8) +
  scale_fill_manual(values = BASE_COL, guide = "none") +
  scale_colour_manual(values = BASE_COL, guide = "none") +
  scale_x_continuous(limits = WGSW + c(-0.5, 0.5), expand = c(0, 0), breaks = SITE + c(-30, -15, 0, 15, 30),
                     labels = scales::comma) +
  scale_y_reverse(expand = c(0.01, 0.01)) +
  labs(x = "chr14 position (±40 bp around the variant)", y = "tumour WGS\nreads",
       subtitle = sprintf(paste0("Intron 10 …CCT[A>G]G | exon 11 (shaded): the AG acceptor (boxed) becomes GG; ",
                                 "tumour WGS %d G / %d A at the site"), tw[["G"]], tw[["A"]])) +
  theme_publication() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.subtitle = element_text(size = SIZE_BASE, hjust = 0.5, family = FONT_FAMILY))

# e: whole gene
gspan <- range(c(exd$start, exd$end))
p_e <- ggplot() +
  annotate("segment", x = gspan[1], xend = gspan[2], y = 0, yend = 0, linewidth = 0.4) +
  geom_rect(data = exd, aes(xmin = start, xmax = end, ymin = -0.35, ymax = 0.35), fill = "#1F4E79") +
  annotate("rect", xmin = WIN[1], xmax = WIN[2], ymin = -0.6, ymax = 0.6, fill = NA, colour = "#D62728",
           linewidth = 0.4) +
  annotate("text", x = gspan[1], y = 0.75, hjust = 0, vjust = 0, label = "YLPM1 (MANE ENST00000325680), + strand",
           size = SIZE_BASE * 0.9 / .pt, family = FONT_FAMILY, fontface = "italic") +
  scale_x_continuous(expand = c(0.01, 0.01), labels = scales::comma) +
  scale_y_continuous(limits = c(-0.8, 1.3)) +
  labs(x = "chr14 position", y = NULL) +
  theme_publication() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

fig <- p_a / p_b / p_c / p_d / p_e +
  plot_layout(heights = c(1.1, 1.6, 0.55, 1.6, 0.45)) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))

ggsave_pub(file.path(OUT_DIR, "figS16_full.pdf"), fig, width = 8.5, height = 10)
cat("[ylpm1] wrote", file.path(OUT_DIR, "figS16_full.pdf"), "\n")
