#!/usr/bin/env Rscript
# Write the DE-junction classification table for a supplementary_plots
# gene + factor (mirrors junction_classifications() in
# ../scripts/lib.R). If config[group][gene]$factor is empty, writes a
# header-only file so snakemake still has the output, and exits cleanly.

script_dir <- dirname(sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))]))
source(file.path(script_dir, "lib_supp.R"))

args        <- commandArgs(trailingOnly = TRUE)
gene        <- args[1]
config_path <- args[2]
out_tsv     <- args[3]

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

empty_header <- function() {
  write.table(data.frame(junction = character(0), Gene_name = character(0),
                          Intron_coord = character(0), Strand = character(0),
                          Annot = character(0), Coding = character(0),
                          UTR = character(0), GencodePC = character(0)),
              out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
}

cfg <- yaml::read_yaml(config_path)
le  <- lookup_entry(cfg, gene)
ent <- le$entry

chr  <- ent$chr
fact <- if (is.null(ent$factor)) "" else as.character(ent$factor)
if (!nzchar(fact)) {
  cat("[de_junctions] factor unset for ", gene, "; writing header-only TSV\n",
      sep = "")
  empty_header(); quit(save = "no", status = 0)
}

res_path    <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "res.RDS")
de_res_path <- file.path(cfg$paths$tcga_chr_dir, chr, gene, "de_res.RDS")
res    <- readRDS(res_path)
de_res <- readRDS(de_res_path)

juncs <- de_junction_names(res, de_res, fact, cfg$thresholds)
cat(sprintf("[de_junctions] %s %s -> %d DE junctions\n",
            gene, fact, length(juncs)))

if (length(juncs) == 0L) { empty_header(); quit(save = "no", status = 0) }

leaf2_class <- fread(cfg$paths$leaf2_class, data.table = FALSE)
leaf2_class$coord <- sub('^.*:', '', leaf2_class$Intron_coord)

jc <- data.frame(junction = juncs) %>%
  mutate(start = str_split_i(junction, '-', 1),
         end   = str_split_i(junction, '-', 2)) %>%
  mutate(start = as.numeric(start) - 1) %>%
  mutate(junction = paste0(start, '-', end)) %>%
  left_join(leaf2_class, by = c("junction" = "coord")) %>%
  dplyr::select(junction, Gene_name, Intron_coord, Strand, Annot,
                 Coding, UTR, GencodePC)

write.table(jc, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[de_junctions] wrote ", out_tsv, "\n", sep = "")
