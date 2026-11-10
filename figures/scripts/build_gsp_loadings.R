# Cache the GSP loadings every variant table keys off: flashier::ldf(out,
# type = "i")$L of the flashier fit (not out$L_pm: ldf rescales).
# Heavy (3.7 GB RDS) -- run on SLURM from quetzal_code/figures/:
#   sbatch -p caslake --account=pi-yangili1 --mem=48G -t 1:00:00 --wrap \
#     "Rscript --vanilla scripts/build_gsp_loadings.R"
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({ library(yaml); library(data.table) })
out <- readRDS(file.path(QROOT, "quetzal_runs/genome_wide/softimpute_flash_300_qc_filtered.RDS"))
L   <- flashier::ldf(out, type = "i")$L
d   <- data.table(rail_id = as.numeric(rownames(L)), L)
setnames(d, c("rail_id", paste0("k", seq_len(ncol(L)))))
fwrite(d, "shared_cache/gsp_loadings.tsv", sep = "\t")
cat("wrote shared_cache/gsp_loadings.tsv:", nrow(d), "samples x", ncol(L), "GSPs\n")
