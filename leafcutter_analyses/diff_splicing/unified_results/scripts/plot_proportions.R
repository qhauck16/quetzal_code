
suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [df_final, gene, tag, output]
df_path  <- args[1]
gene     <- toupper(args[2])
tag      <- args[3]
out_path <- args[4]

df <- read_tsv(df_path, show_col_types = FALSE)

if (nrow(df) == 0) {
  p <- ggplot() + theme_void() +
    labs(title = paste0(gene, " ", tag, " - splice event types (no hits)"))
} else {
  plot_df <- df %>%
    dplyr::count(type) %>%
    mutate(prop = n / sum(n),
           label = scales::percent(prop))

  p <- ggplot(plot_df, aes(x = "", y = prop, fill = type)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_label(aes(label = label),
               position = position_stack(vjust = 0.5),
               label.size = 0,
               label.r = unit(0.15, "lines"),
               label.padding = unit(0.2, "lines")) +
    scale_fill_brewer(palette = "Set2") +
    theme_void(base_size = 14) +
    theme(legend.title = element_blank(), legend.position = "right") +
    labs(title = paste0(gene, " ", tag, " (n = ", nrow(df), ")"))
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(out_path, p, width = 6, height = 4, dpi = 300)
