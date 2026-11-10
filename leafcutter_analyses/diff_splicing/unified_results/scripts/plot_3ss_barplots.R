
suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [df_final, gene, tag,
#        wide_xmin, wide_xmax, zoom_xmin, zoom_xmax,
#        wide_output, zoom_output]
df_path    <- args[1]
gene       <- toupper(args[2])
tag        <- args[3]
wide_xmin  <- as.numeric(args[4])
wide_xmax  <- as.numeric(args[5])
zoom_xmin  <- as.numeric(args[6])
zoom_xmax  <- as.numeric(args[7])
wide_out   <- args[8]
zoom_out   <- args[9]

df <- read_tsv(df_path, show_col_types = FALSE)

plot_hist <- function(d, x_min, x_max, title) {
  if (nrow(d) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste0(title, " (no data)")))
  }
  bin_breaks <- seq(x_min - 0.5, x_max + 0.5, by = 1)
  ggplot(d, aes(x = pos)) +
    geom_histogram(breaks = bin_breaks,
                   color = "white", fill = "steelblue", alpha = 0.85) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    labs(x = "Position of cryptic 3'ss relative to annotated 3'ss",
         y = "Count",
         title = title) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

d_all <- df %>%
  filter(type == "alt 3'ss", !is.na(distance_to_cryptic)) %>%
  mutate(pos = -distance_to_cryptic)

d_wide <- d_all %>% filter(pos >= wide_xmin, pos <= wide_xmax)

p_wide <- plot_hist(d_wide, wide_xmin, wide_xmax,
                    paste0(gene, " ", tag, " - 3'ss position (wide)"))
dir.create(dirname(wide_out), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(wide_out, p_wide, width = 10, height = 4, dpi = 300)

d_zoom <- d_all %>% filter(pos >= zoom_xmin, pos <= zoom_xmax)

zoom_step <- if (zoom_xmax - zoom_xmin <= 30) 1 else 5
p_zoom <- plot_hist(d_zoom, zoom_xmin, zoom_xmax,
                    paste0(gene, " ", tag, " - 3'ss position (zoom)")) +
          scale_x_continuous(breaks = seq(zoom_xmin, zoom_xmax, by = zoom_step))

ggplot2::ggsave(zoom_out, p_zoom, width = 8, height = 4, dpi = 300)
