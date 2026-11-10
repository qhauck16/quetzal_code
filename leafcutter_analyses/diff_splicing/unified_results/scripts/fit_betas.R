
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
# args: [feats, gene, tag, lim_low, lim_high, no_bp_out, with_bp_out]
feats_path  <- args[1]
gene        <- toupper(args[2])
tag         <- args[3]
lim_low     <- as.numeric(args[4])
lim_high    <- as.numeric(args[5])
no_bp_out   <- args[6]
with_bp_out <- args[7]

feats <- read_tsv(feats_path, show_col_types = FALSE)

empty_plot <- function(title) {
  ggplot() + theme_void() +
    labs(title = paste0(title, " (insufficient data)"))
}

beta_tile <- function(fit, title) {
  s <- summary(fit)$coefficients
  df <- as.data.frame(s) %>%
    rownames_to_column("variable") %>%
    rename(estimate = Estimate,
           std_error = `Std. Error`,
           z_value  = `z value`,
           p_value  = `Pr(>|z|)`) %>%
    filter(variable != "(Intercept)") %>%
    mutate(variable = gsub("scale\\((.*)\\)", "\\1", variable),
           sig = case_when(
             p_value < 0.001 ~ "***",
             p_value < 0.01  ~ "**",
             p_value < 0.05  ~ "*",
             TRUE            ~ "ns"),
           label = paste0(round(estimate, 3), " ", sig)) %>%
    arrange(variable)

  ggplot(df, aes(x = variable, y = 1, fill = estimate)) +
    geom_tile(color = "white", height = 0.6) +
    geom_text(aes(y = 1.05, label = label), size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(lim_low, lim_high),
                         name = "Beta") +
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right") +
    ggtitle(title)
}

safe_fit <- function(data, formula) {
  if (nrow(data) < 10 || length(unique(data$condition)) < 2) return(NULL)
  data$condition <- factor(data$condition, levels = c("background", "test"))
  tryCatch(glm(formula, family = binomial, data = data),
           error = function(e) NULL)
}

dir.create(dirname(no_bp_out),   recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(with_bp_out), recursive = TRUE, showWarnings = FALSE)

if (nrow(feats) == 0) {
  ggplot2::ggsave(no_bp_out,
                  empty_plot(paste0(gene, " ", tag, " (no branchpoint)")),
                  width = 5, height = 3, dpi = 300)
  ggplot2::ggsave(with_bp_out,
                  empty_plot(paste0(gene, " ", tag, " (with branchpoint)")),
                  width = 5, height = 3, dpi = 300)
  quit(save = "no", status = 0)
}

# ---- Model 1: without branchpoint ----
d1 <- feats %>%
  filter(!is.na(gc_content), !is.na(length),
         !is.na(score_3ss), !is.na(score_5ss),
         !is.na(downstream_AG_distance), !is.na(upstream_AG_distance))

fit1 <- safe_fit(
  d1,
  condition ~ scale(gc_content) + scale(length) +
              scale(score_3ss)  + scale(score_5ss) +
              scale(downstream_AG_distance) +
              scale(upstream_AG_distance)
)

title1 <- paste0(gene, " ", tag, " (no branchpoint)\n",
                 "test: ", sum(d1$condition == "test"),
                 ", background: ", sum(d1$condition == "background"))

p1 <- if (is.null(fit1)) empty_plot(title1) else beta_tile(fit1, title1)
ggplot2::ggsave(no_bp_out, p1, width = 5, height = 3, dpi = 300)

# ---- Model 2: with branchpoint ----
d2 <- d1 %>%
  filter(!is.na(Branchpoint_score), !is.na(bp_distance_to_3ss),
         bp_distance_to_3ss < 100) %>%
  group_by(subclass) %>%
  filter(n() == 2) %>%
  ungroup()

fit2 <- safe_fit(
  d2,
  condition ~ scale(gc_content) + scale(length) +
              scale(score_3ss)  + scale(score_5ss) +
              scale(downstream_AG_distance) +
              scale(upstream_AG_distance) +
              scale(Branchpoint_score) +
              scale(bp_distance_to_3ss)
)

title2 <- paste0(gene, " ", tag, " (with branchpoint)\n",
                 "test: ", sum(d2$condition == "test"),
                 ", background: ", sum(d2$condition == "background"))

p2 <- if (is.null(fit2)) empty_plot(title2) else beta_tile(fit2, title2)
ggplot2::ggsave(with_bp_out, p2, width = 5, height = 3, dpi = 300)
