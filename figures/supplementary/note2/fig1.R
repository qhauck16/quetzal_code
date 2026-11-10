# Supplementary Note 2, Figure 1 (was S19): TCGA tumours carrying the NKAP neurodevelopmental-
# syndrome substitutions (p.R330C, p.R330H) show the GSP #13 cryptic 3'ss
# phenotype. Data: note2/cache/psi_long.tsv (note2/build_cache.R; ported from
# diff_splicing/nkap_rare_disease.Rmd).
#   a: per-sample ECDFs of PSI over the GSP #13 cryptic 3'ss events (alt 3'ss,
#      -8..+20 nt; NKAP 20v60, dPSI > 0.05, FC > 2.5): every cohort sample in
#      grey, the cohort's per-event median in black, the carrier in red.
#   b: mean PSI over those events, one point per sample; carrier in red with
#      its rank in the cohort.
# Tests (also in figure3/paper_stats.Rmd): empirical rank of the carrier's mean
# PSI within its cohort (P = rank / (n + 1)), and a one-sided paired Wilcoxon
# signed-rank test of carrier PSI vs the cohort median across events.
source("common.R")

long <- fread(file.path("note2", "cache", "psi_long.tsv"), data.table = FALSE,
              colClasses = c(rail_id = "character"))
car  <- fread(file.path("note2", "cache", "carriers.tsv"), data.table = FALSE,
              colClasses = c(rail_id = "character"))
CAR_COL <- "#D62728"
car$plab <- paste0("p.", car$change)                       # p.R330C / p.R330H
panel_lab <- setNames(sprintf("%s (%s)", car$plab, car$cohort), car$cohort)

# ---- statistics ------------------------------------------------------
stats <- lapply(seq_len(nrow(car)), function(i) {
  d   <- long %>% filter(cohort == car$cohort[i])
  med <- d %>% filter(rail_id != car$rail_id[i]) %>% group_by(event) %>%
    summarise(med = median(psi), .groups = "drop")
  cp  <- d %>% filter(rail_id == car$rail_id[i]) %>% select(event, psi) %>%
    inner_join(med, by = "event")
  w   <- wilcox.test(cp$psi, cp$med, paired = TRUE, alternative = "greater")
  sm  <- d %>% group_by(rail_id) %>% summarise(mean_psi = mean(psi), .groups = "drop")
  rk  <- sum(sm$mean_psi >= sm$mean_psi[sm$rail_id == car$rail_id[i]])
  tibble(cohort = car$cohort[i], change = car$change[i], n_events = nrow(cp),
         n_samples = nrow(sm), rank = rk, p_rank = rk / (nrow(sm) + 1),
         p_wilcox = w$p.value,
         carrier_mean = sm$mean_psi[sm$rail_id == car$rail_id[i]],
         cohort_median_mean = median(sm$mean_psi[sm$rail_id != car$rail_id[i]]))
}) %>% bind_rows()
print(as.data.frame(stats), digits = 3)
fwrite(stats, file.path("note2", "cache", "stats.tsv"), sep = "\t")

# ---- a: ECDFs --------------------------------------------------------
cohort_med <- long %>% left_join(car %>% select(cohort, crid = rail_id), by = "cohort") %>%
  filter(rail_id != crid) %>% group_by(cohort, event) %>%
  summarise(psi = median(psi), .groups = "drop")
bg <- long %>% filter(!rail_id %in% car$rail_id)
fg <- long %>% filter(rail_id %in% car$rail_id)
# Panel a carries only the p-value and the legend (Quinn 2026-10-27); the test
# (one-sided paired Wilcoxon, carrier vs cohort median) is in the caption.
ann_a <- stats %>% mutate(lab = sprintf("P = %.1e", p_wilcox))
# Legend ON the plot (never as a subtitle under the title; Quinn 2026-10-27)
KEY <- c(car = "p.R330H / p.R330C tumour", med = "Cohort median per event", oth = "Other tumours in the cohort")
p_a <- ggplot() +
  stat_ecdf(data = bg, aes(psi, group = rail_id, colour = "oth"), linewidth = 0.15, alpha = 0.6) +
  stat_ecdf(data = cohort_med, aes(psi, colour = "med"), linewidth = 0.6) +
  stat_ecdf(data = fg, aes(psi, colour = "car"), linewidth = 0.9) +
  scale_colour_manual(values = c(car = CAR_COL, med = "black", oth = "grey75"), breaks = names(KEY),
                      labels = KEY, name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.9, alpha = 1))) +
  geom_text(data = ann_a, aes(x = 0.4, y = 0.05, label = lab), hjust = 1, vjust = 0,
            size = SIZE_BASE / .pt, family = FONT_FAMILY) +
  facet_wrap(~ cohort, labeller = as_labeller(panel_lab)) +
  coord_cartesian(xlim = c(0, 0.4)) +      # a view: scale limits would drop data
  labs(x = "PSI of GSP #13 cryptic 3'ss events", y = "Cumulative fraction of events") +
  theme_publication() +
  theme(strip.background = element_blank(), strip.text = element_text(size = SIZE_STRIP),
        plot.title = md_title(),
        # middle-right of the right facet is empty (curves have saturated there)
        legend.position = "inside", legend.position.inside = c(0.98, 0.5),
        legend.justification = c(1, 0.5), legend.background = element_blank(),
        legend.text = element_text(size = SIZE_LEGEND, family = FONT_FAMILY),
        legend.key.width = unit(14, "pt"), legend.key.height = unit(9, "pt"))

# ---- b: per-sample mean PSI -----------------------------------------
sm <- long %>% group_by(cohort, rail_id) %>% summarise(mean_psi = mean(psi), .groups = "drop") %>%
  mutate(carrier = rail_id %in% car$rail_id)
ann_b <- sm %>% filter(carrier) %>% left_join(car %>% select(cohort, plab), by = "cohort") %>%
  mutate(lab = plab)                                       # just the substitution
p_b <- ggplot(sm, aes(cohort, mean_psi)) +
  geom_jitter(data = filter(sm, !carrier), width = 0.18, height = 0, size = 0.7,
              colour = "grey55", alpha = 0.7) +
  geom_point(data = filter(sm, carrier), colour = CAR_COL, size = 0.7) +   # same size as the others
  geom_text(data = ann_b, aes(label = lab), hjust = -0.25, vjust = 0.5, colour = CAR_COL,
            size = SIZE_BASE / .pt, family = FONT_FAMILY) +
  scale_x_discrete(labels = function(x) sub("TCGA-", "TCGA-", x)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Mean PSI over cryptic 3'ss events", title = "Carrier rank within its cohort") +
  theme_publication() +
  theme(plot.title = element_text(size = SIZE_TITLE, face = "bold", hjust = 0.5),
        plot.margin = margin(4, 70, 4, 4, "pt"))

fig <- p_a / (p_b | plot_spacer()) + plot_layout(heights = c(1.15, 1)) + plot_annotation_publication() &
  theme(plot.tag = element_text(size = SIZE_PANEL_LABEL, face = "bold", family = FONT_FAMILY))
save_supp(fig, "note2/SN2_Fig1", 8.5, 6.4)
