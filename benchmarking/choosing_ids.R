library(tidyverse)

tcga_metadata <- fread('../inputs/tcga_v2_samples.tsv', data.table = F)
cancer_types <- c('TCGA-LUAD', 'TCGA-TGCT', 'TCGA-KICH')
set.seed(123)  # for reproducibility

tcga_subset <- tcga_metadata %>%
  filter(
    cgc_sample_sample_type == "Primary Tumor",
    gdc_cases.project.project_id %in% cancer_types
  ) %>%
  group_by(gdc_cases.project.project_id) %>%
  slice_sample(n = 20) %>%
  ungroup()




