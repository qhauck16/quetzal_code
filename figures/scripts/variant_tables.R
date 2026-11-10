# Per-gene "full_info" / "full_plus_variants" tables: GSP loading per TCGA
# sample + catalogued mutations + self-called WXS / RNA-seq variants.
#
# Ported 2026-10-05 from tcga/genome_wide/flashier/flash_fact_char.Rmd so
# quetzal_code builds every table its figures read. Logic is the Rmd's,
# chunk by chunk (named below); only the three near-identical
# extract_{sf3b1,nkap,ylpm1}_variants() became extract_called_variants()
# with per-gene arguments. Validate against the Rmd's outputs after any
# change (see build.md files next to the drivers).
#
# Drivers (one per figure, each writes its own cache/):
#   figure2/scripts/build_cache_sf3b1_variants.R  -> SF3B1 (+ SUGP1)
#   figure3/scripts/build_cache_nkap_variants.R   -> NKAP / PRPF8 / FAM50A / SLU7
#   figure4/scripts/build_cache_ylpm1_variants.R  -> YLPM1
# Inputs: ../paths.yaml keys tcga_metadata, mutation_am_merged,
# alphamissense_refs_dir, unified_calls_dir, plus
# figures/shared_cache/gsp_loadings.tsv (scripts/build_gsp_loadings.R).
.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')
if (!exists("QROOT")) QROOT <- local({ d <- normalizePath(getwd())   # quetzal_code root
  while (!file.exists(file.path(d, ".quetzal_root"))) { if (dirname(d) == d) stop("run from inside quetzal_code/"); d <- dirname(d) }
  d })
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(tidyverse)
  library(VariantAnnotation)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
first <- dplyr::first

VT_FIG_ROOT <- file.path(QROOT, "figures")

# ---- shared inputs ("loading re-used files", "correlating flashier
#      factors with mutation status", "per-gene AlphaMissense lookup") ----
vt_inputs <- function() {
  tcga_metadata <- fread(file.path(QROOT, "inputs/tcga_v2_samples.tsv"), data.table = FALSE) %>%
    mutate(sampleId = sapply(strsplit(tcga_barcode, "-"), function(x) {
      x[4] <- sub("[A-Za-z]+$", "", x[4])   # drop trailing letters from 4th piece
      paste(x[1:4], collapse = "-")
    }))
  mut_data_am09 <- fread(file.path(QROOT, "mutation_processing/merged_grch38_with_am.tsv"), data.table = FALSE) %>%
    filter(is.na(am_pathogenicity) | am_pathogenicity > 0.9)
  am_lookup <- bind_rows(lapply(
    c("SLU7", "PRPF8", "NKAP", "FAM50A", "YLPM1", "SF3B1", "SUGP1"),
    function(g) {
      fread(file.path(file.path(QROOT, "mutation_processing/am_refs"),
                      paste0(g, "_alphamissense.tsv"))) %>%
        mutate(SYMBOL = g) %>%
        select(SYMBOL, protein_variant, am_pathogenicity)
    })) %>% distinct(SYMBOL, protein_variant, .keep_all = TRUE)
  L <- fread(file.path(VT_FIG_ROOT, "shared_cache", "gsp_loadings.tsv"),
             data.table = FALSE)
  list(tcga_metadata = tcga_metadata, mut_data_am09 = mut_data_am09,
       am_lookup = am_lookup, L = L)
}

# GSP column + project, as `<GENE>_factor` in the Rmd.
vt_factor <- function(I, k) {
  I$L[, c("rail_id", paste0("k", k))] %>%
    left_join(select(I$tcga_metadata, c(rail_id, gdc_cases.project.project_id)),
              by = "rail_id")
}

# Catalogued mutations for one gene, as `<GENE>_mutations` in the Rmd.
vt_mutations <- function(I, gene) {
  I$mut_data_am09 %>%
    filter(hugoGeneSymbol == gene) %>%
    left_join(select(I$tcga_metadata, c(rail_id, sampleId)), by = "sampleId") %>%
    select(rail_id, proteinChange, hugoGeneSymbol) %>%
    unique()
}

# ---- "Function to process a single VEP VCF" (verbatim) ----------------
# Handles both Mutect2 VCFs (AD/DP in FORMAT) and bcftools VCFs (DP4/DP in INFO)
process_vep_vcf <- function(vcf_file, genome = "hg38", file_id = NULL) {
  vcf <- readVcf(vcf_file, genome = genome)
  alts <- alt(vcf)
  alt_expanded <- unlist(alts)
  n_var <- length(alts)
  if (n_var == 0) {
    empty_final <- data.frame(chrom = character(), pos = integer(),
                              ref = character(), alt = character(),
                              filter = character(), vcf_file = character(),
                              local_variant_index = integer(),
                              stringsAsFactors = FALSE)
    empty_geno <- data.frame(local_variant_index = integer(),
                             sample = character(), GT = character(),
                             DP = integer(), alt_ad = integer(),
                             vcf_file = character(), stringsAsFactors = FALSE)
    return(list(final_vcf_df = empty_final, gt_long_df = empty_geno))
  }
  variant_index <- rep(seq_len(n_var), lengths(alts))
  rr <- rowRanges(vcf)
  filter_col <- fixed(vcf)$FILTER
  base_df <- data.frame(
    chrom = as.character(seqnames(rr))[variant_index],
    pos   = start(rr)[variant_index],
    ref   = as.character(ref(vcf))[variant_index],
    alt   = as.character(alt_expanded),
    filter = rep(filter_col, lengths(alts)),
    vcf_file = file_id,
    local_variant_index = variant_index,
    stringsAsFactors = FALSE)
  csq_header <- info(header(vcf))["CSQ", "Description"]
  if (!is.na(csq_header) && length(csq_header) == 1) {
    csq_fields <- strsplit(sub(".*Format: ", "", csq_header), "\\|")[[1]]
    csq <- info(vcf)$CSQ
    n_fields <- length(csq_fields)
    csq_df <- do.call(rbind, lapply(seq_along(csq), function(i) {
      if (length(csq[[i]]) == 0 || all(is.na(csq[[i]]))) return(NULL)
      entries <- unlist(strsplit(as.character(csq[[i]]), ",", fixed = TRUE))
      rows <- lapply(entries, function(e) {
        fields <- strsplit(e, "\\|", fixed = FALSE)[[1]]
        if (length(fields) < n_fields)
          fields <- c(fields, rep(NA_character_, n_fields - length(fields)))
        fields[seq_len(n_fields)]
      })
      df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
      colnames(df) <- csq_fields
      df$local_variant_index <- i
      df
    }))
    final_vcf_df <- cbind(base_df[csq_df$local_variant_index, ], csq_df)
  } else {
    final_vcf_df <- base_df
  }
  rownames(final_vcf_df) <- NULL
  gt <- geno(vcf)
  samples <- colnames(gt$GT)
  has_AD <- "AD" %in% names(gt)
  has_GT_DP <- "DP" %in% names(gt)
  if (!has_AD || !has_GT_DP) {
    info_data <- info(vcf)
    if (!has_GT_DP) info_dp <- info_data$DP
    if (!has_AD) info_dp4 <- info_data$DP4
  }
  gt_long <- lapply(samples, function(samp) {
    gt_vec <- as.character(gt$GT[, samp])
    n_alt <- lengths(alts)
    dp_vec <- if (has_GT_DP) gt$DP[, samp] else info_dp
    if (has_AD) {
      ad_list <- gt$AD[, samp]
      alt_ad <- unlist(lapply(seq_len(n_var), function(i) {
        if (is.null(ad_list[[i]])) return(rep(NA, n_alt[i]))
        ad_vals <- ad_list[[i]][-1]
        ad_vals[seq_len(n_alt[i])]
      }))
    } else {
      alt_ad <- unlist(lapply(seq_len(n_var), function(i) {
        dp4 <- info_dp4[[i]]
        if (is.null(dp4) || length(dp4) < 4) return(rep(NA, n_alt[i]))
        rep(dp4[3] + dp4[4], n_alt[i])
      }))
    }
    data.frame(local_variant_index = rep(seq_len(n_var), n_alt),
               sample = samp, GT = rep(gt_vec, n_alt),
               DP = rep(dp_vec, n_alt), alt_ad = alt_ad,
               vcf_file = file_id, stringsAsFactors = FALSE)
  })
  list(final_vcf_df = final_vcf_df, gt_long_df = bind_rows(gt_long))
}

# hgvs_protein_simple for MODERATE/HIGH calls (SF3B1 + NKAP chunks).
.hgvs_missense_indel <- function(d) {
  d %>% mutate(hgvs_protein_simple = case_when(
    Consequence == "missense_variant" & !is.na(Amino_acids) &
      str_detect(Amino_acids, "/") ~
      paste0("p.", str_split_fixed(Amino_acids, "/", 2)[, 1],
             Protein_position, str_split_fixed(Amino_acids, "/", 2)[, 2]),
    Consequence == "inframe_insertion" & !is.na(Amino_acids) &
      str_detect(Amino_acids, "/") ~ {
        aa <- str_split_fixed(Amino_acids, "/", 2)
        paste0("p.", aa[, 1], Protein_position, "ins",
               str_replace(aa[, 2], paste0("^", aa[, 1]), ""))
      },
    Consequence == "inframe_deletion" & !is.na(Amino_acids) &
      str_detect(Amino_acids, "/") ~ {
        aa <- str_split_fixed(Amino_acids, "/", 2)
        paste0("p.", aa[, 1], Protein_position, "del")
      },
    !is.na(HGVSp) & HGVSp != "" ~ str_replace(HGVSp, "^p\\.", "p."),
    TRUE ~ NA_character_))
}

# extract_{sf3b1,nkap,ylpm1}_variants(), merged. `mode`:
#   "moderate_high" (SF3B1, NKAP): IMPACT MODERATE/HIGH, then `region`
#   "ptc" (YLPM1): SYMBOL == ptc_symbol and a PTC consequence
# `transcript_match`: "prefix" = Feature without version %in% transcripts
# (SF3B1, NKAP); "grepl" = grepl(transcripts, Feature) (YLPM1).
extract_called_variants <- function(I, vcf_dir, source_label, apply_filter_pass,
                                    transcripts, transcript_match = "prefix",
                                    mode = "moderate_high", region = NULL,
                                    ptc_symbol = NULL, ptc_consequences = NULL) {
  vcf_files <- list.files(vcf_dir, pattern = "\\.vep\\.vcf$", full.names = TRUE)
  col_name <- paste0(source_label, "_mutation")
  if (length(vcf_files) == 0)
    return(tibble(rail_id = numeric(), !!col_name := character()))
  vcfs <- lapply(seq_along(vcf_files), function(i)
    process_vep_vcf(vcf_files[i], file_id = basename(vcf_files[i])))
  joined <- bind_rows(lapply(vcfs, function(el) {
    temp <- el$final_vcf_df
    if (nrow(temp) == 0) return(NULL)
    temp <- temp[, -ncol(temp)]
    j <- full_join(temp, el$gt_long_df, by = "local_variant_index")
    if (transcript_match == "prefix")
      filter(j, str_split_i(Feature, "\\.", 1) %in% transcripts)
    else filter(j, grepl(transcripts, Feature))
  }))
  passing <- joined %>%
    { if (apply_filter_pass) filter(., filter %in% c("PASS", "multiallelic")) else . }
  if (mode == "moderate_high") {
    passing <- passing %>%
      filter(IMPACT %in% c("MODERATE", "HIGH")) %>%
      filter(alt_ad / DP > 0.1) %>%
      filter(DP > 10) %>%
      .hgvs_missense_indel() %>%
      filter(region(.))
  } else {
    passing <- passing %>%
      filter(SYMBOL == ptc_symbol) %>%
      filter(Consequence %in% ptc_consequences) %>%
      filter(alt_ad / DP > 0.1) %>%
      filter(DP > 10) %>%
      distinct() %>%
      mutate(hgvs_protein_simple = case_when(
        Consequence %in% c("stop_gained", "frameshift_variant") &
          !is.na(Amino_acids) ~
          paste0("p.", substr(Amino_acids, 1, 1), Protein_position,
                 ifelse(Consequence == "frameshift_variant", "fs", "*")),
        !is.na(HGVSp) & HGVSp != "" ~ str_replace(HGVSp, "^p\\.", "p."),
        TRUE ~ NA_character_))
  }
  passing <- passing %>%
    mutate(am_key = if_else(
      Consequence == "missense_variant" & !is.na(Amino_acids) &
        str_detect(Amino_acids, "/"),
      paste0(str_split_fixed(Amino_acids, "/", 2)[, 1], Protein_position,
             str_split_fixed(Amino_acids, "/", 2)[, 2]),
      NA_character_)) %>%
    left_join(I$am_lookup, by = c("SYMBOL", "am_key" = "protein_variant")) %>%
    filter(is.na(am_pathogenicity) | am_pathogenicity > 0.9)
  prefix <- if (mode == "ptc") paste0(ptc_symbol, ":") else NULL
  passing %>%
    mutate(rail_id  = as.numeric(str_split_i(vcf_file.x, "\\.", 1)),
           mutation = if (is.null(prefix))
                        paste0(SYMBOL, ":", str_split_i(hgvs_protein_simple, "\\.", 2))
                      else paste0(prefix, str_split_i(hgvs_protein_simple, "\\.", 2))) %>%
    filter(!is.na(rail_id), !is.na(hgvs_protein_simple)) %>%
    distinct(rail_id, mutation) %>%
    group_by(rail_id) %>%
    summarise(!!col_name := paste(mutation, collapse = ","), .groups = "drop")
}

# Collapse per-(rail, gene) mutation strings into one "GENE:a,b; GENE2:c"
# column -- the two-level summarise every gene chunk uses.
.collapse_gene_mutation <- function(d, kcol, empty_value) {
  d %>%
    group_by(rail_id, hugoGeneSymbol) %>%
    summarise(
      !!kcol := dplyr::first(.data[[kcol]]),
      project_id = dplyr::first(`gdc_cases.project.project_id`),
      hotspot = any(hotspot),
      mutations = ifelse(all(is.na(proteinChange)), NA,
                         paste(unique(proteinChange[!is.na(proteinChange)]), collapse = ",")),
      .groups = "drop_last") %>%
    summarise(
      !!kcol := dplyr::first(.data[[kcol]]),
      project_id = dplyr::first(project_id),
      hotspot = any(hotspot),
      gene_mutation = paste(
        ifelse(is.na(hugoGeneSymbol), NA, paste0(hugoGeneSymbol, ":", mutations)),
        collapse = "; "),
      .groups = "drop") %>%
    mutate(gene_mutation = ifelse(gene_mutation == "NA; ", empty_value, gene_mutation)) %>%
    arrange(desc(.data[[kcol]]))
}

.add_calls <- function(full_info, wxs_vars, rna_vars, backfill_col) {
  full_info %>%
    left_join(wxs_vars, by = "rail_id") %>%
    left_join(rna_vars, by = "rail_id") %>%
    mutate(across(where(is.character), ~ ifelse(. == "NA", NA, .))) %>%
    mutate(!!backfill_col := case_when(
      !is.na(.data[[backfill_col]]) & .data[[backfill_col]] != "" ~ .data[[backfill_col]],
      !is.na(wxs_mutation)    ~ wxs_mutation,
      !is.na(rnaseq_mutation) ~ rnaseq_mutation,
      TRUE                    ~ .data[[backfill_col]]))
}

# ---- SF3B1: "SF3B1 cdf", "SUGP1 CDFs", "adding variants to SF3B1 df" ----
build_sf3b1_tables <- function(I, k = 10) {
  kcol <- paste0("k", k)
  sf <- left_join(vt_factor(I, k), vt_mutations(I, "SF3B1"), by = "rail_id")
  sf$hotspot <- sf$proteinChange %in% c(
    'R625H', 'R625C', 'K700E', 'K666T', 'G740E', 'K666T', 'K666E', 'K666N',
    'K741E', 'K741N', 'K741Q', 'G740V')
  sf <- .collapse_gene_mutation(sf, kcol, "")
  mutation_positions <- lapply(str_extract_all(sf$gene_mutation, "\\d+"), as.numeric)
  near_hotspots <- vapply(mutation_positions,
                          function(pos) any(pos > 620 & pos < 800, na.rm = TRUE), logical(1))
  near_hotspots[is.na(near_hotspots)] <- FALSE
  sf$plot_color <- "No SF3B1 Mutation"
  has_sf3b1_mut <- sf$gene_mutation != 'NA'
  sf$plot_color[has_sf3b1_mut] <- "Other SF3B1 Mutation"
  sf$plot_color[near_hotspots & !sf$hotspot] <- "Mutation Near SF3B1 Hotspot"
  sf$plot_color[sf$hotspot] <- "SF3B1 hotspot mutation"
  sf[[paste0("logk", k)]] <- log(sf[[kcol]])
  SF3B1_full_info <- sf

  sf$proteinChange <- str_split_i(sf$gene_mutation, ':', 2)
  SUGP1_mutations <- vt_mutations(I, "SUGP1")
  combined <- bind_rows(vt_mutations(I, "SF3B1"), SUGP1_mutations) %>%
    group_by(rail_id, hugoGeneSymbol) %>%
    summarise(mutations = ifelse(all(is.na(proteinChange)), NA,
                                 paste(unique(proteinChange[!is.na(proteinChange)]), collapse = ",")),
              .groups = "drop_last") %>%
    summarise(gene_mutation = paste(
      ifelse(is.na(hugoGeneSymbol), NA, paste0(hugoGeneSymbol, ":", mutations)),
      collapse = "; "), .groups = "drop") %>%
    mutate(gene_mutation = ifelse(gene_mutation == "NA; ", NA, gene_mutation))
  sugp1_lookup <- SUGP1_mutations %>%
    group_by(rail_id) %>%
    summarise(SUGP1_mutation = ifelse(all(is.na(proteinChange)), NA_character_,
                paste(unique(proteinChange[!is.na(proteinChange)]), collapse = ",")),
              .groups = "drop")
  SF3B1_SUGP1_full_info <- sf %>%
    select(-any_of("gene_mutation")) %>%
    left_join(combined, by = "rail_id") %>%
    left_join(sugp1_lookup, by = "rail_id")

  region <- function(d) (d$SYMBOL == "SF3B1" & d$Protein_position > 500 &
                           d$Protein_position < 800) | d$SYMBOL == "SUGP1"
  args <- list(I = I, transcripts = c("NM_012433", "NM_172231"), region = region)
  vd <- file.path(file.path(QROOT, "variant_recalling/unified"), "sf3b1", "vep")
  wxs <- do.call(extract_called_variants, c(args, list(vcf_dir = file.path(vd, "wxs"),
                 source_label = "wxs", apply_filter_pass = TRUE)))
  rna <- do.call(extract_called_variants, c(args, list(vcf_dir = file.path(vd, "rna"),
                 source_label = "rnaseq", apply_filter_pass = FALSE)))
  list(SF3B1_full_info = SF3B1_full_info,
       SF3B1_SUGP1_full_info = SF3B1_SUGP1_full_info,
       SF3B1_full_plus_variants = .add_calls(SF3B1_SUGP1_full_info, wxs, rna, "gene_mutation"))
}

# ---- NKAP: "processing variants" --------------------------------------
build_nkap_tables <- function(I, k = 13) {
  kcol <- paste0("k", k)
  all_mutations <- unique(bind_rows(lapply(c("NKAP", "PRPF8", "FAM50A", "SLU7"),
                                           function(g) vt_mutations(I, g))))
  nk <- left_join(vt_factor(I, k), all_mutations, by = "rail_id")
  nk$hotspot <- nk$hugoGeneSymbol == "NKAP" &
    nk$proteinChange %in% c('R368Q', 'R402*', 'R402Sfs*13')
  k_full_info <- .collapse_gene_mutation(nk, kcol, "")
  NKAP_full_info <- k_full_info %>%
    mutate(
      genes = str_extract_all(gene_mutation, "[A-Z0-9]+(?=:)"),
      positions = str_extract_all(gene_mutation, "(?<=[:,])[A-Z]*([0-9]+)") %>%
        lapply(function(x) as.numeric(str_extract(x, "[0-9]+"))),
      plot_color = case_when(
        map_lgl(seq_along(genes), ~ "NKAP" %in% genes[[.x]] &&
                  any(positions[[.x]] %in% c(368, 402))) ~ "NKAP hotspot mutation",
        map_lgl(seq_along(genes), ~ "NKAP" %in% genes[[.x]] &&
                  any(positions[[.x]] >= 320 & positions[[.x]] <= 420)) ~ "Near NKAP Hotspot",
        map_lgl(seq_along(genes), ~ "PRPF8" %in% genes[[.x]] &&
                  any(positions[[.x]] >= 1570 & positions[[.x]] <= 1900)) ~ "PRPF8-RH mutation",
        map_lgl(genes, ~ "FAM50A" %in% .x) ~ "FAM50A",
        map_lgl(genes, ~ any(.x %in% c("NKAP", "PRPF8", "FAM50A"))) ~
          "Other Mutation in NKAP/PRPF8/FAM50A",
        TRUE ~ "No Mutation")) %>%
    select(-genes, -positions)
  NKAP_full_info[[paste0("logk", k)]] <- log(NKAP_full_info[[kcol]])

  region <- function(d) d$SYMBOL == "NKAP" & d$Protein_position > 300 |
    d$SYMBOL == "PRPF8" & d$Protein_position >= 1570 & d$Protein_position <= 1900 |
    d$SYMBOL == "FAM50A" | d$SYMBOL == "SLU7"
  args <- list(I = I, transcripts = c("NM_024528", "NM_004699", "NM_006445", "NM_006425"),
               region = region)
  vd <- file.path(file.path(QROOT, "variant_recalling/unified"), "nkap", "vep")
  wxs <- do.call(extract_called_variants, c(args, list(vcf_dir = file.path(vd, "wxs"),
                 source_label = "wxs", apply_filter_pass = TRUE)))
  rna <- do.call(extract_called_variants, c(args, list(vcf_dir = file.path(vd, "rna"),
                 source_label = "rnaseq", apply_filter_pass = FALSE)))
  list(k13_full_info = k_full_info,
       NKAP_full_plus_variants = .add_calls(NKAP_full_info, wxs, rna, "gene_mutation"))
}

# ---- YLPM1: "YLPM1 parsing vcf" ---------------------------------------
build_ylpm1_tables <- function(I, k = 30) {
  kcol <- paste0("k", k)
  yl <- left_join(vt_factor(I, k), vt_mutations(I, "YLPM1"), by = "rail_id")
  yl$frameshift <- grepl('\\*', yl$proteinChange)
  yl <- yl %>%
    group_by(rail_id) %>%
    summarise(!!kcol := dplyr::first(.data[[kcol]]),
              project_id = dplyr::first(`gdc_cases.project.project_id`),
              hugoGeneSymbol = dplyr::first(hugoGeneSymbol),
              frameshift = any(frameshift),
              proteinChange = paste(unique(proteinChange), collapse = "; "))
  yl$plot_color <- 'No YLPM1 Mutation'
  yl$plot_color[(yl$proteinChange != 'NA') & !yl$frameshift] <- 'Other YLPM1 Mutation'
  yl[yl$frameshift, "plot_color"] <- 'YLPM1 Frameshift/Stop Gained'
  yl[[paste0("logk", k)]] <- log(yl[[kcol]])

  args <- list(I = I, transcripts = "NM_019589", transcript_match = "grepl",
               mode = "ptc", ptc_symbol = "YLPM1",
               ptc_consequences = c("stop_gained", "frameshift_variant",
                                    "splice_acceptor_variant", "splice_donor_variant"))
  vd <- file.path(file.path(QROOT, "variant_recalling/unified"), "ylpm1", "vep")
  wxs <- do.call(extract_called_variants, c(args, list(vcf_dir = file.path(vd, "wxs"),
                 source_label = "wxs", apply_filter_pass = TRUE)))
  rna <- do.call(extract_called_variants, c(args, list(vcf_dir = file.path(vd, "rna"),
                 source_label = "rnaseq", apply_filter_pass = FALSE)))
  list(YLPM1_full_info = yl,
       YLPM1_full_plus_variants = .add_calls(yl, wxs, rna, "proteinChange"))
}

# Write each table to <cache_dir>/<name>.tsv (readr, as the Rmd did).
write_variant_tables <- function(tabs, cache_dir) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  for (n in names(tabs)) {
    f <- file.path(cache_dir, paste0(n, ".tsv"))
    readr::write_tsv(tabs[[n]], f)
    cat(sprintf("[variant_tables] wrote %s (%d rows)\n", f, nrow(tabs[[n]])))
  }
}

# Read one back with the mutation columns forced to character: they are
# ~all NA and readr would guess logical and silently NA every call.
read_variant_table <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE,
                  col_types = readr::cols(.default = readr::col_guess(),
                    gene_mutation = "c", wxs_mutation = "c", rnaseq_mutation = "c",
                    proteinChange = "c", SUGP1_mutation = "c", hugoGeneSymbol = "c"))
}
