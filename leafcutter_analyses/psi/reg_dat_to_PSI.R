#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file_in  <- args[1]
file_out <- args[2]

# file_in <- '../../genome_wide/250916_filters/chr6/FastTopics_output/PCMT1/dat.tsv'

suppressPackageStartupMessages({
  library(tidyverse)
  library(fastTopics)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(purrr)
})

# -----------------------------
# Check for tiny input file (<100B)
# -----------------------------
if (file.size(file_in) < 100) {
  message("Input file too small: ", file_in, " -> writing empty PSI table")
  write.table(data.frame(), file_out, sep = "\t", quote = FALSE)
  quit(save = "no", status = 0)
}

# -----------------------------
# Load leafcutter clusters
# -----------------------------
clu <- fread("../leaf2_clustering/reg.junction_counts.unclassified.gz", select = 1)
clu <- clu[-1, ]
clu[, c("chr", "start", "end", "dot", "cluster", "strand") := tstrsplit(V1, "[:_]", fixed = FALSE)]
clu$end <- as.numeric(clu$end) + 1
clu[, cluster := sub("^clu_", "", cluster)]
clu_map <- clu[, .(chr, strand, junction = paste0(start, "-", end), cluster)]

# -----------------------------
# Load GENCODE annotation
# -----------------------------
gencode <- readRDS("../../inputs/gencode_v46_granges.RDS")

# -----------------------------
# Functions
# -----------------------------
# This makes a junction-cluster set for each gene based only on TCGA data
# and adds GTEx or TCGA additional junctions back to the TCGA-defined clusters
map_matrices <- function(file, clu_map, gencode) {
  library(data.table)
  library(dplyr)
  library(igraph)
  
  chrom <- sub(".*(chr[0-9XY]+).*", "\\1", file)
  gene_name <- basename(dirname(file))
  
  gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c("exon")]
  gene_strand <- unique(as.character(strand(gene_data)))
  if (length(gene_strand) != 1) {
    writeLines("ambiguous strand", file_out)
    quit(save="no", status = 0)
  }
  
  dat <- fread(file, data.table = F)
  if('V1' %in% colnames(dat)){
    dat <- column_to_rownames(dat, 'V1')
  }else{
    dat <- column_to_rownames(dat, 'rownames')
  }
  
  rail_ids <- rownames(dat)
  
  
  if(is_null(dim(dat))){
    writeLines("minimal dat input", file_out)
    quit(save="no", status = 0)
  }
  
  juncs <- names(dat)
  map <- data.table(junction = juncs)
  psi_table <- clu_map %>%
    filter(strand == gene_strand) %>%
    right_join(map, by = 'junction') %>%
    mutate(start = str_split_i(junction, '-', 1),
           end   = str_split_i(junction, '-', 2))
  
  
  final_cluster <- data.frame(psi_table %>%
                                filter(!is.na(cluster)))
  to_cluster <- data.frame(psi_table %>%
                             filter(is.na(cluster)))
  
  #add junctions to established clusters
  #WARNING: This can lead to junctions in up to two clusters
  matching_start <- to_cluster$start %in% final_cluster$start
  matching_end   <- to_cluster$end   %in% final_cluster$end
  
  to_add_idx <- which(matching_start | matching_end)
  
  if (length(to_add_idx) > 0) {
    
    new_rows <- lapply(to_add_idx, function(i) {
      
      clusters_that_match <- unique(c(
        final_cluster$cluster[final_cluster$start == to_cluster$start[i]],
        final_cluster$cluster[final_cluster$end   == to_cluster$end[i]]
      ))
      
      data.frame(
        chr     = chrom,                       # replace with real column names
        strand     = NA,
        junction = to_cluster$junction[i],
        cluster  = clusters_that_match,
        start    = to_cluster$start[i],
        end      = to_cluster$end[i],
        stringsAsFactors = FALSE
      )
    })
    
    final_cluster <- rbind(final_cluster, do.call(rbind, new_rows))
  }
  
  map_long <- final_cluster %>%
    select(junction, cluster)
  list(map_long = unique(map_long), rail_ids = rail_ids)
}

#This turns each matrix into a PSI matrix
compute_PSI_matrix <- function(file, rail_ids, map_long) {
  dat <- fread(file, data.table = F)
  
  # handle GTEx "rownames" column case
  if ("rownames" %in% colnames(dat)) {
    new_rail_ids <- dat$rownames
    dat <- dat[, colnames(dat) %in% map_long$junction, drop = FALSE]
    missing_cols <- setdiff(map_long$junction, colnames(dat))
    for (col in missing_cols) {
      dat[[col]] <- 0
    }
  }
  
  mat <- as.matrix(dat)
  mode(mat) <- "numeric"
  
  if (exists("new_rail_ids")) {
    rownames(mat) <- new_rail_ids
  } else if (!is.null(rail_ids)) {
    rownames(mat) <- rail_ids
  }
  
  psi_list <- list()
  colnames_list <- c()
  
  # loop over clusters
  for (cl in unique(map_long$cluster)) {
    cl_juncs <- filter(map_long, cluster == cl)$junction
    idx <- which(colnames(mat) %in% cl_juncs)
    if (length(idx) == 0) next
    
    cluster_counts <- rowSums(mat[, idx, drop = FALSE], na.rm = TRUE)
    
    # compute PSI for each junction in this cluster
    for (junc in cl_juncs) {
      j_idx <- which(colnames(mat) == junc)
      if (length(j_idx) == 0) next
      
      vals <- mat[, j_idx]
      psi_vals <- ifelse(cluster_counts > 0, vals / cluster_counts, NA_real_)
      
      # column name as "cluster:start-end"
      junc_coords <- paste0(cl, ":", junc)
      
      psi_list[[junc_coords]] <- psi_vals
      colnames_list <- c(colnames_list, junc_coords)
    }
  }
  
  if (length(psi_list) == 0) {
    return(matrix(nrow = 0, ncol = 0))
  }
  
  psi_mat <- do.call(cbind, psi_list)
  colnames(psi_mat) <- colnames_list
  
  if (exists("new_rail_ids")) {
    rownames(psi_mat) <- new_rail_ids
  } else if (!is.null(rail_ids)) {
    rownames(psi_mat) <- rail_ids
  }
  
  psi_mat
}

# -----------------------------
# Main
# -----------------------------
map_info <- map_matrices(file_in, clu_map, gencode)
psi_mat <- compute_PSI_matrix(file_in, map_info$rail_ids, map_info$map_long)

write.table(psi_mat, file_out, sep = "\t", quote = FALSE)

