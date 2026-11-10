library(tidyverse)
library(readr)
library(fastTopics)
library(nnls)
library(GenomicRanges)

args <- commandArgs(TRUE)

input_file <- args[1]
meta_file <- args[2]
threads <- as.numeric(args[3])
min_samples_per_junc <- as.numeric(args[4])
outdir <- args[5]
outprefix <- args[6]
elbow_cutoff <- as.numeric(args[7])
max_factors <- as.numeric(args[8])
avg_clust_read_count_min <- as.numeric(args[9])

#to use to remove junctions that are far away from annotated gene model
gencode <- readRDS(args[10])
avg_sample_reads_per_cluster <- as.numeric(args[11])
# maximum allowed junction coordinate away from gencode established gene boundaries
gene_range_bound <- as.numeric(args[12])


# input_file <- 'all_genes/chr1/snaptron_output/ACAP3_snaptron.tsv'
# meta_file <- 'gtexv2_metadata.tsv'
# threads <- 22
# min_samples_per_junc <- 5
# outdir <- 'testing/'
# outprefix <- 'ACAP3'
# elbow_cutoff <- 0.0001
# max_factors <- 32
# avg_clust_read_count_min <- 10
# gene_range_bound <- 2000
# #to use to remove junctions that are far away from annotated gene model
# gencode <- readRDS('hg38_granges.RDS')
# avg_sample_reads_per_cluster <- 5


gene_name <- str_split_i(basename(input_file), '_', 1)
gene_data <- gencode[gencode$gene_name == gene_name & gencode$type %in% c("exon")]
t_models <- rtracklayer::split(gene_data, gene_data$transcript_id)

lower_bound <- min(t_models@unlistData@ranges@start)
upper_bound <- max(as.numeric(t_models@unlistData@ranges@start) + as.numeric(t_models@unlistData@ranges@width))

if(!dir.exists(outdir)){
  dir.create(outdir)
}

if (!dir.exists(paste0(outdir, '/', outprefix))){
  dir.create(paste0(outdir, '/', outprefix))
}

if (file.size(input_file) < 100L){
  #do nothing
  
  saveRDS(list(), paste0(outdir, '/', outprefix, '/', 'dat.tsv'))
  
}else{
  
  tcga_metadata <- read_tsv(meta_file)
  
  min_junctions_for_pca <- 3
  
  working_dir <- getwd()
  
  gene_table <- read_tsv(input_file)
  
  #at least one gene has duplicates for some reason
  #Get rid of any junctions that probably don't actually correspond to the gene
  #get rid of junctions on opposite strand
  gene_table <- gene_table %>%  
    filter(samples_count > min_samples_per_junc) %>%
    unique() %>%
    filter(start > (lower_bound - gene_range_bound)) %>%
    filter(end < (upper_bound + gene_range_bound)) %>% 
    filter(strand %in% unique(as.character(strand(gene_data))))
  
  
  
  if(nrow(gene_table) < min_junctions_for_pca){
    
    saveRDS(list(), paste0(outdir, '/', outprefix, '/', 'dat.tsv'))
  }else{
    
    #functions to parse the snaptron samples column
    split_to_ids <- function(tem){
      return(str_split(tem, ':')[[1]][1])
    }
    split_to_counts <- function(tem){
      return(str_split(tem, ':')[[1]][2])
    }
    
    #setup matrix of samples by junction counts
    just_counts <- matrix(0, nrow = nrow(gene_table), ncol = length(tcga_metadata$rail_id))
    colnames(just_counts) <- tcga_metadata$rail_id
    rownames(just_counts) <- as.character(1:nrow(just_counts))
    
    #updating just_counts with counts corresponding to junctions of interest
    for (i in 1:nrow(gene_table)){
      pre_split <- str_split(gene_table$samples[i], ',')
      pre_split <- pre_split[[1]][2:length(pre_split[[1]])]
      
      rail_ids <- as.vector(sapply(pre_split, split_to_ids))
      counts <- as.vector(sapply(pre_split, split_to_counts))
      
      just_counts[cbind(i, rail_ids)] <- counts
    }
    
    #ensure proper behavior of matrix as numeric
    colnames_to_keep <- colnames(just_counts)
    just_counts <- matrix(as.numeric(just_counts), nrow = nrow(just_counts), ncol = ncol(just_counts))
    colnames(just_counts) <- colnames_to_keep
    
    
    #filter out junctions that previously passed due to many normal samples having counts of them
    gene_table <- gene_table[!rowSums(just_counts > 0) < min_samples_per_junc, ]
    just_counts <- just_counts[!rowSums(just_counts > 0) < min_samples_per_junc, ]
    
    
    junction_by_count <- cbind(gene_table, just_counts)
    
    
    # leafcutter-esque clustering ------------------------------------------------------
    
    clustering <- junction_by_count %>%
      group_by(start) %>%
      mutate(ID = cur_group_id())
    
    
    clusters <- clustering$ID
    
    for (i in 2:max(clusters)){
      ends <- clustering$end[clusters==i]
      if (sum(ends %in% clustering$end[clusters < i]) > 0){
        index_of_matching_end <- c(1:length(clusters))[(clustering$end %in% ends & clustering$ID < i)][1]
        clusters[clusters == i] = clusters[index_of_matching_end]
      }
    }
    
    junction_by_count$cluster <- clusters
    
    intron_cluster <- junction_by_count %>%
      group_by(cluster)%>%
      mutate(across(colnames(just_counts), ~sum(.x), .names = "{.col}_clust_sum"))
    
    #end of leafcutter-esque clustering--------------------------------------------
    
    #remove clusters that do not have an average of at least x reads 
    #retain only junction info columns, counts and cluster column
    clust_sums <- intron_cluster[,((ncol(intron_cluster)-1)-ncol(just_counts)+1):ncol(intron_cluster)] %>%
      unique()
    #kept cluster to avoid throwing away clusters w/same dist of junctions (unlikely), now need to remove it
    clust_sums <- subset(clust_sums, select = -cluster)
    clust_to_remove <- (rowSums(clust_sums)/ncol(clust_sums) < avg_clust_read_count_min)
    filtered_intron_cluster <- intron_cluster %>% 
      filter(!cluster %in% unique(clusters)[clust_to_remove]) %>% 
      .[,(1:ncol(junction_by_count))]
    
    #remove any very low coverage clusters prior to removing low coverage samples
    clust_sums <- clust_sums[!clust_to_remove, ]
    
    #need to filter out genes that all clusters are thrown away
    if(nrow(clust_sums) == 0){
      saveRDS(list(), paste0(outdir, '/', outprefix, '/', 'dat.tsv'))
    }else{
      
      junction_names <- paste0(filtered_intron_cluster$start, '-', filtered_intron_cluster$end)
      
      #remove samples that do not have an average of at least 5 reads per cluster
      samples_to_remove <- (colSums(clust_sums)/nrow(clust_sums) < avg_sample_reads_per_cluster)
      filtered_intron_cluster <- filtered_intron_cluster[,c(rep(F, 18),!as.vector(samples_to_remove), T)]
      
      
      #We may have no samples that fit this criterion
      if(ncol(filtered_intron_cluster) == 1){
        saveRDS(list(), paste0(outdir, '/', outprefix, '/', 'dat.tsv'))
      }else{
      
        #remove any junctions that in final dat have no support, not counting cluster column in RowSums
        
        junction_names <- junction_names[(rowSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]) >  0)]
        filtered_intron_cluster <- filtered_intron_cluster[(rowSums(filtered_intron_cluster[1:(ncol(filtered_intron_cluster)-1)]) > 0),]
                
        #want something actually relevant to perform PCA on
        if (nrow(filtered_intron_cluster) < min_junctions_for_pca || ncol(filtered_intron_cluster)
            < min_junctions_for_pca){

          saveRDS(list(), paste0(outdir, '/', outprefix, '/', 'dat.tsv'))
          
        }else{
          
          #want samples to at least have a read, highly unlikely this happens but possible if
          #a sample only has reads in an unsupported junction
          filtered_intron_cluster <- filtered_intron_cluster[, colSums(filtered_intron_cluster) > 0]
          
        
          dat <- t(as.matrix(filtered_intron_cluster[,1:(ncol(filtered_intron_cluster)-1)]))
          
          
          dat <- as.data.frame(dat)
          colnames(dat) <- junction_names
          df_out <- cbind(rownames = rownames(dat), dat)
          write_tsv(df_out, paste0(outdir, outprefix, '/', 'dat.tsv'))
        }
      }
    }
  }
}
  
  
  
  
