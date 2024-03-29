# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Defines functions to perform a single MOFAcell
#' run
#' 
#' This pipeline starts from 2 elements:
#' 1. A count matrix with samples as columns and genes in rows
#' 2. A colData dataframe with at least 3 columns (cell_type, cell_counts, donor_id)
#' 
#' 

library(tidyverse)
library(scater)
library(scran)
library(uwot)
library(edgeR)

# Creates summarized experiment
create_init_exp <- function(counts, coldata) {
  
  pb_dat <- SummarizedExperiment(assays = list("counts" = counts), colData = DataFrame(coldata))
  
  return(pb_dat)
}


# Filter profiles 
# cts must be a vector with cell_type names
filt_profiles <- function(pb_dat, ncells = 50, cts) {
  
  # by n of cells 
  
  ix <- which(colData(pb_dat)[,"cell_counts"] >= ncells)
  pb_dat <- pb_dat[, ix]

  # by views of interest
  
  if(is.null(cts)) {
    
    cts <- set_names(colData(pb_dat)[,"cell_type"] %>%
                       unique()) 
    
  } else {
    
    cts <- purrr::set_names(cts)
    
  }
  
  pb_dat_list <- map(cts, function(ctype) { 
    
    ix <- which(colData(pb_dat)[,"cell_type"] == ctype)
    
    return(pb_dat[,ix])
    
  })
  
  return(pb_dat_list)
  
}

# Performs filtering of genes (lowly expressed)
filt_gex_byexpr <- function(pb_dat_list, min.count, min.prop) {
  
  pb_dat_red <- map(pb_dat_list, function(x) {
    
    useful_genes <- edgeR::filterByExpr(x, min.count = min.count, min.prop = min.prop)
    
    return(x[useful_genes, ])
  })
  
  return(pb_dat_red)
  
}

# Performs filtering of highly variable genes (after data transformation)
# If your prior lacks some cells then hvgs are estimated
filt_gex_byhvg <- function(pb_dat_list, prior_hvg = NULL, var.threshold = 1) {
  
  if(is.null(prior_hvg)) {
    
    pb_dat_red <- map(pb_dat_list, function(x) {
      hvg <- getTopHVGs(x,var.threshold = var.threshold)
      return(x[hvg, ])
    }) 
    
    return(pb_dat_red)
    
  } else {
    
    cts_in_data <- set_names(names(pb_dat_list))
    cts_in_prior <- set_names(names(prior_hvg))
    
    in_cts <- cts_in_data[cts_in_data %in% cts_in_prior]
    out_cts <- cts_in_data[!cts_in_data %in% cts_in_prior]
    
    in_cts_data <- pb_dat_list[in_cts]
    
    for(ct in in_cts) {
      ct_genes <- in_cts_data[[ct]] %>% rownames()
      ct_genes <- ct_genes[ct_genes %in% prior_hvg[[ct]]]
      in_cts_data[[ct]] <- in_cts_data[[ct]][ct_genes,]
    }
    
    if(length(out_cts) == 0) {
      
      return(in_cts_data)
      
    } else {
      
      out_cts_data <- pb_dat_list[out_cts]
      
      out_cts_data <- map(out_cts_data, function(x) {
        hvg <- getTopHVGs(x,var.threshold = var.threshold)
        return(x[hvg, ])
      }) 

      return(c(in_cts_data, out_cts_data))
      
    }
  }
}

# Performs filtering of highly variable genes (after data transformation)
# This is based on marker genes
# The assumption is that background gene expression can be traced
# by expression of cell type marker genes in cell types which shouldn't
# express the gene.
# In @prior_mrks one must provide a named list with marker genes defined by the user
# We will keep only marker genes in the hvgs if they are expressed in the expected cell type
filt_gex_bybckgrnd <- function(pb_dat_list, prior_mrks) {
  
  # Current genes per view
  ct_genes <- map(pb_dat_list, rownames) %>%
    enframe("view","gene") %>%
    unnest()
  
  prior_mrks_df <- prior_mrks %>%
    enframe("view_origin","gene") %>%
    unnest() %>%
    dplyr::mutate(marker_gene = TRUE)
  
  # Here are genes that aren't cell type markers
  ok_genes <- ct_genes %>%
    left_join(prior_mrks_df, by = "gene") %>%
    dplyr::filter(is.na(marker_gene)) %>%
    dplyr::select(view, gene)
  
  # Here are genes selected as HVG that are marker
  # genes, we will keep only genes if they appear
  # in the right cell
  not_bckground_genes <- ct_genes %>%
    left_join(prior_mrks_df, by = "gene") %>%
    na.omit() %>%
    unnest() %>%
    dplyr::filter(view == view_origin) %>%
    dplyr::select(view, gene)
  
  clean_hvgs <- bind_rows(ok_genes, 
                          not_bckground_genes) %>%
    group_by(view) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x[[1]])) %>%
    deframe()
  
  pb_dat_list <- pb_dat_list %>%
    filt_gex_byhvg(pb_dat_list = .,
                   prior_hvg = clean_hvgs,
                   var.threshold = NULL)
  
  return(pb_dat_list)
}


# Performs normalization via TMM
tmm_trns <- function(pb_dat_list, scale_factor = 1000000) {

  pb_dat_red <- map(pb_dat_list, function(x) {
    all_nf <- edgeR::calcNormFactors(x, method = "TMM")
    sfs <- all_nf$samples$lib.size * all_nf$samples$norm.factors
    pb <- sweep(assay(x, "counts"), MARGIN = 2, sfs, FUN = "/")
    assay(x, "logcounts") <- log1p(pb * scale_factor)
    
    return(x)
    
  })
  
  return(pb_dat_red)
  
}

# Makes a MOFA ready data set
pb_dat2MOFA <- function(pb_dat_list) {
  
  pb_red <- map(pb_dat_list, function(x) {
    
    dat <- assay(x, "logcounts") 
    
    colnames(dat) <- colData(x)[,"donor_id"]
    
    dat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "sample", values_to = "value")
    
  }) %>% 
    enframe(name = "view") %>%
    unnest() %>%
    dplyr::mutate(feature = paste0(view, "_", feature))
  
  return(pb_red)
  
}

pb_dat2long <- function(pb_dat_list) {
  
  pb_red <- map(pb_dat_list, function(x) {
    
    dat <- assay(x, "logcounts") 
    
    rest_info <- colData(x) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id")
    
    colnames(dat) <- rest_info[,"sample_id"]
    
    dat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "sample_id", values_to = "value") %>%
      left_join(rest_info, by = "sample_id")
    
  }) %>% 
    enframe(name = "view") %>%
    unnest()
  
  return(pb_red)
  
}



                       
                       