# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script I will generalize different tasks
#' that allow to explore an individual Factor within
#' a MOFAcell run

library(MOFA2)
library(tidyverse)
library(ggpubr)

#
model_outfile <- "./results/MOFA_mcell/MI_model.hdf5"
model <- MOFA2::load_model(model_outfile)

annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

batch_info <- read_csv("./data_MI/snrna_batch_ann.csv") %>%
  select(orig.ident, batch) %>%
  unique()

sample_dict <- read_csv("./data_MI/rna_patient_anns_revisions.csv",show_col_types = F) %>%
  left_join(annotation_names, by = "patient_group") %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))) %>%
  #dplyr::select(-sample_id) %>%
  unique()

sample_dict_red <- sample_dict %>%
  left_join(batch_info, by = c("sample_id" = "orig.ident")) %>%
  dplyr::select(patient_region_id, major_labl, patient_group, batch) %>%
  unique() %>%
  dplyr::rename("sample" = patient_region_id)

meta <- sample_dict_red

factor <- "Factor2"

# 0. Get scores and loadings of factor of interest

#' @param model = MOFAcell model
#' @param meta = a data frame with sample + any other colums
#' @param factor = Factor# label
#' 
#' returns a df with the factor scores
get_fscores <- function(model, meta, factor) {
  
  factor_scores <- get_factors(model, factors = "all")[[1]][,factor, drop= F] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample")
  
  return(factor_scores)
  
  
}

#' @param model = MOFAcell model
#' @param factor = Factor# label
#' 
#' returns a df with the factor loadings

get_floadings <- function(model, factor) {
  
  factor_loadings <- get_weights(model, as.data.frame = T) %>%
    as.data.frame() %>%
    dplyr::mutate(feature = strsplit(as.character(feature), "_") %>%
                    map_chr(., ~ .x[[2]]),
                  ctype = strsplit(as.character(view), "_") %>%
                    map_chr(., ~ .x[[1]])) %>%
    dplyr::rename("factors" = factor) %>%
    dplyr::select(-view) %>%
    dplyr::filter(factors == factor) %>%
    dplyr::select(-factors)
  
  return(factor_loadings)
  
}

# 1. Associate with sample grouping variable

#' @param factor_scores = get_fscores output
#' @param factor = Factor# label
#' @param covar = string with column to be used for pairwise comparisons
#' 
#' returns a df with the factor loadings

get_pw_scorecomp <- function(factor_scores, factor, covar) {
  
  combinations <- combn(factor_scores[,covar] %>%
                          unique,2, simplify = F)
  
  t_comps <- map(combinations, function(comb) {
    
    # Identify proper comparison for unification
    test_term <- levels(comb)
    test_term <- test_term[test_term %in% comb]
    reference_term <- test_term[2]
    test_term <- test_term[1]
    
    comb_scores <- factor_scores %>%
    dplyr::filter(if_any(.cols = all_of(covar),
                         .fns = ~ .x %in% as.character(comb)))
    
    res_tibble <- t.test(as.formula(paste0(factor, " ~ ", covar)),
                      comb_scores) %>%
             broom::tidy() %>%
      dplyr::select(statistic, p.value) %>%
      dplyr::mutate(reference = reference_term,
                    test = test_term)

  }) %>%
    enframe() %>%
    unnest(cols = c(value)) %>%
    dplyr::select(-name) %>%
    dplyr::mutate(p_adj = p.adjust(p.value))
  
  return(t_comps)
  
}

# 2. Correlate gene expression to factors

get_filtered_loadings <- function(factor_scores, get_floadings, pb_data, ) {
  
  
  
  
  
}

gene_cor <- mi_pb_red %>%
  left_join(factor_score) %>%
  group_by(view, feature) %>%
  nest() %>%
  mutate(cor_res = map(data, function(dat) {
    cor.test(dat$value, dat$Factor2) %>%
      broom::tidy()
  })) %>%
  dplyr::select(view, feature, cor_res) %>%
  unnest(c(cor_res)) %>%
  ungroup() %>%
  dplyr::mutate(adj_p = p.adjust(p.value))





















