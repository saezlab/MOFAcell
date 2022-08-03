# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will generate utilities to deal with data
#' The basic structure is a feature/value/patient_region_id DF 

library(tidyverse)
library(compositions)

# Builds a view matrix
view_to_matrix <- function(df) {
  
  df %>%
    pivot_wider(names_from = patient_region_id, 
                values_from = value) %>%
    column_to_rownames("feature") %>%
    as.matrix()
}

# Builds a view from a matrix
# Features in rows, samples as columns
matrix_to_view <- function(mat, view_name) {
  mat %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "patient_region_id") %>%
    mutate(view = view_name)
}

# Runs PCA and returns a PCA list with PCs to patient_region_id
# And a explained variance vector

do_pca <- function(df, meta = NULL, features = NULL) {
  
  if(!is.null(features)) {
    sel_feats <- features[features %in% rownames(df)]
  } else {
    sel_feats <- rownames(df)
  }
  
  pca_obj <- prcomp(t(df)[, sel_feats], center = T, scale. = T)
  
  coords <- pca_obj$x %>%
    data.frame(check.names = F, stringsAsFactors = F) %>%
    rownames_to_column("patient_region_id") %>%
    as_tibble() %>%
    pivot_longer(-patient_region_id,
                 names_to = "PC", values_to = "value")
  
  if (!is.null(meta)) {
    coords <- left_join(coords, meta, by = "patient_region_id")
  }
  
  var <- round(summary(pca_obj)$importance[2, ] * 100, 2)
  
  res <- list()
  res$coords <- coords
  res$var <- tibble("PC" = names(var),
                    "Expl_var" = var)
  
  return(res)
}


# Runs PCA regression to a categorical variable included (patient_region_id, category)
regress_PCA <- function(pca_obj, cat_var = "patient_group") {
  
  pc_df <- pca_obj$coords %>%
    select_at(c("PC", "value", cat_var))
  
  pc_df %>% 
    group_by(PC) %>%
    nest() %>%
    mutate(anova_res = map(data, function(dat) {
      broom::tidy(aov(value ~ .,data = dat)) %>%
        dplyr::filter(term == cat_var)
    })) %>%
    dplyr::select(PC, anova_res) %>%
    unnest() %>%
    left_join(pca_obj$var) %>%
    ungroup() %>%
    dplyr::mutate(p_adj = p.adjust(p.value)) %>%
    dplyr::select(PC, p.value, Expl_var, p_adj)
  
}

# Builds recipes of ILR, CLR, ALR

build_compositions <- function(df, flavor = "clr") {
  
  # Features in columns:
  # Completes compositions
  comps <- compositions::acomp(t(df)) 
  
  if(flavor == "clr") {
    
    lcomp <- clr(comps) %>% as.matrix()
    colnames(lcomp) <- paste0("clr_", colnames(lcomp))
    
  } else if (flavor == "ilr") {
    
    lcomp <- ilr(comps) %>% as.matrix()
    colnames(lcomp) <- paste0("ilr_", seq(1, ncol(lcomp)))
    
  } else if (flavor == "alr") {
    
    lcomp <- alr(comps) %>% as.matrix()
    colnames(lcomp) <- paste0("alr_", colnames(lcomp))
    
  }
  
  return(t(lcomp))
  
}

