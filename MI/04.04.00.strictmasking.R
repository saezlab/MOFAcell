# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we mask the decoupler normalized scores based on compositions

library(tidyverse)

compositions <- read_csv("./data_MI/cell_type_compositions.csv")

compositions <- compositions %>%
  dplyr::mutate(spot_id = gsub("[..]","_", spot_id)) %>%
  dplyr::mutate(spot_id = gsub("__","_", spot_id))

compositions <- compositions %>%
  dplyr::select(spot_id, name, value) %>%
  pivot_wider(names_from = name,values_from = value) %>%
  column_to_rownames("spot_id") %>%
  as.matrix()

# Mask strict version
# Only std masking for not-active scores
process_modules <- function(module_file, 
                            sd_thrsh = 2) {
  
  # This recovers the module scores
  module_mat <- read_csv(module_file,show_col_types = F) %>%
    as.data.frame()
  
  fit_module_file <- gsub("[.]csv", "_msk_strict.csv", module_file)
  
  rownames(module_mat) <- module_mat[,1]
  module_mat <- module_mat[, -1]
  
  # Filter compositions
  slide_comps <- compositions[rownames(module_mat),]
  
  # Mask module matrix
  module_mat[module_mat < sd_thrsh] = 0
  
  # Get ct order
  ct_order <- colnames(module_mat) %>% 
    strsplit(., "_") %>%
    map_chr(., ~.x[[1]])
  
  slide_comps <- slide_comps[rownames(module_mat), ct_order]
  colnames(slide_comps) <- colnames(module_mat)
  
  # Get masked module score matrix
  mask_module_mat <- module_mat * slide_comps
  
  write_csv(mask_module_mat %>% 
              as.data.frame() %>%
              rownames_to_column("spot_id"),
            fit_module_file)
  
  return(NULL)
  
}

dir <- "./results/MI/MOFA_mcell/factor_desc/Factor1_char/decoupler_ct/"

samples <- list.files(dir)
samples <- samples[!grepl("msk",samples)]
samples <- set_names(paste0(dir,samples), gsub("[.]csv", "", samples))

walk(samples, process_modules, sd_thrsh = 2)
