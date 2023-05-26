# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Run MOFA+ models
library(tidyverse)
library(compositions)
library(MOFA2)

reticulate::use_condaenv(condaenv = "/Users/ricardoramirez/opt/miniconda3/envs/MOFA")

# The first layer of views in gene expression
gex_pb <- readRDS("./data_MI/mi_pb_red.rds")

# The second layer of views are cell compositions ---------------------------

# Builds a view matrix
view_to_matrix <- function(df) {
  
  df %>%
    pivot_wider(names_from = patient_region_id, 
                values_from = value) %>%
    column_to_rownames("feature") %>%
    as.matrix()
}

# Builds a view df
matrix_to_view <- function(mat, view_name) {
  mat %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "patient_region_id") %>%
    mutate(view = view_name)
}

# Make clr or ilr transformations from a matrix
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

ct_comps <- read_csv("./data_MI/momics_ctcomps.csv", 
                     show_col_types = FALSE) %>% 
  dplyr::select(-view) %>% 
  view_to_matrix() %>%
  build_compositions() %>%
  matrix_to_view(view_name = "clr_comps") %>%
  dplyr::rename("sample" = "patient_region_id") %>%
  dplyr::select(view, feature, sample, value)

# The third layer of views are interactions
spatial_ints <- read_csv("./data_MI/spatial_nicheints.csv",
                         show_col_types = FALSE) %>%
  dplyr::rename("sample" = "patient_region_id") %>%
  dplyr::select(view, feature, sample, value)

# Put all data in a single data frame

multiview_dat <- bind_rows(gex_pb, ct_comps, spatial_ints)

# Fit a MOFA model as before

MOFAobject <- MOFA2::create_mofa(multiview_dat)

data_opts <- MOFA2::get_default_data_options(MOFAobject)
train_opts <- MOFA2::get_default_training_options(MOFAobject)
model_opts <- MOFA2::get_default_model_options(MOFAobject)

# This avoids the regularization of multicellular programs per cell type.
# This avoids less sparse gene weights
model_opts$spikeslab_weights <- FALSE 

# Define the number of factors needed
model_opts$num_factors <- 6

# Prepare MOFA model:
MOFAobject <- MOFA2::prepare_mofa(object = MOFAobject,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)

outfile <- file.path("./results/MI/mofaplusmodel.hdf5")

model <- MOFA2::run_mofa(MOFAobject, outfile)

# Defining meta-data for downstream analysis
annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))

batch_info <- read_csv("./data_MI/snrna_batch_ann.csv") %>%
  select(orig.ident, batch) %>%
  unique()

meta_data <- read_csv("./data_MI/rna_patient_anns_revisions.csv",show_col_types = F) %>%
  left_join(annotation_names, by = "patient_group") %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic-enriched", "ischemic-enriched", "fibrotic-enriched"))) %>%
  unique() %>%
  left_join(batch_info, by = c("sample_id" = "orig.ident")) %>%
  dplyr::select(patient_region_id, major_labl, patient_group, batch) %>%
  unique()

UMAP_embedding <- MOFAcellulaR::plot_sample_2D(model = model,
                                               method = "UMAP",
                                               metadata = meta_data,
                                               sample_id_column = "patient_region_id",
                                               color_by = "patient_group")

group_assoc <- MOFAcellulaR::get_associations(model = model,
                                             metadata = meta_data,
                                             sample_id_column = "patient_region_id",
                                             test_variable = "patient_group",
                                             test_type = "categorical",
                                             group = FALSE)

batch_assoc <- MOFAcellulaR::get_associations(model = model,
                                             metadata = meta_data,
                                             sample_id_column = "patient_region_id",
                                             test_variable = "batch",
                                             test_type = "categorical",
                                             group = FALSE)

# Make final model plot

assoc_list <- list("group" = group_assoc,
                   "batch" = batch_assoc)

cplx_hmap <- plot_MOFA_hmap(model = model,
                            group = FALSE,
                            metadata = meta_data,
                            sample_id_column = "patient_region_id",
                            sample_anns = c("patient_group"),
                            assoc_list = assoc_list)

draw(cplx_hmap)

pdf("./results/lupus/lupus_summary.pdf", height = 5, width = 4)

draw(cplx_hmap)

dev.off()

