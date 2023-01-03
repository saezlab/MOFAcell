# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will use scITD to
#' use PAGODA's implementation of highly variable 
#' genes to be able to compare MOFAcell with
#' the decomposition of scITD

library(tidyverse)
library(scITD)

# counts matrix
counts <- readRDS('./scITDdata/pb_snRNA_pats_mat.rds')

# Defining the patient annotations -----------------------------------------------------------------
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
  unique()

# meta data matrix
meta <- readRDS('./scITDdata/pb_snRNA_pats_meta.rds')[, c("patient_region_id", "cell_type")] %>%
  rownames_to_column("cell_id") %>%
  left_join(sample_dict_red, by = "patient_region_id") %>%
  dplyr::rename("donors" = patient_region_id, 
                "ctypes" = cell_type) %>%
  column_to_rownames("cell_id")

# set up project parameters
# I will exclude lowly abundant cell-types
param_list <- initialize_params(ctypes_use = c("Fib", "CM", "Endo", "Myeloid", "PC", "vSMCs", "Lymphoid"),
                                ncores = 4, 
                                rand_seed = 10)

# create project container
container <- make_new_container(count_data=counts, 
                                meta_data=meta,
                                params=param_list,
                                label_donor_sex = FALSE)


# Do scITD manual processing

container <- parse_data_by_ctypes(container)

container <- clean_data(container, donor_min_cells = 25)

container <- get_pseudobulk(container)

container <- normalize_pseudobulk(container, method="trim", scale_factor=1000000)

container <- get_normalized_variance(container)

ct_hvgs <- map(set_names(container$experiment_params$ctypes_use), function(ct) {
  
  norm_variances <- container$scMinimal_ctype[[ct]]$norm_variances
  norm_variances <- norm_variances[order(norm_variances,decreasing=TRUE)]
  # limit to overdispersed genes
  norm_variances <- norm_variances[norm_variances > 1.5] # 1.75
  names(norm_variances)
  
} )

saveRDS(ct_hvgs, file = "./scITDdata/hvg_list.rds")




