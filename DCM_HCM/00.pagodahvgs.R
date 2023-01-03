library(tidyverse)
library(scITD)

# counts matrix
counts <- readRDS('./scITDdata/pb_snRNA_pats_mat.rds')

# meta data matrix
meta <- readRDS('./scITDdata/pb_snRNA_pats_meta.rds')[, c("patient_region_id", "cell_type", 
                                                                         "patient_group", "major_labl")] %>%
  dplyr::rename("donors" = patient_region_id, 
                "ctypes" = cell_type)

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

container <- clean_data(container, donor_min_cells = 5)

container <- get_pseudobulk(container)

container <- normalize_pseudobulk(container, method="trim", scale_factor=10000)

container <- get_normalized_variance(container)

ct_hvgs <- map(set_names(container$experiment_params$ctypes_use), function(ct) {
  
  norm_variances <- container$scMinimal_ctype[[ct]]$norm_variances
  norm_variances <- norm_variances[order(norm_variances,decreasing=TRUE)]
  # limit to overdispersed genes
  norm_variances <- norm_variances[norm_variances > 1.75]
  names(norm_variances)
  
} )

saveRDS(ct_hvgs, file = "./scITDdata/hvg_list.rds")




