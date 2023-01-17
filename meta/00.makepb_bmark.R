# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generate pseudobulk and compositions info for 
#' Deconvolution benchmark
#' 

library(tidyverse)
library(scuttle)

include_cts <- c("CM", "Endo", "Fib", "Myeloid", "Lymphoid", "vSMCs", "PC")
data <- list()

# Broads dataset: Chaffin2022 --------------------------------------------------

# Importing pb data
chaffin_data <- read_csv("./data_HCMDCM_Nature/pb_data.csv",
                    show_col_types = FALSE)

colnames(chaffin_data)[1] <- "sample_id"

chaffin_data <- chaffin_data %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()

# Importing coldata of the matrices
chaffin_coldat <- read_csv("./data_HCMDCM_Nature/pb_coldata.csv",
                   show_col_types = FALSE)[,-1] %>%
  column_to_rownames("colname") %>%
  dplyr::rename(ncells = "counts")

# Filtering for cells of interest
chaffin_coldat <- chaffin_coldat[chaffin_coldat$cell_type %in% include_cts, ] 

# Filtering for profiles of interest
chaffin_data <- chaffin_data[, rownames(chaffin_coldat)]

# Generate new pseudobulk

chaffin_patient_pb <- scuttle::summarizeAssayByGroup(x = chaffin_data, 
                                                     ids = chaffin_coldat$donor_id, 
                                                     statistics = "sum") %>%
  assay(., "sum")

# Generate props
chaffin_patient_props <- chaffin_coldat %>%
  group_by(donor_id) %>%
  dplyr::mutate(pat_cells = sum(ncells)) %>%
  dplyr::mutate(cell_prop = ncells/pat_cells) %>%
  dplyr::select(donor_id, cell_type, cell_prop) %>%
  pivot_wider(names_from = cell_type, values_from = cell_prop)

chaffin_patient_pb <- chaffin_patient_pb[, chaffin_patient_props$donor_id]

data[["Chaffin2022"]] <- list("counts" = chaffin_patient_pb, 
                            "props" = chaffin_patient_props)

# Science dataset:Reichart2022 --------------------------------------------------

# Importing pb data
reichart_data <- readRDS("./data_DCMACM_Science/pb_red_mat.rds")

# Importing coldata of the matrices
reichart_coldat <- read_csv("./data_DCMACM_Science/pb_coldata.csv",
                   show_col_types = FALSE)[,-1] %>%
  dplyr::mutate(cell_type = strsplit(colname,"_") %>%
                  map_chr(., ~ .x %>% last())) %>%
  column_to_rownames("colname") %>%
  dplyr::rename(ncells = "counts")

# Filtering for cells of interest
reichart_coldat <- reichart_coldat[reichart_coldat$cell_type_uni %in% include_cts, ] 

# Filtering for profiles of interest
reichart_data <- reichart_data[, rownames(reichart_coldat)]

# Generate new pseudobulk

reichart_patient_pb <- scuttle::summarizeAssayByGroup(x = reichart_data, 
                                                     ids = reichart_coldat$donor_id, 
                                                     statistics = "sum") %>%
  assay(., "sum")

# Generate props
reichart_patient_props <- reichart_coldat %>%
  group_by(donor_id) %>%
  dplyr::mutate(pat_cells = sum(ncells)) %>%
  dplyr::mutate(cell_prop = ncells/pat_cells) %>%
  dplyr::select(donor_id, cell_type_uni, cell_prop) %>%
  pivot_wider(names_from = cell_type_uni, values_from = cell_prop,values_fill = 0)

reichart_patient_pb <- reichart_patient_pb[, reichart_patient_props$donor_id]

data[["Reichart2022"]] <- list("counts" = reichart_patient_pb, 
                              "props" = reichart_patient_props)

saveRDS(data, "./data_META/HF_studiespb.rds")

