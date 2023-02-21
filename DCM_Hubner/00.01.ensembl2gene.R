# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Transform ensenmbl IDs to gene names for the MOFA model

library(MOFA2)
library(tidyverse)
library(HDF5Array)
library(ComplexHeatmap)
library(scater)
library(scran)
library(uwot)
library(edgeR)
library(circlize)
library(compositions)
library(biomaRt)
source("./MOFAcell/code/MOFAcellprep.R")
source("./MOFAcell/code/factorutils.R")

# Importing pb data
pb_data <- read_csv("./data_DCMACM_Science/pb_data.csv",
                    show_col_types = FALSE)

colnames(pb_data)[1] <- "sample_id"

pb_data <- pb_data %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()


# Importing mart

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(pb_data)
gene_dict <- getBM(filters= "ensembl_gene_id", 
                   attributes= c("ensembl_gene_id","hgnc_symbol"),
                   values=genes,mart= mart)

# Transforming

pb_dat_red <- pb_data %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id) %>%
  left_join(gene_dict,  by = "ensembl_gene_id") %>%
  na.omit() %>%
  dplyr::select(-ensembl_gene_id) %>%
  group_by(name, hgnc_symbol) %>%
  summarise(sum_value = sum(value)) %>%
  pivot_wider(names_from = name, values_from = sum_value)

# Saving

pb_data <- pb_dat_red %>%
  dplyr::filter(hgnc_symbol != "") %>%
  column_to_rownames("hgnc_symbol") %>%
  as.matrix()

saveRDS(pb_data, "./data_DCMACM_Science/pb_red_mat.rds")



