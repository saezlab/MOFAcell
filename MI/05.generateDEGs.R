# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Generation of DEGs for each major cell type
#' in the MOFAcell model

library(MOFAcellulaR)
library(tidyverse)
library(edgeR)

# Read counts
mi_pb <- readRDS("./data_MI/mi_pb.rds")
counts <- assay(mi_pb, "counts")

# Get coldata
coldata <- colData(mi_pb) %>%
  as.data.frame() %>%
  dplyr::rename("donor_id" = "patient_region_id",
                "cell_counts" = "ncells") %>%
  dplyr::mutate(colid = paste0(cell_type, "_", donor_id)) %>%
  tibble::column_to_rownames("colid")

# Add patient info

annotation_names <- tibble(patient_group = c("group_1", "group_2", "group_3"),
                           patient_group_name = c("myogenic", "ischemic", "fibrotic"))

sample_dict <- read_csv("./data_MI/rna_patient_anns_revisions.csv",show_col_types = F) %>%
  left_join(annotation_names, by = "patient_group") %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic", "ischemic", "fibrotic"))) %>%
  dplyr::select(patient_region_id, patient_group) %>%
  unique()

coldata <- left_join(coldata, sample_dict, by = c("donor_id" = "patient_region_id"))

# Define same cell-types as in MOFAcell run
exclude_ct <- c("Mast", "Neuronal", "prolif", "Adipo")

cts <- coldata$cell_type %>%
  unique()

cts <- cts[! cts %in% exclude_ct]

# Make MOFA cell object
multiview_dat <- MOFAcellulaR::create_init_exp(counts = counts,  
                                               coldata = coldata) %>%
  MOFAcellulaR::filt_profiles(pb_dat = .,
                              cts = cts,
                              ncells = 25,
                              counts_col = "cell_counts", # This refers to the column name in testcoldata where the number of cells per profile was stored
                              ct_col = "cell_type") %>%
  MOFAcellulaR::filt_gex_byexpr(pb_dat_list = .,
                                min.count = 100,
                                min.prop = 0.25)

# Run edgeR
getDEGs <- function(view_dat) {
  
  dat <- DGEList(SummarizedExperiment::assay(view_dat, "counts"), 
                 samples = colData(view_dat),
                 group = colData(view_dat)[,"patient_group"])
  
  dat <- calcNormFactors(dat)
  
  design <- model.matrix(~ 0 + group, dat$samples)
  
  colnames(design) <- levels(dat$samples$group)
  
  dat <- estimateDisp(dat, design, robust = TRUE)
  
  fit <- glmQLFit(dat, design, robust = TRUE)
  
  my_contrasts <- makeContrasts(ischVSmyog = ischemic-myogenic, 
                                fibVSmyog = fibrotic-myogenic, 
                                ischVSfib = ischemic-fibrotic, 
                                levels = design)
  
  DEGs <- list(ischVSmyog = glmQLFTest(fit, contrast = my_contrasts[,"ischVSmyog"]),
               fibVSmyog = glmQLFTest(fit, contrast = my_contrasts[,"fibVSmyog"]),
               ischVSfib = glmQLFTest(fit, contrast = my_contrasts[,"ischVSfib"]))
  
  DEGs <- map(DEGs, function(deg) {
    
    topTags(deg, n = Inf) %>%
      as.data.frame() %>%
      rownames_to_column("gene")
    
  }) %>%
    enframe() %>%
    unnest(cols = c(value))
  
  return(DEGs)

}

# Run the model for all views and save the results

edgeR_allcts <- map(multiview_dat, getDEGs)

edgeR_allcts %>%
  enframe(name = "view") %>%
  unnest(cols = c(value)) %>%
  write_csv("./results/MI/classic_degs.csv")