# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' I will map corrected states to bulk
#' in a collected disease score

library(decoupleR)
library(tidyverse)
library(scater)
library(scran)
library(edgeR)

# Import loadings and make them direction

loadings <- read_csv("./results/DCM_broad/MOFA_mcell/factor_desc/Factor1_char/loadings.csv",
                     show_col_types = F) %>%
  dplyr::mutate(sign = strsplit(celltype, "_") %>%
                  map_chr(., ~.x[[2]])) %>%
  dplyr::mutate(celltype = strsplit(celltype, "_") %>%
                  map_chr(., ~.x[[1]])) %>%
  dplyr::mutate(value = ifelse(sign == "pos",
                               value, value * -1))

# Check pseudobulk of visium

# Patient annotations ------------------------------------------------------
sample_dict <- read_csv("./data_MI/visium_patient_anns_revisions.csv")

# Read pseudobulk data of visium
pseudobulk_data <- readRDS("./data_MI/ps_integrated_slides.rds")[[1]][["gex"]]
assay(pseudobulk_data, "counts") <- assay(pseudobulk_data, "sum")

# Filtering
useful_genes <- edgeR::filterByExpr(pseudobulk_data, min.count = 10, min.prop = 0.85)
pseudobulk_data <- pseudobulk_data[useful_genes, ]

# Normalizing
all_nf <- edgeR::calcNormFactors(pseudobulk_data, method = "TMM")
sfs <- all_nf$samples$lib.size * all_nf$samples$norm.factors
pseudobulk_data <- scater::logNormCounts(pseudobulk_data, size_factors = sfs)
cols_pb_dat <- colData(pseudobulk_data)[,1]
pseudobulk_data <- assay(pseudobulk_data, "logcounts")
colnames(pseudobulk_data) <- cols_pb_dat

scaled_pb <- pseudobulk_data %>%
  t() %>%
  scale() %>%
  t()

# Decoupler

dR_run <- decoupleR::run_wmean(mat = scaled_pb, 
                               network = loadings,
                               .source = celltype, 
                               .target = gene,
                               .mor = value) %>%
  dplyr::filter(statistic == "wmean")

dR_run <- dR_run  %>%
  left_join(sample_dict %>%
              dplyr::select(sample_id, patient_group),
            by = c("condition" = "sample_id"))

ggplot(dR_run, aes(x = patient_group, y = score)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~source, ncol = 3)
  
# Do ReHeat

# Calculate disease scores for all samples ----------------------------------
reheat <- readRDS("./data_reheat/METAheart.rds")

dr_res <- map(reheat, function(x) {
  
  mat <- x$GEX %>%
    t() %>%
    scale() %>%
    t()
  
  dR_run <- decoupleR::run_wmean(mat = mat, 
                                 network = loadings,
                                 .source = celltype, 
                                 .target = gene,
                                 .mor = value) %>%
    dplyr::filter(statistic == "wmean")
  
  dR_run <- dR_run  %>%
    left_join(x$TARGETS %>%
                dplyr::select(Sample, HeartFailure),
              by = c("condition" = "Sample"))
  
}) %>%
  enframe() %>%
  unnest()

# Complete disease score
multicell_ds <- dr_res %>%
  group_by(name, source) %>%
  dplyr::mutate(scale_score = scale(score)[,1]) %>%
  dplyr::mutate(HeartFailure = factor(HeartFailure, 
                                      levels = c("yes", "no"))) %>%
  ungroup() %>%
  group_by(name, condition, HeartFailure) %>%
  summarise(mcell_ds = mean(scale_score)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(ds_t = map(data, function(dat) {
    
    t.test(mcell_ds ~ HeartFailure, data = dat) %>%
    broom::tidy()
    
  }))

# Summary stats

multicell_ds_stat <- multicell_ds %>%
  dplyr::select(name, ds_t) %>%
  unnest() %>%
  arrange(-statistic)


# Single cell type

ct_ds <- dr_res %>%
  group_by(name, source) %>%
  dplyr::mutate(scale_score = scale(score)[,1]) %>%
  dplyr::mutate(HeartFailure = factor(HeartFailure, 
                                      levels = c("yes", "no"))) %>%
  nest() %>%
  mutate(ds_t = map(data, function(dat) {
    
    t.test(scale_score ~ HeartFailure, data = dat) %>%
      broom::tidy()
    
  }))

ct_ds_stat <- ct_ds %>%
  dplyr::select(name, ds_t) %>%
  unnest() %>%
  arrange(-statistic)


ct_ds_stat %>%
  dplyr::select(source, name, statistic) %>%
  pivot_wider(values_from = statistic, names_from = name)
  



multicell_ds %>%
  dplyr::filter(name == "Kong10") %>%
  unnest(data) %>%
  ggplot(aes(x = HeartFailure, y = mcell_ds)) +
  geom_boxplot() +
  geom_point()












