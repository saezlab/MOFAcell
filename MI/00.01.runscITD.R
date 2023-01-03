# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will use scITD to
#' fit multicellular factor analysis

library(tidyverse)
library(scITD)
library(ComplexHeatmap)

# counts matrix
pb_counts <- readRDS('./scITDdata/pb_snRNA_pats_mat.rds')

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
  rownames_to_column("cellid") %>%
  left_join(sample_dict_red, by = "patient_region_id") %>%
  dplyr::rename("donors" = patient_region_id, 
                "ctypes" = cell_type) %>%
  column_to_rownames("cellid")

# hvgs from PAGODA
# pagoda_hvgs <- readRDS("./scITDdata/hvg_list.rds") %>% unlist() %>% unique()
pagoda_hvgs <- readRDS("./data_MI/mi_pb_red.rds") %>% 
  pull(feature) %>% 
  strsplit(.,"_") %>%
  map_chr(., ~ .x[[2]]) %>%
  unique()
  
# set up project parameters
# I will exclude lowly abundant cell-types
param_list <- initialize_params(ctypes_use = c("Fib", "CM", "Endo", 
                                               "Myeloid", "PC", "vSMCs", 
                                               "Lymphoid"),
                                ncores = 4, 
                                rand_seed = 10)

# create project container
container <- make_new_container(count_data=pb_counts, 
                                meta_data=meta,
                                params=param_list,
                                label_donor_sex = FALSE)


rm(pb_counts)

# form the tensor from the data
container <- form_tensor(container,
                         custom_genes = pagoda_hvgs,
                         donor_min_cells=25,
                         norm_method='trim', 
                         scale_factor=10000,
                         vargenes_method='norm_var_pvals', 
                         vargenes_thresh=.1,
                         scale_var = TRUE, 
                         var_scale_power = 2)

# Run tensor decomposition
container <- run_tucker_ica(container, ranks=c(6,8),
                            tucker_type = 'regular',
                            rotation_type = 'hybrid')

# I need to get the output

# 1) Explained variance factors

expl_var <- tibble(factor_label = paste0("Factor_", seq(1, length(container$exp_var))),
                   var = container$exp_var)

sum(expl_var$var)

write_csv(expl_var,"./results/MI/scITD/factor_totalR2.csv")

# 2) Explained variance per cell-type per factor

get_ctype_exp_var <- function(container, factor_use, ctype) {
  tnsr <- rTensor::as.tensor(container$tensor_data[[4]])
  donor_mat <- container$tucker_results[[1]]
  ldngs <- container$tucker_results[[2]]
  
  
  ctype_ndx <- which(container$tensor_data[[3]]==ctype)
  recon1 <- donor_mat[,factor_use,drop=FALSE] %*% ldngs[factor_use,,drop=FALSE]
  recon1 <- rTensor::k_fold(recon1,m=1,modes=tnsr@modes)
  
  # The reconstruction should be the original tensor with reconstruction for the one ctype
  recon2 <- tnsr
  recon2[,,ctype_ndx] <- recon2[,,ctype_ndx] - recon1[,,ctype_ndx]
  
  
  # calculate error from using just a single factor
  unexp_var <- (rTensor::fnorm(recon2)**2) / (rTensor::fnorm(tnsr)**2)
  exp_var <- (1 - unexp_var) * 100
  
  return(exp_var)
}

cts <- c("Fib", "CM", "Endo", "Myeloid", "PC", "vSMCs", "Lymphoid") %>%
  set_names()

expl_var_ct_fact <- map(cts, function(ct) {
  
  map_dbl(set_names(seq(1,6,1)), function(i){
    
    get_ctype_exp_var(container, i, ct)
    
  }) %>% 
    enframe(name = "feature") %>%
    dplyr::mutate(feature = paste0("Factor_", feature))
  
}) %>%
  enframe(name = "cell_type") %>%
  unnest() %>%
  pivot_wider(names_from = cell_type, values_from = value)

write_csv(expl_var_ct_fact,"./results/MI/scITD/factor_ctR2.csv")

# 3) Get the factor scores per sample

factor_scores <- container[["tucker_results"]][[1]] %>%
  as.data.frame()

colnames(factor_scores) <- paste0("Factor_", seq(1, ncol(factor_scores)))

write_csv(factor_scores %>%
            rownames_to_column("sample_id"), "./results/MI/scITD/factor_scores.csv")

factor_scores_long <- factor_scores %>% rownames_to_column("patient_region_id") %>%
  pivot_longer(-patient_region_id, names_to = "feature", values_to = "value")


# 4) Get loadings of all factors and all genes

loadings <- map(seq(1,6), ~ get_one_factor(container, factor_select=.x)[[2]] %>%
                  as.data.frame() %>%
                  rownames_to_column("gene") %>%
                  pivot_longer(-gene,names_to = "ctype")) %>%
  enframe() %>%
  mutate(name = paste0("Factor_", name)) %>%
  unnest()

write_csv(loadings, "./results/MI/scITD/loadings.csv")

# Generate associations of factors to covariates

meta <- sample_dict_red

get_associated_factors <- function(meta, factor_scores, predicted_label) {
  
  expl_var <- factor_scores %>%
    left_join(meta) %>%
    dplyr::select_at(c("patient_region_id", predicted_label, "feature", "value")) %>%
    group_by(feature) %>%
    nest() %>%
    mutate(pvalue = map(data, function(dat) {
      
      gene_aov <- aov(as.formula(paste0("value ~ ", predicted_label)), data = dat) %>%
        broom::tidy() %>%
        dplyr::filter(term == predicted_label) %>%
        dplyr::select(term, p.value)
      
      return(gene_aov)
    })) %>%
    tidyr::unnest(pvalue) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(adj_pvalue = p.adjust(p.value)) %>%
    dplyr::select(feature, term, p.value, adj_pvalue)
  
  return(expl_var)
  
}

# First with patient group

expl_var_pgroup <- get_associated_factors(meta = meta, 
                       factor_scores = factor_scores_long, 
                       predicted_label = "patient_group")

biol_factors <- expl_var_pgroup %>%
  dplyr::filter(adj_pvalue < 0.05) %>%
  pull(feature)

expl_var %>%
  dplyr::filter(factor_label %in% biol_factors) %>%
  pull(var) %>%
  mean()

factor_scores_long %>%
  left_join(meta) %>%
  dplyr::filter(feature %in% c("Factor_1", "Factor_3")) %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  ggplot(., aes(x = Factor_1,
                y = Factor_3,
                color = patient_group,
                shape = batch)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size =12))

# Then batch

expl_var_batch <- get_associated_factors(meta = meta, 
                       factor_scores = factor_scores_long, 
                       predicted_label = "batch")

batch_factors <- expl_var_batch %>%
  dplyr::filter(adj_pvalue < 0.05) %>%
  pull(feature)

expl_var %>%
  dplyr::filter(factor_label %in% batch_factors) %>%
  pull(var) %>%
  mean()

factor_scores_long %>%
  left_join(meta) %>%
  dplyr::filter(feature %in% c("Factor_6")) %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  ggplot(., aes(x = batch,
                y = Factor_6,
                color = batch)) +
  geom_boxplot() +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size =12),
        legend.position = "none") +
  xlab("") +
  scale_color_manual(values = set_names(c("black", "darkgrey"), c("A","B")))


# Plot heatmap of scores --------------------------------------------------------


ht_opt$ROW_ANNO_PADDING <- unit(2.5, "mm")
ht_opt$COLUMN_ANNO_PADDING <- unit(2.5, "mm")

# Association colors
col_fun_assoc <-  circlize::colorRamp2(seq(0, 20, length = 20), hcl.colors(20,"Purples",rev = T))
col_fun_r2 <- circlize::colorRamp2(seq(0, 20, length = 50), hcl.colors(50,"Oranges",rev = T))
col_fun_fact <- circlize::colorRamp2(seq(-1, 1, length = 50), hcl.colors(50,"Green-Brown",rev = T))


row_anns <- factor_scores %>% 
  rownames() %>% 
  enframe(value = "sample") %>% 
  dplyr::select(sample) %>% 
  left_join(meta, by = c("sample" = "patient_region_id")) %>%
  dplyr::select(patient_group, batch) %>%
  dplyr::mutate(patient_group = gsub("-enriched", "", patient_group)) %>%
  dplyr::mutate(patient_group = factor(patient_group,
                                       levels = c("myogenic", "fibrotic", "ischemic")))

row_ha <- rowAnnotation(condition = row_anns$patient_group,
                       batch = row_anns$batch,
                       col = list(condition = c("fibrotic" = "darkblue",
                                                   "myogenic" = "darkred", 
                                                   "ischemic" = "darkgreen"),
                                  batch = c("A" = "black",
                                            "B" = "darkgrey")),
                       gap = unit(2.5, "mm"),
                       border = TRUE)

# Add explain variance per cell-type
r2_per_factor <- expl_var_ct_fact %>%
  column_to_rownames("feature") %>%
  as.matrix()

assoc_pvals <- bind_rows("condition" = -log10(expl_var_pgroup$adj_pvalue),
                         "batch" = -log10(expl_var_batch$adj_pvalue)) %>%
  as.matrix()

column_ha <- HeatmapAnnotation("r2" = r2_per_factor,
                              "assocs" = assoc_pvals,
                              gap = unit(2.5, "mm"),
                              border = TRUE,
                              col = list(r2 = col_fun_r2,
                                         assocs = col_fun_assoc))


scores_hmap <- Heatmap(factor_scores, 
                       name = "factor_scores", 
                       right_annotation = row_ha,
                       top_annotation = column_ha,
                       cluster_columns = FALSE, 
                       show_row_dend = TRUE,
                       show_row_names = FALSE,
                       border = TRUE,
                       gap = unit(2.5, "mm"),
                       col = col_fun_fact)

pdf("./results/MI/scITD/factor_scores.pdf", height = 5, width = 3.7)

scores_hmap

dev.off()



