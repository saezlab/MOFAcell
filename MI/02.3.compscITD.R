# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we will compare the factors
#' and loadings generated from MOFA and scITD

library(MOFA2)
library(tidyverse)
library(ComplexHeatmap)
library(cluster)
source("./MOFAcell/factorutils.R")

# First, identify and subset factors of interest from the MOFA model

model_outfile <- "./results/MI/MOFA_mcell/MI_model6factors.hdf5"
model <- MOFA2::load_model(model_outfile)
pb_data <- readRDS("./data_MI/mi_pb_red.rds")
out_alias <- "./results/MI/MOFA_mcell/factor_desc"

# Meta data

annotation_names <- tibble(patient_group = c("group_1", 
                                             "group_2", 
                                             "group_3"),
                           patient_group_name = c("myogenic", 
                                                  "ischemic", 
                                                  "fibrotic"))

batch_info <- read_csv("./data_MI/snrna_batch_ann.csv") %>%
  select(orig.ident, batch) %>%
  unique()

sample_dict <- read_csv("./data_MI/rna_patient_anns_revisions.csv",show_col_types = F) %>%
  left_join(annotation_names, by = "patient_group") %>%
  dplyr::select(-patient_group) %>%
  dplyr::rename("patient_group" = patient_group_name) %>%
  dplyr::mutate(patient_group = factor(patient_group, 
                                       levels = c("myogenic", "ischemic", "fibrotic"))) %>%
  #dplyr::select(-sample_id) %>%
  unique()

sample_dict_red <- sample_dict %>%
  left_join(batch_info, by = c("sample_id" = "orig.ident")) %>%
  dplyr::select(patient_region_id, major_labl, patient_group, batch) %>%
  unique() %>%
  dplyr::rename("sample" = patient_region_id)

meta <- sample_dict_red

# Get factor scorees
factor_scores <- get_fscores(model = model,
                             meta = meta) %>%
  dplyr::select(-c(major_labl, patient_group, batch)) %>%
  column_to_rownames("sample")

colnames(factor_scores) <- paste0("MOFA", "_", 
                                        colnames(factor_scores))

# Get loadings
factor_loadings <-  get_weights(model, as.data.frame = T) %>%
  as.data.frame() %>%
  dplyr::mutate(feature = strsplit(as.character(feature), "_") %>%
                  map_chr(., ~ .x[[2]]),
                ctype = strsplit(as.character(view), "_") %>%
                  map_chr(., ~ .x[[1]])) %>%
  dplyr::rename("factors" = factor) %>%
  dplyr::select(-view) %>%
  dplyr::group_by(factors) %>%
  nest() %>%
  mutate(data = map(data, function(dat) {
    dat %>%
      pivot_wider(names_from = ctype, values_from = value,values_fill = 0) %>%
      column_to_rownames("feature") %>%
      as.matrix()
  }))

# scITD -----------------------------------------------------

scITD_factor_scores <- read_csv("./results/MI/scITD/factor_scores.csv",
                                show_col_types = F) %>%
  #dplyr::select_at(c("sample_id","Factor_2","Factor_3", "Factor_4")) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

colnames(scITD_factor_scores) <- paste0("scITD", "_", 
                                        colnames(scITD_factor_scores))

scITD_loadings <- read_csv("./results/MI/scITD/loadings.csv")  %>%
  as.data.frame() %>%
  #dplyr::filter(name %in% c("Factor_2","Factor_3", "Factor_4")) %>%
  dplyr::group_by(name) %>%
  nest() %>%
  mutate(data = map(data, function(dat) {
    dat %>%
      pivot_wider(names_from = ctype, values_from = value,values_fill = 0) %>%
      column_to_rownames("gene") %>%
      as.matrix()
  }))

# Added sillouhette width comparison

calculate_sw <- function(scores, meta, test_label) {
  
  meta_info <- meta %>%
    dplyr::select_at(c("sample", test_label)) %>%
    column_to_rownames("sample") %>%
    dplyr::mutate(label_ix = .data[[test_label]] %>% 
                    as.factor() %>% 
                    as.integer())
  
  meta_info <- meta_info[rownames(scores), ]
  
  patient_dists <- dist(scores,
                     method = "euclidean")
  
  si <- (silhouette(x = meta_info$label_ix, patient_dists)) %>%
    as.data.frame() %>%
    dplyr::select(cluster, sil_width) %>%
    left_join(unique(meta_info), by = c("cluster" = "label_ix"))

  return(si)
}

# Compare batch clustering

batch_sw <- bind_rows(calculate_sw(factor_scores,meta = meta, "batch") %>%
  mutate(method = "MOFA"),
calculate_sw(scITD_factor_scores,meta = meta, "batch") %>%
  mutate(method = "scITD"))

batch_sw %>% group_by(batch) %>%
  nest() %>%
  dplyr::mutate(comp_sw = map(data, function(dat) {
    
    t.test(sil_width ~ method, dat) %>%
      broom::tidy() %>%
      dplyr::select(p.value)
    
  })) %>%
  dplyr::select(batch, comp_sw) %>%
  unnest()

batch_sw_plt <- ggplot(batch_sw,
       aes(x = batch, y = sil_width, color = method)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(0.5,0,0,0.5), "cm")) + #
  scale_color_manual(values = c("black", "grey")) +
  xlab("") +
  ylab("silhouette \n width")

pdf("./results/MI/comparison/sw_batch.pdf", height = 2, width = 3.5)
plot(batch_sw_plt)
dev.off()


# Compare patient clustering

patient_sw <- bind_rows(calculate_sw(factor_scores,
                                     meta = meta, "patient_group") %>%
                        mutate(method = "MOFA"),
                      calculate_sw(scITD_factor_scores,meta = meta, "patient_group") %>%
                        mutate(method = "scITD"))

patient_sw_plt <- ggplot(patient_sw,
       aes(x = patient_group, y = sil_width, color = method)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
  theme(plot.margin = unit(c(0.5,0,0,0.5), "cm")) + #
  scale_color_manual(values = c("black", "grey")) +
  xlab("") +
  ylab("silhouette \n width")

pdf("./results/MI/comparison/sw_patient.pdf", height = 2.5, width = 3.5)
plot(patient_sw_plt)
dev.off()

patient_sw %>% group_by(patient_group) %>%
  nest() %>%
  dplyr::mutate(comp_sw = map(data, function(dat) {
    
    t.test(sil_width ~ method, dat) %>%
      broom::tidy() %>%
      dplyr::select(p.value)
    
  })) %>%
  dplyr::select(patient_group, comp_sw) %>%
  unnest() 

# Paper panel
pdf("./results/MI/comparison/sw_patientbatch.pdf", height = 2.5, width = 6.5)

cowplot::plot_grid(patient_sw_plt, batch_sw_plt, align = "hv", rel_widths = c(0.9,0.75))

dev.off()
# Comparisons

pat_order <- rownames(scITD_factor_scores)

factor_cor <- cor(factor_scores[pat_order,],
    scITD_factor_scores[pat_order,], method = "spearman")
  

pdf(height = 5, width = 5, "./results/MI/comparison/corr_factors.pdf")

corrplot::corrplot(factor_cor, method = "color", 
                   addCoef.col = 'black', tl.col = 'black')

dev.off()


# Comparison of loadings

comp_loadings <- function(MOFA_factor, scITD_factor) {
  
  scITDdat <- scITD_loadings %>%
    dplyr::filter(name == scITD_factor) %>%
    pull(data)
  
  scITDdat <- scITDdat[[1]]
  
  colnames(scITDdat) <- paste0("scITD", "_", 
                               colnames(scITDdat))
  
   MOFAdat <- factor_loadings %>%
    dplyr::filter(factors == MOFA_factor) %>%
    pull(data)
   
   MOFAdat <- MOFAdat[[1]]
  
   colnames(MOFAdat) <- paste0("MOFA", "_", 
                                colnames(MOFAdat))
   
   shared_genes <- intersect(rownames(MOFAdat), 
                             rownames(scITDdat))
   
   corres <- cor(MOFAdat[shared_genes, ], 
              scITDdat[shared_genes,], method = "spearman")
   
   return(corres * diag(7))
  

}

# MOFA 1 vs scITD 1,3,4

comp_loadings("Factor1", "Factor_1") %>%
  corrplot::corrplot(method = "color", addCoef.col = 'black', tl.col = 'black')

diag(comp_loadings("Factor1", "Factor_1") ) %>% mean()
diag(comp_loadings("Factor1", "Factor_3") ) %>% mean()
diag(comp_loadings("Factor1", "Factor_4") ) %>% mean()


load_corr <- set_names(diag(comp_loadings("Factor1", "Factor_1")),
                       rownames(comp_loadings("Factor1", "Factor_1")) 
                       %>% strsplit(.,"_") 
                       %>% map_chr(., ~.x %>% last())) %>%
  enframe() %>%
  arrange(-value)

pdf(height = 3.5, width = 2.5, "./results/MI/comparison/corr_loadings.pdf")

ggplot(load_corr, aes(x = factor(name,
                                 levels = .data[["name"]]), y = value)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust =0.5),
        axis.text.y = element_text(size =12)) +
  xlab("") +
  ylab("Spearman Correlation \n scITD Fact. 1 and MOFA Fact. 1")

dev.off()

# Checking how scITD added background noise

scITD_fscores_plt <- scITD_factor_scores %>% 
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(meta) %>%
  dplyr::mutate(patient_group = gsub("-enriched", "", patient_group)) %>%
  dplyr::mutate(patient_group = factor(patient_group,
                                       levels = c("myogenic", "fibrotic", "ischemic")))  %>%
  pivot_longer(-c(sample, major_labl, patient_group, batch)) %>%
  dplyr::filter(name %in% c("scITD_Factor_1",
                            "scITD_Factor_3",
                            "scITD_Factor_4")) %>%
  ggplot(aes(x = patient_group, y = value, color =  patient_group)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(.~ name, ncol = 3) +
  theme_classic() +
  xlab("") +
  ylab("Score") +
  scale_color_manual(values = set_names(c("darkred", "darkgreen", "darkblue"), 
                                        c("myogenic","ischemic", "fibrotic"))) +
  theme(axis.text.x = element_text(size =12, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size =12),
        legend.position = "none")

pdf("./results/MI/scITD/condition_scores.pdf", height = 2.3, width = 4.0)
plot(scITD_fscores_plt)
dev.off()

scITDdat <- scITD_loadings %>%
  dplyr::filter(name == "Factor_1") %>%
  pull(data)

MOFAdat <- factor_loadings %>%
  dplyr::filter(factors == "Factor1") %>%
  pull(data)

background_plt <- cbind("scITD" = scITDdat[[1]]["POSTN",], 
      "MOFA" = MOFAdat[[1]]["POSTN",]) %>%
  data.frame() %>%
  rownames_to_column("cell_type") %>%
  ggplot(aes(x = MOFA, y = scITD, color = cell_type)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm")) + 
  ggtitle("POSTN") +
  xlab("MOFA Factor 1") +
  ylab("scITD Factor 1")

pdf("./results/MI/comparison/background_comparison_POSTN.pdf", height = 2.5, width = 3.5)
plot(background_plt)
dev.off()


background_plt <- cbind("scITD" = scITDdat[[1]]["TTN",], 
                        "MOFA" = MOFAdat[[1]]["TTN",]) %>%
  data.frame() %>%
  rownames_to_column("cell_type") %>%
  ggplot(aes(x = MOFA, y = scITD, color = cell_type)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  theme(plot.margin = unit(c(0.5,0,0.5,0.5), "cm")) + 
  ggtitle("TTN") +
  xlab("MOFA Factor 1") +
  ylab("scITD Factor 1")

pdf("./results/MI/comparison/background_comparison_TTN.pdf", height = 2.5, width = 3.5)
plot(background_plt)
dev.off()

