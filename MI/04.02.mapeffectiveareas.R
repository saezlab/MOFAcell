# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' I will map corrected states to visium slides
#' and correlate them in space

library(Seurat)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)

# Get individual slide info ---------------------------------------------
visium_folder <- "/Users/ricardoramirez/Dropbox/PhD/Research/mi_atlas/processed_visium/objects/"
visium_files <- list.files(visium_folder, full.names = F)
visium_samples <- gsub("[.]rds", "", visium_files)

visium_df <- tibble(visium_file = paste0(visium_folder, 
                                         visium_files),
                    slide = visium_samples)

# Meta information of slides
annotation_names <- tibble(patient_group = c("group_1", 
                                             "group_2", 
                                             "group_3"),
                           patient_group_name = factor(c("myogenic-enriched", 
                                                         "ischemic-enriched", 
                                                         "fibrotic-enriched"),
                                                       levels = c("myogenic-enriched", 
                                                                  "ischemic-enriched", 
                                                                  "fibrotic-enriched")))

meta <- read_csv("./data_MI/visium_patient_anns_revisions.csv") %>%
  left_join(annotation_names)


# Get module info -------------------------------------------------------
module_folder <- "./results/MI/MOFA_mcell/factor_desc/Factor1_char/decoupler_ct/"

param_df <- list.files(module_folder) %>%
  enframe(value = "file") %>%
  dplyr::select(-name) %>%
  dplyr::filter(grepl("msk", file)) %>%
  dplyr::mutate(slide = gsub("_msk.csv", "", file)) %>%
  left_join(visium_df, by = "slide") %>%
  dplyr::mutate(file = paste0(module_folder,file))



# We will plot in slides
# the scores
extract_visium_scores <- function(file, slide, visium_file) {
  
  print(slide)
  
  scores <- read_csv(file) %>%
    dplyr::mutate(spot_id = gsub(slide, "", spot_id) %>%
                    strsplit(., "_") %>%
                    map_chr(., ~.x[[2]])) %>%
    column_to_rownames("spot_id") %>%
    as.matrix()
  
  visium_slide <- readRDS(visium_file)
  
  scores <- scores[colnames(visium_slide),] %>%
    t()
  
  visium_slide[["MOFA_mcell"]] <- Seurat::CreateAssayObject(scores)
  
  DefaultAssay(visium_slide) <- "MOFA_mcell"
  
  pdf(paste0("./results/MOFA_mcell/visium_plts/", slide, ".pdf"), height = 4, width = 6)
  
  for(i in 1:length(rownames(visium_slide))) {
    
  slide_plt <- SpatialFeaturePlot(visium_slide, 
                                  features = rownames(visium_slide)[i],
                                  stroke = 0) +
      viridis::scale_fill_viridis(option = "B")
    
    plot(slide_plt)
    
  }
  
  dev.off()
  
  return(NULL)
}

# Correlate scores:
# Whenever is 0 in both places, don't correlate

extract_correlation <- function(file, slide) {
  
  scores <- read_csv(file) %>%
    dplyr::mutate(spot_id = gsub(slide, "", spot_id) %>%
                    strsplit(., "_") %>%
                    map_chr(., ~.x[[2]])) %>%
    column_to_rownames("spot_id") %>%
    as.matrix()
  
  module_comb <- combn(colnames(scores),2) %>%
    t() %>%
    as.data.frame() 
  
  cors <- map2(module_comb$V1, module_comb$V2, function(x,y) {
    
    mat <- scores[,c(x,y)]
    ix <- which(rowSums(mat)!=0)
    mat <- mat[ix,]
    
    cor_res <- cor(mat) %>%
      as.data.frame() %>%
      rownames_to_column("moduleA") %>%
      pivot_longer(-moduleA, names_to = "moduleB",values_to = "corr")
    
    return(cor_res)
    
  }) %>%
    bind_rows()
  
  cors_mat <- cors %>%
    dplyr::filter(moduleA != moduleB) %>%
    pivot_wider(names_from = moduleB, 
                values_from = corr,
                values_fill = 1) %>%
    column_to_rownames("moduleA") %>%
    as.matrix()
  
  #pdf(paste0("./results/MOFA_mcell/cor_plts/", slide, ".pdf"), height = 5, width = 5)
  
  #draw(ComplexHeatmap::Heatmap(cors_mat,name = "cor"))
  
  #dev.off()
  
  return(cors %>%
           dplyr::filter(moduleA != moduleB))
}


# Main ----------------------------------------------------------------------------------------

# Plot scores
pmap(param_df, extract_visium_scores)

# Correlate scores
slide_cors <- pmap(param_df[,c("file", "slide")],
     extract_correlation)

names(slide_cors) <- param_df$slide

slide_cors <- enframe(slide_cors) %>% unnest()

slide_cors <- slide_cors %>% na.omit()

aov_res <- slide_cors %>%
  left_join(meta, by = c("name" = "sample_id")) %>%
  group_by(moduleA, moduleB) %>%
  nest() %>%
  dplyr::mutate(anova_ints = map(data, function(dat){
    
    aov(corr ~ patient_group_name,data = dat) %>%
      broom::tidy() %>%
      dplyr::filter(term == "patient_group_name")
    
    
  })) %>%
  unnest(anova_ints) %>%
  dplyr::mutate(mean_corrs = map(data, function(dat){
    
    dat %>%
      group_by(patient_group_name) %>%
      summarize(mean_cor = mean(corr)) %>%
      pivot_wider(names_from = patient_group_name,
                  values_from = mean_cor)
    
    
  })) %>%
  unnest(mean_corrs) %>%
  dplyr::mutate(moduleA_ct = strsplit(moduleA, "_") %>%
                  map_chr(., ~.x[[1]]),
                moduleA_sign = strsplit(moduleA, "_") %>%
                  map_chr(., ~.x[[2]]),
                moduleB_ct = strsplit(moduleB, "_") %>%
                  map_chr(., ~.x[[1]]),
                moduleB_sign = strsplit(moduleB, "_") %>%
                  map_chr(., ~.x[[2]]))

# Within the same cell-type, correspondant programs 
# are usually anticorrelated in space
aov_res %>%
  dplyr::filter(moduleA_ct == moduleB_ct &
                  moduleA_sign == "pos" &
                  moduleB_sign == "neg")

# Between the cell-types, correspondant programs 
# are usually anticorrelated in space
aov_res %>%
  dplyr::filter(moduleA_ct != moduleB_ct &
                  moduleA_sign == "pos" &
                  moduleB_sign == "neg") %>%
  print(n = 50)

# Between the cell-types, same direction programs 
# are usually anticorrelated in space
aov_res %>%
  dplyr::filter(moduleA_ct != moduleB_ct &
                  moduleA_sign == "pos" &
                  moduleB_sign == "pos") %>%
  print(n = 50)

aov_res %>%
  dplyr::select(moduleA_ct,moduleB_ct,moduleA_sign,moduleB_sign, data) %>%
  unnest() %>%
  dplyr::filter(moduleA_ct != moduleB_ct &
                  moduleA_sign == "neg" &
                  moduleB_sign == "neg") %>%
  print(n = 50) %>%
  arrange(-corr)



aov_res %>%
  dplyr::select(moduleA_ct,moduleB_ct,moduleA_sign,moduleB_sign, data) %>%
  unnest() %>%
  dplyr::filter(moduleA_ct == "Endo" &
                  moduleA_sign == "pos" &
                  moduleB_sign == "pos",
                patient_group == "group_2") %>%
  arrange(-corr)


aov_res <- aov_res %>%
  dplyr::filter(p.value < 0.05) %>%
  arrange(p.value)





