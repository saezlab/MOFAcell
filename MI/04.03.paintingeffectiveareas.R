# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' I will map corrected states to visium slides
#' and correlate them in space

library(Seurat)
library(tidyverse)
source("./MOFAcell/paintR.R")

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
  dplyr::filter(grepl("msk.csv", file)) %>%
  dplyr::mutate(slide = gsub("_msk.csv", "", file)) %>%
  left_join(visium_df, by = "slide") %>%
  dplyr::mutate(file = paste0(module_folder,file))

# Here we will define a unified scoring across slides -------------------

extract_max_vals <- function(file_path) {
  
  scores <- read_csv(file_path, show_col_types = F) %>%
    pivot_longer(-spot_id) %>%
    dplyr::select(-spot_id) %>%
    #dplyr::mutate(name = strsplit(name, "_") %>%
    #              map_chr(., ~.x[[1]])) %>%
    group_by(name) %>%
    dplyr::summarise(max_val = max(value))
  
  return(scores)
  
}

norm_fact_df <- param_df %>%
  dplyr::mutate(max_scores = map(file, extract_max_vals)) %>%
  dplyr::select(slide, max_scores) %>%
  unnest() %>%
  group_by(name) %>%
  summarize(norm_fact = max(max_val))
  
norm_fact <- set_names(norm_fact_df$norm_fact,
                       norm_fact_df$name)

# We will plot in slides
# the scores
extract_visium_scores <- function(file, slide, visium_file) {
  
  print(slide)
  
  scores <- read_csv(file, show_col_types = F) %>%
    dplyr::mutate(spot_id = gsub(slide, "", spot_id) %>%
                    strsplit(., "_") %>%
                    map_chr(., ~.x[[2]])) %>%
    column_to_rownames("spot_id") %>%
    as.matrix()
  
  ix <- colnames(scores) #%>%
    #strsplit(., "_") %>%
    #map_chr(., ~.x[[1]])
  
  
  # Here we make the scores RGB friendly
  #scores <- sweep(scores, MARGIN = 2, STATS = apply(scores, 2, max), FUN = "/")
  
  scores <- sweep(scores, MARGIN = 2, STATS = norm_fact[ix], FUN = "/")
  
  scores[is.nan(scores)] <- 0
  
  visium_slide <- readRDS(visium_file)
  
  geometry <- GetTissueCoordinates(visium_slide) %>%
    as.data.frame() %>%
    rownames_to_column("spot_id")
  
  return(list("geometry" = geometry, "scores" = scores))

}

score_list <- param_df %>%
  dplyr::filter(slide %in% c("Visium_1_CK279", "Visium_14_CK292" ,"Visium_15_CK293")) %>%
  pmap(., extract_visium_scores)

names(score_list) <- param_df %>%
  dplyr::filter(slide %in% c("Visium_1_CK279", "Visium_14_CK292" ,"Visium_15_CK293")) %>%
  pull(slide)


# Trying paintR

# For a given cell-type paint the healthy and disease areas

paintCTs <- function(score_list_obj, ct, col_select = c("blue", "red")) {
  
  ct_names <- colnames(score_list_obj$scores)
  ct_select <- ct_names[grepl(ct, ct_names)] %>%
    sort() # so as to have neg always first
  
  comb_comps_cols <-  matrix_to_hex(spot_matrix = score_list_obj$scores[,ct_select] %>%
                                      t(), colors = col_select)
  
  coordinate_data <- score_list_obj$geometry %>%
    left_join(comb_comps_cols$data_hex, by = c("spot_id" = "sample_id"))
  
  
  slide_plt <- ggplot(coordinate_data, aes(x = imagecol, y = imagerow * -1, color = spot_id)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = set_names(coordinate_data$hex, coordinate_data$spot_id)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title =element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black",
                                          colour = "black")) +
    xlab("") + ylab("")
  
  return(slide_plt)
}

# Paint for CMs
cts <- c("CM","Endo",  "Fib", "Myeloid") %>%
  set_names() 

myogenic_column <- map(cts, paintCTs, 
                       score_list_obj = score_list$Visium_1_CK279,
                       col_select = c("blue", "red")) %>%
  cowplot::plot_grid(plotlist = ., ncol = 1, align = "hv")

fibrotic_column <- map(cts, paintCTs, 
                       score_list_obj = score_list$Visium_14_CK292,
                       col_select = c("blue", "red")) %>%
  cowplot::plot_grid(plotlist = ., ncol = 1, align = "hv")

ischemic_column <- map(cts, paintCTs, 
                       score_list_obj = score_list$Visium_15_CK293,
                       col_select = c("blue", "red")) %>%
  cowplot::plot_grid(plotlist = ., ncol = 1, align = "hv")

spatial_map <- cowplot::plot_grid(myogenic_column, ischemic_column, fibrotic_column, ncol = 3)

pdf("./results/MI/paintR/spatial_map.pdf", height = 6, width = 5)

plot(spatial_map)

dev.off()
