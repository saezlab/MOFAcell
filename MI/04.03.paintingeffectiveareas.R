# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' I will map corrected states to visium slides
#' and correlate them in space

library(Seurat)
library(tidyverse)
source("./MOFAcell/code/paintR.R")

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
  
  # Here we make the scores RGB friendly
  
  scores <- sweep(scores, MARGIN = 2, STATS = apply(scores, 2, max), FUN = "/")
  
  scores[is.nan(scores)] <- 0
  
  visium_slide <- readRDS(visium_file)
  
  geometry <- GetTissueCoordinates(visium_slide) %>%
    as.data.frame() %>%
    rownames_to_column("spot_id")
  
  return(list("geometry" = geometry, "scores" = scores))

}

score_list <- pmap(param_df, extract_visium_scores)
names(score_list) <- param_df$slide

# Trying paintR

# For a given cell-type paint the healthy and disease areas

paintCTs <- function(score_list_obj, ct, col_select = c("red", "blue")) {
  
  ct_names <- colnames(score_list_obj$scores)
  ct_select <- ct_names[grepl(ct, ct_names)]
  comb_comps_cols <-  matrix_to_hex(spot_matrix = score_list_obj$scores[,ct_select] %>%
                                      t(), colors = col_select)
  
  coordinate_data <- score_list_obj$geometry %>%
    left_join(comb_comps_cols$data_hex, by = c("spot_id" = "sample_id"))
  
  
  slide_plt <- ggplot(coordinate_data, aes(x = imagecol, y = imagerow * -1, color = spot_id)) +
    geom_point() +
    scale_color_manual(values = set_names(coordinate_data$hex, coordinate_data$spot_id)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title =element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("") + ylab("")
  
  #Legend
  legend <- ggplot(comb_comps_cols$legend_key, aes(y = name, x = 1, fill = name)) +
    geom_tile() +
    scale_fill_manual(values = set_names(comb_comps_cols$legend_key$hex, comb_comps_cols$legend_key$name)) +
    ylab("") +
    xlab("") +
    coord_equal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("")
  
  f_plot <- cowplot::plot_grid(slide_plt, legend, rel_widths = c(1,0.4))
  
  return(f_plot)
}

# Paint for CMs

fibrotic_CM <- paintCTs(score_list_obj = score_list$AKK001_157785,
                        ct = "CM",col_select = c("red","blue"))

pdf("./results/MI/paintR/CM_fibrotic.pdf", height = 4, width = 5.5)

plot(fibrotic_CM)

dev.off()

healthy_CM <- paintCTs(score_list_obj = score_list$AKK006_157771,
                        ct = "CM",col_select = c("red","blue"))

pdf("./results/MI/paintR/CM_healthy.pdf", height = 4, width = 5.5)

plot(healthy_CM)

dev.off()

ischemic_CM <- paintCTs(score_list_obj = score_list$Visium_19_CK297,
                       ct = "CM",col_select = c("red","blue"))

pdf("./results/MI/paintR/CM_ischemic.pdf", height = 4, width = 5.5)

plot(ischemic_CM)

dev.off()

# Paint for Fibs

fibrotic_Fib <- paintCTs(score_list_obj = score_list$AKK001_157785,
                        ct = "Fib",col_select = c("red","blue"))

pdf("./results/MI/paintR/Fib_fibrotic.pdf", height = 4, width = 5.5)

plot(fibrotic_Fib)

dev.off()

healthy_Fib <- paintCTs(score_list_obj = score_list$AKK006_157771,
                       ct = "Fib",col_select = c("red","blue"))

pdf("./results/MI/paintR/Fib_healthy.pdf", height = 4, width = 5.5)

plot(healthy_Fib)

dev.off()

ischemic_Fib <- paintCTs(score_list_obj = score_list$Visium_19_CK297,
                        ct = "Fib",col_select = c("red","blue"))

pdf("./results/MI/paintR/Fib_ischemic.pdf", height = 4, width = 5.5)

plot(ischemic_Fib)

dev.off()
