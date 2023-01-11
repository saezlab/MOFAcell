# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' I will map corrected states to visium slides
#' and correlate them in space, particularly disease states

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
    dplyr::mutate(name = strsplit(name, "_") %>%
                    map_chr(., ~.x[[1]])) %>%
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


# We will plot in slides the scores -------------------------
# Normalizing by the maximum per cell-type ------------------

extract_visium_scores <- function(file, slide, visium_file) {
  
  print(slide)
  
  scores <- read_csv(file, show_col_types = F) %>%
    dplyr::mutate(spot_id = gsub(slide, "", spot_id) %>%
                    strsplit(., "_") %>%
                    map_chr(., ~.x[[2]])) %>%
    column_to_rownames("spot_id") %>%
    as.matrix()
  
  # Here we make the scores RGB friendly
  
  ix <- colnames(scores) %>%
    strsplit(., "_") %>%
    map_chr(., ~.x[[1]])
  
  norm_fact_max = norm_fact
  norm_fact_max[names(norm_fact_max)] = max(norm_fact_max)

  #Recipe 2
  scores <- sweep(scores, MARGIN = 2, STATS = norm_fact_max[ix], FUN = "/")
  
  scores[is.nan(scores)] <- 0
  
  visium_slide <- readRDS(visium_file)
  
  geometry <- GetTissueCoordinates(visium_slide) %>%
    as.data.frame() %>%
    rownames_to_column("spot_id")
  
  return(list("geometry" = geometry, "scores" = scores))
  
}

score_list <- pmap(param_df, extract_visium_scores)
names(score_list) <- param_df$slide

# For a given factor, plot mcell relationships

paintCTs <- function(score_list_obj, 
                     sign_type,
                     cts,
                     col_select = NULL) {
  
  ct_names <- colnames(score_list_obj$scores)
  ct_select <- ct_names[grepl(sign_type, ct_names)]
  
  score_mat <- score_list_obj$scores[, ct_select, drop = F]
  colnames(score_mat) <- gsub(paste0("_",sign_type), "", colnames(score_mat))
  
  ct_names <- colnames(score_mat)
  ct_select <- ct_names[ct_names %in% cts]
  
  comb_comps_cols <-  matrix_to_hex(spot_matrix = score_mat[,ct_select] %>%
                                      t(), colors = col_select)
  
  coordinate_data <- score_list_obj$geometry %>%
    left_join(comb_comps_cols$data_hex, by = c("spot_id" = "sample_id"))
  
  
  slide_plt <- ggplot(coordinate_data, aes(x = imagecol, y = imagerow * -1, color = spot_id)) +
    geom_point(size = 2) +
    scale_color_manual(values = set_names(coordinate_data$hex, coordinate_data$spot_id)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title =element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.background = element_rect(fill = "black",
                                          colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          ) +
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
          axis.ticks.x=element_blank()) +
    xlab("")
  
  f_plot <- cowplot::plot_grid(slide_plt, legend, rel_widths = c(1,0.4))
  
  return(f_plot)
}

# Main: Plot examples in best in class slides

best_in_class <- read_csv("./results/MI/MOFA_mcell/factor_desc/Factor1_char/misty/bestinclass.csv")
best_in_class$sample %>% table() %>% sort(decreasing = T)
best_in_class$target %>% table() %>% sort()

# Fibrotic
fib_paint <- paintCTs(score_list_obj = score_list$Visium_13_CK291, 
                      cts = c("CM", "Fib", "Myeloid"),
                      sign_type = "pos", 
                      col_select = c("red", "blue", "green"))

pdf("./results/MI/paintR/fibrotic_paintR3.pdf", height = 4.5, width = 6.5)

plot(fib_paint)

dev.off()

# Ischemic
ischemic_paint <- paintCTs(score_list_obj = score_list$Visium_18_CK296, 
                           cts = c("CM", "Fib", "Myeloid"),
                           sign_type = "pos", 
                           col_select = c("red", "blue", "green"))

pdf("./results/MI/paintR/ischemic_paintR3.pdf", height = 4.5, width = 6.5)

plot(ischemic_paint)

dev.off()



myogenic_paint <- paintCTs(score_list_obj = score_list$Visium_1_CK279, 
                           cts = c("CM", "Fib", "Myeloid"),
                           sign_type = "pos", 
                           col_select = c("red", "blue", "green"))

pdf("./results/MI/paintR/myogenic_paintR3.pdf", height = 4.5, width = 6.5)

plot(myogenic_paint)

dev.off()


myogenic_paint <- paintCTs(score_list_obj = score_list$Visium_1_CK279, 
                           cts = c("CM", "Fib", "Myeloid"),
                           sign_type = "pos", 
                           col_select = c("red", "blue", "green"))

ischemic_paint <- paintCTs(score_list_obj = score_list$Visium_15_CK293, 
                           cts = c("CM", "Fib", "Myeloid"),
                           sign_type = "pos", 
                           col_select = c("red", "blue", "green"))

ischemic_paint <- paintCTs(score_list_obj = score_list$Visium_16_CK294, 
                           cts = c("CM", "Fib", "Myeloid"),
                           sign_type = "pos", 
                           col_select = c("red", "blue", "green"))

