# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Here we explore if the + programs of a given cell-type
#' depend on the + programs of other cell-types
#' in space

library(Seurat)
library(tidyverse)
library(mistyR)
source("./MOFAcell/misty_utilities.R")

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
misty_folder <- "./results/MI/MOFA_mcell/factor_desc/Factor1_char/misty/"

param_df <- list.files(module_folder) %>%
  enframe(value = "file") %>%
  dplyr::select(-name) %>%
  dplyr::filter(grepl("msk.csv", file)) %>%
  dplyr::mutate(slide = gsub("_msk.csv", "", file)) %>%
  left_join(visium_df, by = "slide") %>%
  dplyr::mutate(file = paste0(module_folder,file),
                misty_file = paste0(misty_folder, slide)) %>%
  left_join(meta, by = c("slide" = "sample_id")) %>%
  dplyr::filter(patient_group != "group_1")

# Get normalization scores

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
 
# Fit a MISTy model with a masked cell loading matrix --------------------

multicell_framework <- function(file, slide, visium_file, misty_file) {
  
  print(slide)
  
  scores <- read_csv(file) %>%
    dplyr::mutate(spot_id = gsub(slide, "", spot_id) %>%
                    strsplit(., "_") %>%
                    map_chr(., ~.x[[2]])) %>%
    column_to_rownames("spot_id") %>%
    as.matrix()
  
  ix <- colnames(scores) %>%
    strsplit(., "_") %>%
    map_chr(., ~.x[[1]])
  
  # Here we make the scores RGB friendly
  
  norm_fact_max = norm_fact
  norm_fact_max[names(norm_fact_max)] = max(norm_fact_max)
  
  #Recipe 2
  scores <- sweep(scores, MARGIN = 2, STATS = norm_fact_max[ix], FUN = "/")
  
  scores[is.nan(scores)] <- 0
  
  # Here we also make the scores, Seurat ready
  
  visium_slide <- readRDS(visium_file)
  
  scores <- scores[colnames(visium_slide),]
  
  scores <- scores[,grepl("pos", colnames(scores))]
  
  visium_slide[['ct_loads']] = CreateAssayObject(data = t(scores))
  
  # Create MISTy pipeline
  
  #Defining constant parameter of the pipeline
  view_assays <- list("main" = "ct_loads",
                      "para_ct" = "ct_loads")
  
  view_features <- list("main" = NULL, #uses all
                        "para_ct" = NULL)
  
  view_types <- list("main" = "intra", #uses all
                     "para_ct" = "para")
  
  view_params = list("main" = NULL,
                     "para_ct" = 5)
  
  misty_out <- misty_file
  
  run_misty_seurat(visium.slide = visium_slide,
                   # Seurat object with spatial transcriptomics data.
                   view.assays = view_assays,
                   # Named list of assays for each view.
                   view.features = view_features,
                   # Named list of features/markers to use.
                   # Use all by default.
                   view.types = view_types,
                   # Named list of the type of view to construct
                   # from the assay.
                   view.params = view_params,
                   # Named list with parameters (NULL or value)
                   # for each view.
                   spot.ids = NULL,
                   # spot IDs to use. Use all by default.
                   out.alias = misty_out,
                   bypass.intra = FALSE)
  # folder name for output
  
}

# Run MISTy models

fast_param_df = param_df %>%
  dplyr::select(file, slide, visium_file, misty_file) %>%
  dplyr::filter(slide != "Visium_13_CK291")

pmap(fast_param_df, multicell_framework)

# Generate plots of the interactions

all_misty <- mistyR::collect_results(param_df %>% 
                                       dplyr::filter(patient_group == "group_3",
                                                     slide != "Visium_13_CK291") %>% 
                                       pull(misty_file))

all_misty$improvements.stats

mistyR::plot_improvement_stats(all_misty, "multi.R2")

all_misty$contributions %>% dplyr::filter(view == "para_ct_5") %>% group_by(target) %>% summarize(median(value))

all_misty$improvements %>% dplyr::filter(measure == "multi.R2") %>% group_by(target) %>% summarize(median(value))


R2_plt <- last_plot() +
  theme(axis.text.x = element_text(size = 12, angle = 90, 
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size =12))

pdf("./results/MI/MOFA_mcell/factor_desc/Factor1_char/misty/R2_metric.pdf", height = 3, width = 2.3)

plot(R2_plt)

dev.off()

mistyR::plot_interaction_heatmap(all_misty, "intra",cutoff = 0)

intra_plt <- last_plot() +
  theme(axis.text.x = element_text(size = 12, angle = 90, 
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size =12)) +
  ggtitle("colocalization")

pdf("./results/MI/MOFA_mcell/factor_desc/Factor1_char/misty/intra_imp.pdf", height = 4, width = 4)

plot(intra_plt)

dev.off()

mistyR::plot_interaction_heatmap(all_misty, "para_ct_5",cutoff = 0)

para_plt <- last_plot() +
  theme(axis.text.x = element_text(size = 12, angle = 90, 
                                   hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size =12)) +
  ggtitle("local neighborhood")

pdf("./results/MI/MOFA_mcell/factor_desc/Factor1_char/misty/para_imp.pdf", height = 4, width = 4)

plot(para_plt)

dev.off()

# Save best in class
all_misty$improvements %>%
  dplyr::filter(measure == "multi.R2",
                value > 30) %>%
  arrange(desc(value)) %>%
  dplyr::mutate(sample = strsplit(sample, "/") %>%
                  map_chr(., ~ .x %>% last())) %>%
  write_csv("./results/MI/MOFA_mcell/factor_desc/Factor1_char/misty/bestinclass.csv")


















