# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script I will generalize different tasks
#' that allow to explore an individual Factor within
#' a MOFAcell run

library(MOFA2)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)


get_allfscores <-  function(model, meta, factor, group = FALSE) {
  
  if(group) {
    
    factor_scores <- get_factors(model, factors = "all") %>%
      do.call(rbind, .)
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      left_join(meta, by = "sample")  %>%
      pivot_longer(-colnames(meta), names_to = "Factor")
    
  } else { 
    
    factor_scores <- get_factors(model, factors = "all")[[1]] %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      left_join(meta, by = "sample")  %>%
      pivot_longer(-colnames(meta), names_to = "Factor")
    
    }
  
  return(factor_scores)
  
}

# 0. Get scores and loadings of factor of interest

#' @param model = MOFAcell model
#' @param meta = a data frame with sample + any other colums
#' @param factor = Factor# label
#' 
#' returns a df with the factor scores
get_fscores <- function(model, meta, factor, group = FALSE) {
  
  if(group) {
    
    factor_scores <- (get_factors(model, factors = "all") %>%
      do.call(rbind, .))[,factor, drop= F] %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      left_join(meta, by = "sample")
    
  } else { 
    
    factor_scores <- get_factors(model, factors = "all")[[1]][,factor, drop= F] %>%
      as.data.frame() %>%
      rownames_to_column("sample") %>%
      left_join(meta, by = "sample")
    
  }
  
  return(factor_scores)
  
}

#' @param model = MOFAcell model
#' @param factor = Factor# label
#' 
#' returns a df with the factor loadings

get_floadings <- function(model, factor) {
  
  factor_loadings <- get_weights(model, as.data.frame = T) %>%
    as.data.frame() %>%
    dplyr::mutate(feature = strsplit(as.character(feature), "_") %>%
                    map_chr(., ~ .x[[2]]),
                  ctype = strsplit(as.character(view), "_") %>%
                    map_chr(., ~ .x[[1]])) %>%
    dplyr::rename("factors" = factor) %>%
    dplyr::select(-view) %>%
    dplyr::filter(factors == factor) %>%
    dplyr::select(-factors)
  
  return(factor_loadings)
  
}

# 1. Associate with sample grouping variable

#' @param factor_scores = get_fscores output
#' @param factor = Factor# label
#' @param covar = string with column to be used for pairwise comparisons
#' 
#' returns a df with the factor loadings

get_pw_scorecomp <- function(factor_scores, factor, covar) {
  
  combinations <- combn(factor_scores[,covar] %>%
                          unique,2, simplify = F)
  
  t_comps <- map(combinations, function(comb) {
    
    # Identify proper comparison for unification
    test_term <- levels(comb)
    test_term <- test_term[test_term %in% comb]
    reference_term <- test_term[2]
    test_term <- test_term[1]
    
    comb_scores <- factor_scores %>%
    dplyr::filter(if_any(.cols = all_of(covar),
                         .fns = ~ .x %in% as.character(comb)))
    
    res_tibble <- t.test(as.formula(paste0(factor, " ~ ", covar)),
                      comb_scores) %>%
             broom::tidy() %>%
      dplyr::select(statistic, p.value) %>%
      dplyr::mutate(reference = reference_term,
                    test = test_term)

  }) %>%
    enframe() %>%
    unnest(cols = c(value)) %>%
    dplyr::select(-name) %>%
    dplyr::mutate(p_adj = p.adjust(p.value))
  
  return(t_comps)
  
}

# 2. Correlate gene expression to factors

get_assocgenes <- function(factor_scores,
                           factor,
                           factor_loadings, 
                           pb_data,
                           p_thrs = 0.05, 
                           cor_thrs = 0.5) {
  
  gene_cor <- pb_data %>%
    left_join(factor_scores) %>%
    group_by(view, feature) %>%
    nest() %>%
    mutate(cor_res = map(data, function(dat) {
      cor.test(dat$value, dat[[factor]]) %>%
        broom::tidy()
    })) %>%
    dplyr::select(view, feature, cor_res) %>%
    unnest(c(cor_res)) %>%
    ungroup() %>%
    dplyr::mutate(adj_p = p.adjust(p.value))
  
  sign_genes <- gene_cor %>%
    dplyr::filter(adj_p < p_thrs, 
                  abs(estimate) >= cor_thrs) %>%
    dplyr::mutate(gene = strsplit(feature, "_") %>%
                    map_chr(., ~.x[[2]]))
  
  return(sign_genes)
  
}


get_filtered_loadings <- function(sign_genes,
                                  factor_loadings) {
  genes <- sign_genes %>%
    pull(gene) %>%
    unique()
  
  filt_factor_loadings <- factor_loadings %>% 
    dplyr::filter(feature %in% genes) %>%
    pivot_wider(names_from = ctype, values_from = value,values_fill = 0) %>%
    column_to_rownames("feature") %>%
    as.matrix() 
  
  return(filt_factor_loadings)
  
}

# 3. Create modules -----------------------

get_genemodules <- function(filt_factor_loadings, k = 10) {
  
  gene_clust <- hclust(dist(filt_factor_loadings))
  
  gene_modules <- cutree(gene_clust, k = k)
  
  gene_modules <- gene_modules %>% 
    enframe %>% 
    mutate(value = paste0("module_", value)) %>%
    group_by(value) %>%
    dplyr::rename("feature" = name,
                  "module" = value)
    
  
  return(gene_modules)
}

# 3. Get module cell_type proportion -----------------------

get_genemodules_summ <- function(filt_factor_loadings, gene_modules) {
  
  n_genes <- gene_modules %>%
    group_by(module) %>%
    dplyr::summarise(n_genes = n())
  
  ann_loadings <- filt_factor_loadings %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "ctype") %>%
    dplyr::filter(value != 0) %>%
    left_join(gene_modules, by = "feature") %>%
    left_join(n_genes, by = "module")
  
  loading_ct_comp <- ann_loadings %>%
    group_by(ctype, module) %>%
    mutate(n_cell = n()) %>%
    ungroup() %>%
    dplyr::select(ctype, module, n_genes, n_cell) %>%
    unique() %>%
    dplyr::mutate(prop_cell = n_cell/n_genes) %>%
    dplyr::select(ctype, module, prop_cell) %>%
    pivot_wider(names_from = ctype, 
                values_from = prop_cell,
                values_fill = 0) %>%
    column_to_rownames("module") %>%
    as.matrix() %>%
    t()
  
  mean_loadings <- ann_loadings %>%
    group_by(module) %>%
    summarise(mean_loading = mean(value))
    
    
  return(list("n_genes" = n_genes,
              "loading_ct_comp" = loading_ct_comp,
              "mean_loadings" = mean_loadings))
  
}


summarize_factor <- function(model, 
                             factor, 
                             pb_data, 
                             meta, 
                             covar,
                             out_alias,
                             k = 10,
                             load_p_thrs = 0.05, 
                             load_cor_thrs = 0.5) {
  
  # Create results directory
  res_dir <- paste0(out_alias, "/", factor,"_char/")
  dir_creation <- paste0("mkdir ", res_dir)
  system(dir_creation)
  
  # Get factor scorees
  factor_scores <- get_fscores(model = model,
                               meta = meta,
                               factor = factor)
  # Get loadings
  factor_loadings <- get_floadings(model = model,
                                   factor = factor)
  # Compare the factor scores
  pw_scorecomp <- get_pw_scorecomp(factor_scores = factor_scores,
                   factor = factor,
                   covar = covar)
  
  bplot <- ggplot(factor_scores, aes_string(x = covar, y = factor)) +
    theme_classic() +
    geom_boxplot() +
    geom_point() +
    theme(axis.text = element_text(size =12),
          axis.text.x = element_text(angle = 90))
  
  plot(bplot)
  print(pw_scorecomp)
  
  sign_genes <- get_assocgenes(factor_scores = factor_scores ,
                               factor_loadings = factor_loadings,
                               pb_data = pb_data,
                               factor = factor,
                               p_thrs = load_p_thrs,
                               cor_thrs = load_cor_thrs)
  
  sign_genes_summ <- sign_genes %>% 
    group_by(view) %>% summarise(view_genes = n()) %>%
    arrange(-view_genes) %>%
    dplyr::mutate(view = factor(view))
  
  sign_genes_plt <- sign_genes_summ %>%
    ggplot(aes(x = factor(view,
                          levels = sign_genes_summ$view), 
                          y = view_genes)) +
    geom_bar(stat = "identity", color = "black") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5)) +
    ylab("N. assoc genes") +
    xlab("")
  
  plot(sign_genes_plt)
  
  #associate gene expression with factor scores
  filt_factor_loadings <- get_filtered_loadings(sign_genes = sign_genes,
                                                factor_loadings = factor_loadings)
  
  gene_counter <- (filt_factor_loadings != 0) %>%
    rowSums() %>%
    enframe(name = "gene", value = "n_times")
  
  gene_counter_summ <- gene_counter %>%
    group_by(n_times) %>%
    summarize(n_class = n())
  
  gene_counter_summ_plt <- ggplot(gene_counter_summ, 
         aes(x = n_times, y = n_class)) +
    geom_bar(stat = "identity") +
    theme(axis.text = element_text(size =18),
          axis.title = element_text(size =18))
  
  plot(gene_counter_summ_plt)
  
  
  #modularize
  gene_modules <- get_genemodules(filt_factor_loadings = filt_factor_loadings,
                  k = k)
  
  
  genemodules_summ <- get_genemodules_summ(filt_factor_loadings = filt_factor_loadings,
                                           gene_modules = gene_modules)
  
  
  # Print complex
  col_fun <- circlize::colorRamp2(c(-1,0, 1), c("blue", "white", "red"))
  col_fun_main <- circlize::colorRamp2(c(0, 1), c("white", "red"))
  
  n_genes_module <- set_names(genemodules_summ$n_genes$n_genes,
                              genemodules_summ$n_genes$module)
  
  n_genes_module <- n_genes_module[colnames(genemodules_summ$loading_ct_comp)]
  
  module_loading <- set_names(genemodules_summ$mean_loadings$mean_loading,
                              genemodules_summ$mean_loadings$module)
  
  module_loading <- module_loading[colnames(genemodules_summ$loading_ct_comp)]
  
  ha <- HeatmapAnnotation("size" = anno_barplot(n_genes_module), 
                          "mean_load" = module_loading,
                          col = list("size" = "black",
                                     "mean_load" = col_fun))
  
  hmap <- Heatmap(genemodules_summ$loading_ct_comp, 
          top_annotation = ha, 
          show_column_dend = F, 
          col = col_fun_main,name = "cell_prop")
  
  draw(hmap)
  
  # Just plot the 
  
  #Make a complex heatmap of loadings grouped by modules
  
  filt_factor_loadings <- filt_factor_loadings[gene_modules$feature, ] 
  
  loads_hmap <- Heatmap(filt_factor_loadings, 
          right_annotation = rowAnnotation(module = gene_modules$module),
          row_split = gsub("module_", "",gene_modules$module),
          cluster_rows = TRUE,
          cluster_columns = TRUE, 
          show_column_names = TRUE,
          show_row_names = FALSE,
          show_column_dend = FALSE,
          show_row_dend = FALSE)
  
  draw(loads_hmap)
  
  loads_hmap_simple <- Heatmap(filt_factor_loadings,
                               cluster_rows = TRUE,
                               cluster_columns = TRUE, 
                               show_column_names = TRUE,
                               show_row_names = FALSE,
                               show_column_dend = FALSE,
                               show_row_dend = FALSE,
                               name = "loading")
  
  draw(loads_hmap_simple)
  
  # Save stuff
  
  write_csv(genemodules_summ$mean_loadings,
            paste0(res_dir, "mean_module_loadings.csv"))
  
  write_csv(genemodules_summ$mean_loadings,
            paste0(res_dir, "mean_module_loadings.csv"))
  
  write_csv(genemodules_summ$loading_ct_comp %>%
              as.data.frame(),
            paste0(res_dir, "module_comp.csv"))
  
  write_csv(gene_modules,
            paste0(res_dir, "gene_modules.csv"))
  
  write_csv(pw_scorecomp,
            paste0(res_dir, "pw_comp.csv"))
  
  filt_factor_loadings %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, 
                 names_to = "celltype",
                 values_to = "value") %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(dir = ifelse(value > 0 , "pos", "neg")) %>%
    dplyr::mutate(celltype = paste0(celltype, "_", dir)) %>%
    dplyr::select(-dir) %>%
    dplyr::mutate(value = abs(value)) %>%
    write_csv(paste0(res_dir, "loadings.csv")) 
  
  pdf(paste0(res_dir, "pw_comp.pdf"), height = 6, width = 4)
  
  plot(bplot)
  
  dev.off()
  
  pdf(paste0(res_dir, "module_summ.pdf"), height = 5, width = 8)
  
  draw(hmap)
  
  dev.off()
  
  pdf(paste0(res_dir, "n_genes.pdf"), height = 3, width = 3)
  
  plot(sign_genes_plt)
  
  dev.off()
  
  pdf(paste0(res_dir, "loadings_module.pdf"), height = 5, width = 8)
  
  draw(loads_hmap)
  
  dev.off()
  
  pdf(paste0(res_dir, "loadings.pdf"), height = 5, width = 8)
  
  draw(loads_hmap_simple)
  
  dev.off()
  
  
  
}


# Check that by default loadings are put to 0
summarize_factor_lite <- function(model, 
                             factor, 
                             meta, 
                             out_alias,
                             load_cor_thrs = 0.5) {
  
  # Create results directory
  res_dir <- paste0(out_alias, "/", factor,"_char/")
  dir_creation <- paste0("mkdir ", res_dir)
  system(dir_creation)
  
  # Get factor scorees
  factor_scores <- get_fscores(model = model,
                               meta = meta,
                               factor = factor)
  
  # Get loadings
  factor_loadings <- get_floadings(model = model,
                                   factor = factor) %>%
    dplyr::filter(value != 0)
  
  write_csv(factor_loadings, paste0(res_dir, "loadings_all.csv"))
  
  
  loading_hist <- ggplot(factor_loadings, aes(x = value, color = ctype, fill = ctype)) +
    geom_histogram(bins = 20, alpha =0.5) + theme_classic() +
    theme(axis.text = element_text(size = 9)) +
    xlab("loading") +
    facet_wrap(.~ctype)
  
  pdf(paste0(res_dir, "loadings_hist.pdf"), height = 4, width = 6)
  
  plot(loading_hist)
  
  dev.off()
  
  # Get number of tested genes
  
  print("Unique genes tested")
  
  print(factor_loadings$feature %>%
          unique() %>%
          length())
  
  factor_loadings_filt <- factor_loadings %>%
    dplyr::filter(abs(value) >= load_cor_thrs)
    
  # Filter genes based on a loading cut-off
  
  model_size <- factor_loadings_filt %>%
    dplyr::mutate(direction = ifelse(value > 0 , "Ischemic-like", "Myogenic-like")) %>%
    group_by(ctype, direction) %>%
    summarise(used_genes = length(feature)) %>%
    arrange(-used_genes)
  
  model_size_bar <- model_size %>%
    ggplot(aes(x = ctype, y = log10(used_genes), fill =direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_classic() +
    theme(axis.text = element_text(size = 9),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = c("#FF6666", "#4169E1")) +
    xlab("") +
    ylab("log10(Number of genes)")
  
  pdf(paste0(res_dir, "model_size_bar.pdf"), height = 2.3, width = 3.2)
  
  plot(model_size_bar)
  
  dev.off()
  
  model_size_bar
  
  write_csv(model_size, paste0(res_dir, "model_size.csv"))
  
  # Calculate the expression of genes across cells
  
  shared_genes <- factor_loadings_filt %>%
    group_by(feature) %>%
    summarize(n_cells = length(feature)) %>%
    arrange(-n_cells) %>%
    group_by(n_cells) %>%
    summarize(count = length(n_cells))
  
  shared_genes_plt <- ggplot(shared_genes, aes(x = log10(count), y = paste(n_cells, "types"))) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12)) +
    theme(plot.margin = unit(c(0,1,0.2,0.2), "cm")) +
    ylab("number of cell types") +
    xlab("log10(number of genes)")
  
  pdf(paste0(res_dir, "shared_genes_plt.pdf"), height = 2.3, width = 2.7)
  
  plot(shared_genes_plt)
  
  dev.off()
    
  shared_genes %>%
    dplyr::mutate(multicellular = ifelse(n_cells >1, "Mcell", "Ucell")) %>%
    group_by(multicellular) %>%
    summarise(count = sum(count)) %>%
    dplyr::mutate(percent = count/sum(count)) %>%
    print()
  
  # Here we calculate the jaccard index between shared genes across cell - types
  
  disease_genes <- factor_loadings_filt %>%
    dplyr::mutate(direction = ifelse(value > 0 , "Disease", "Healthy")) %>%
    dplyr::filter(., 
                  direction == "Disease") %>%
    dplyr::select(-direction) %>%
    group_by(ctype) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x[[1]])) %>%
    deframe()
  
  healthy_genes <- factor_loadings_filt %>%
    dplyr::mutate(direction = ifelse(value > 0 , "Disease", "Healthy")) %>%
    dplyr::filter(., 
                  direction == "Healthy") %>%
    dplyr::select(-direction) %>%
    group_by(ctype) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x[[1]])) %>%
    deframe()
  
  disease_jaccard <- sapply(disease_genes, function(x){
    sapply(disease_genes, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  
  disease_jaccard[upper.tri(disease_jaccard, diag = T)] = NA
  disease_jaccard <- disease_jaccard %>%
    as.data.frame() %>%
    rownames_to_column("cell_a") %>%
    pivot_longer(-cell_a, names_to = "cell_b", values_to = "jaccard_ix") %>%
    na.omit() %>%
    dplyr::mutate(direction = "Disease")
  
  healthy_jaccard <- sapply(healthy_genes, function(x){
    sapply(healthy_genes, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  
  healthy_jaccard[lower.tri(healthy_jaccard, diag = T)] = NA
  healthy_jaccard <- healthy_jaccard %>%
    as.data.frame() %>%
    rownames_to_column("cell_a") %>%
    pivot_longer(-cell_a, names_to = "cell_b", values_to = "jaccard_ix") %>%
    na.omit() %>%
    dplyr::mutate(direction = "Healthy")
  
  all_jaccard <- bind_rows(disease_jaccard, healthy_jaccard)
  
  all_jaccard_plt <- all_jaccard %>%
    ggplot(aes(x = cell_b, y = cell_a, alpha = jaccard_ix, fill = direction)) +
    geom_tile() +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.border = element_rect(colour = "black", fill=NA)) +
    scale_fill_manual(values = c("#FF6666", "#4169E1")) +
    coord_equal() +
    xlab("") + ylab("")
  
  plot(all_jaccard_plt)
  
    
  pdf(paste0(res_dir, "jaccard_ix.pdf"), height = 3, width = 3.5)
  
  plot(all_jaccard_plt)
  
  dev.off()
  
  write_csv(all_jaccard, paste0(res_dir, "jaccard_ix.csv"))
  
  # Select genes to be plotted
  
  sel_genes <- factor_loadings_filt %>%
    pull(feature) %>%
    unique()
  
  feature_matrix <- factor_loadings %>%
    dplyr::filter(feature %in% sel_genes) %>%
    pivot_wider(names_from = ctype, values_from = value,values_fill = 0) %>%
    column_to_rownames("feature") %>%
    as.matrix()
  
  #Make a complex heatmap of loadings grouped by modules
  
  col_fun_fact <- circlize::colorRamp2(seq(-1.5, 1.5, length = 30), hcl.colors(30,"RdBu",rev = T))
  
  loads_hmap <- Heatmap(feature_matrix,
                        name = "loading",
                        cluster_rows = TRUE,
                        cluster_columns = TRUE, 
                        show_column_names = TRUE,
                        show_row_names = FALSE,
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        col = col_fun_fact)
  
  draw(loads_hmap)
  
  # Save stuff
  
  factor_loadings_filt %>%
    dplyr::rename("celltype" = "ctype",
                  "gene" = "feature") %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(dir = ifelse(value > 0 , "pos", "neg")) %>%
    dplyr::mutate(celltype = paste0(celltype, "_", dir)) %>%
    dplyr::select(-dir) %>%
    dplyr::mutate(value = abs(value)) %>%
    write_csv(paste0(res_dir, "loadings.csv")) 
 
  pdf(paste0(res_dir, "loadings.pdf"), height = 3.5, width = 2.5)
  
  draw(loads_hmap)
  
  dev.off()
  
}










