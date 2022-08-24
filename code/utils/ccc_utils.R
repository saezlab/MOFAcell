#' Wrapper function to run MOFA /w liana output
#' @param context_df_dict
#'
#' @param score_col name of the column
#'
#' @param data_opts passed to MOFA
#' @param model_opts passed to mofa
#' @param train_opts passed to mofa
#' @inheritDotParams liana_cc2mofa
#'
liana_mofa <- function(context_df_dict,
                       score_col,
                       data_options = NULL,
                       model_options = NULL,
                       train_options = NULL,
                       remove_inactive_factors = FALSE,
                       ...){

    # Prepare ligand-receptor scores in a way that fit mofa
    scores <- liana:::liana_cc2mofa(context_df_dict,
                                    score_col = score_col,
                                    ...)

    # Create object
    liana.mofa <- create_mofa(scores)

    # Prep
    liana.mofa <- prepare_mofa(
        object = liana.mofa,
        data_options = .modify_options(
            get_default_data_options(liana.mofa),
            data_options
        ),
        model_options = .modify_options(
            get_default_model_options(liana.mofa),
            model_options
        ),
        training_options = .modify_options(
            get_default_training_options(liana.mofa),
            train_options
        ))

    # Run
    outfile = file.path(getwd(), "model_simple.hdf5")
    liana.mofa <- run_mofa(liana.mofa,
                           outfile)

    # reload model /w all factors
    liana.mofa <- load_model("model_simple.hdf5",
                             remove_inactive_factors = FALSE)

    liana.mofa@samples_metadata$group <-
        as.factor(gsub("[|].*", # To be fixed
                       "",
                       liana.mofa@samples_metadata$sample))


    return(liana.mofa)
}


#' Helper function to modify MOFA default options
#'
.modify_options <- function(default.opts, passed.opts){
    if(is.null(passed.opts)) return(default.opts)

    # Elements to modify
    to_modify <- intersect(names(default.opts), names(passed.opts))

    # modify default opts
    modifyList(default.opts, passed.opts[to_modify])

}


#' Helper function to filter outlier samples in terms of total counts
#'
#' @param sce SingleCellExperiment
#' @param sample_col Name of the column with sample ids
filter_samples <- function(sce, sample_col, zscore_tresh=-3, plot_dir = NULL){

    pb <- scuttle::aggregateAcrossCells(sce, ids=sce[[sample_col]])

    sample_dist <- scale(colSums(pb@assays@data$counts)) %>%
        as_tibble(rownames="sample",
                  .name_repair = "universal") %>%
        dplyr::rename(zscore = `...1`) %>%
        mutate(keep_total = zscore > zscore_tresh)


    if(!is.null(plot_dir)){
        cairo_pdf(filename = plot_dir,
                  width = 6,
                  height = 5
                  )
        print(sample_dist %>%
            ggplot(aes(x = zscore, fill=keep_total)) +
            geom_histogram(colour="darkgrey") +
            scale_fill_manual(values = c('TRUE' = 'black', 'FALSE' = 'red'),
                              guide = "none") +
            geom_vline(aes(xintercept=zscore_tresh),
                       color="red", linetype="dashed") +
            theme_bw())
        dev.off()
    }


    keep_total <- sample_dist %>% select(sample, keep_total) %>% deframe()
    keep_total <- colData(sce) %>%
        as_tibble(rownames="barcode") %>%
        filter(.data[[sample_col]] %in% names(keep_total[keep_total])) %>%
        pull(barcode)
    # Keep only samples with sum counts within 3 z-scores
    sce <- sce[, keep_total]

    return(sce)
}


#' Function to PCA Regress
#'
#' @param scores scores to be used to separate the PCs from conditions /w ANOVA
#' @param norm whether to normlize the scores or not
#'
pca_regress <- function(scores, norm, rank. = NULL){

    # Check if a row is only 0s
    all_0s <- rowVars(as.matrix(scores))!=0
    if(any(!all_0s)){
        warning("Rows with variance of 0 detected!")
        scores <- scores[all_0s, ]
    }

    # pca
    pca_obj <- prcomp(t(scores), center = norm, scale = norm)

    pca_obj$coords <- pca_obj$x %>%
        data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
        rownames_to_column("patient_region_id") %>%
        as_tibble() %>%
        pivot_longer(-patient_region_id,
                     names_to = "PC", values_to = "value")
    pca_obj$var <- round(summary(pca_obj)$importance[2, ], 2) %>%
        enframe(name = "PC", value = "Expl_var")

    # Regress PCA
    pc_df <- pca_obj$coords %>%
        select_at(c("PC", "value", "patient_region_id")) %>%
        mutate(patient_region_id = gsub("[|].*", "", patient_region_id)) # to cond
    pc_df

    pca_obj$anova_res <- pc_df %>%
        group_by(PC) %>%
        mutate(patient_region_id = as.factor(patient_region_id)) %>%
        nest() %>%
        mutate(anova_res = map(data, function(dat) {
            broom::tidy(aov(value ~ ., data = dat)) %>%
                dplyr::filter(term == "patient_region_id")
        })) %>%
        dplyr::select(PC, anova_res) %>%
        unnest(anova_res) %>%
        left_join(pca_obj$var, by="PC") %>%
        ungroup() %>%
        dplyr::mutate(p_adj = p.adjust(p.value)) %>%
        dplyr::select(PC, p.value, Expl_var, p_adj)

    pca_obj$total_var <- pca_obj$anova_res %>%
        filter(p.value <= 0.05) %>%
        summarise(total_var = sum(Expl_var)) %>%
        pull(total_var)

    return(pca_obj)

}


#' Function to train random forest on PCA coords
run_rf <- function(pca_obj){
    # scores
    rf_data <- pca_obj$coords %>%
        pivot_wider(names_from = PC, values_from = value) %>%
        as.data.frame() %>%
        mutate(patient_group = factor(gsub("[|].*", "", patient_region_id))) %>%
        column_to_rownames("patient_region_id")


    set.seed(1)
    rf <- (ranger(patient_group ~ ., data= rf_data,
                  classification = TRUE,
                  replace = FALSE, seed = 1))


    return(rf)
}

##' Ugly plot_variance explained need to rework to be a heatmap with
##' source and target names on both sides
# plot_variance <- function(mofa.obj){
#     ### Variance explained
#     r2 <- calculate_variance_explained(mofa.obj)
#         r2 <- r2$r2_per_factor$group1
#     # format and order
#     r2_tib <- r2 %>% t() %>%
#         as_tibble(rownames="cellpair")
#
#     # Use these for right and left side of the heatmap
#     # sources <- gsub("&.*", replacement = "", r2_tib$cellpair)
#     # targets <- gsub(".*&", replacement = "", r2_tib$cellpair)
#
#
#     # Sep into heatmap
#     r2_tib %>%
#         # separate(cellpair, into = c("source", "target"), sep = "&") %>%
#         # as.data.frame() %>%
#         column_to_rownames("cellpair") %>%
#         as.matrix() %>%
#         t() %>%
#         ComplexHeatmap::Heatmap(cluster_rows = FALSE,
#                                 cluster_columns = FALSE,
#                                  heatmap_legend_param =
#                                     list(title = "Variance \nExplained (%)",
#                                          fontsize = 20,
#                                          grid_height = unit(10, "mm"),
#                                          grid_width = unit(10, "mm"))
#                                 )
# }
