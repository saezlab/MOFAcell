context_df_dict <- readRDS("/media/dbdimitrov/SSDDimitrov/Repos/biopsies/context_df_dict.RDS")
# context_df_dict$`FSGS|K68` <- NULL

scores <- liana_cc2mofa(context_df_dict,
                        lr_prop = 0.2,
                        lr_min = 20)

# MOFA
require(MOFA2)

# Create object
liana.mofa <- create_mofa(scores)

# MOFA options
data_opts <- get_default_data_options(liana.mofa)
model_opts <- get_default_model_options(liana.mofa)
model_opts$num_factors <- 15
train_opts <- get_default_training_options(liana.mofa)
train_opts$maxiter <- 250

# Prep
liana.mofa <- prepare_mofa(
    object = liana.mofa,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts)

# Run
outfile = file.path(getwd(), "model.hdf5")
MOFAobject.trained <- run_mofa(liana.mofa,
                               outfile)


MOFAobject.trained@samples_metadata$group <-
    as.factor(gsub("[|].*",
                   "",
                   liana.mofa@samples_metadata$sample))


# Var explained across views by factor
var_explained <- plot_variance_explained(MOFAobject.trained,
                                         x="view",
                                         y="factor") +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust=1))



# Factor loadings
factor_loadings <- plot_factor(MOFAobject.trained,
                               factors = c(1:(MOFA2::get_factors(MOFAobject.trained)[[1]] %>% ncol())),
                               color_by = "group",
                               dot_size = 5, # change dot size
                               dodge = TRUE, # dodge points with different colors
                               legend = FALSE, # remove legend
                               add_violin = TRUE, # add violin plots,
                               violin_alpha = 0.25  # transparency of violin plots
                               )
MOFA2::get_factors(MOFAobject.trained)[[1]] %>%
    as_tibble(rownames="Sample") %>%
    arrange(desc(abs(Factor6)))


overview <- patchwork::wrap_plots(list(factor_loadings,
                                       var_explained
                                       ),
                                  heights = c(2, 1.3),
                                  ncol = 1,
                                  nrow = 2) +
    patchwork::plot_layout(guides = "collect")

grid::grid.draw(overview)





