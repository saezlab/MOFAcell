require(liana)

context_df_dict <- readRDS("data/context_df_dict.RDS")

tensor <- liana_tensor_c2c(context_df_dict,
                           score_col = "LRscore",
                           conda_env = "cell2cell",
                           rank=15)

factors <- format_c2c_factors(tensor)

saveRDS(factors, "data/tensor_factors.RDS")
