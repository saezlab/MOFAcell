require(liana)

context_df_dict <- readRDS("data/context_df_dict.RDS")

tensor <- liana_tensor_c2c(context_df_dict,
                           score_col = "LRscore",
                           conda_env = "/net/data.isilon/ag-saez/bq_ddimitrov/SOFTWARE/miniconda3/envs/cell2cell/",
                           rank=10)

factors <- format_c2c_factors(tensor)

saveRDS(factors, "data/tensor_factors.RDS")
