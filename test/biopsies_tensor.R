require(liana)

context_df_dict <- readRDS("../data/context_df_dict.RDS")

tensor <- liana_tensor_c2c(context_df_dict)

saveRDS(liana::format_c2c_factors(tensor), "data/tensor_factors.RDS")
