context_df_dict <- readRDS("data/context_df_dict.RDS")

tensor <- liana_tensor_c2c(context_df_dict,
                           score_col = "LRscore",
                           conda_env = "cell2cell",
                           how = 'outer',
                           rank=20)

factors <- liana::format_c2c_factors(tensor$factors)

saveRDS(factors, "data/tensor_factors.RDS")

# saveRDS(factors, "/media/dbdimitrov/SSDDimitrov/Repos/biopsies/tensor_factors.RDS")

### Simplest Run
context_df_dict <- readRDS("/media/dbdimitrov/SSDDimitrov/Repos/biopsies/context_df_dict.RDS")[
    str_detect(names(context_df_dict), pattern = "CTRL|IgA|MN")
    ]

tensor <- liana_tensor_c2c(context_df_dict,
                           score_col = "LRscore",
                           conda_env = "cell2cell",
                           how = 'outer',
                           rank=10)
saveRDS(liana::format_c2c_factors(tensor$factors),
        "/media/dbdimitrov/SSDDimitrov/Repos/biopsies/tensor_factors_simple.RDS")
