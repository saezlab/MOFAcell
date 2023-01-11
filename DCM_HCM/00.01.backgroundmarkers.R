# Copyright (c) [2023] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' In this script we calculate markers of cells
#' using edgeR and pseudobulk profiles of all samples

library(SingleCellExperiment)
library(scater)
library(edgeR)
library(tidyverse)

# Importing pb data
pb_data <- read_csv("./data_HCMDCM_Nature/pb_data.csv",
                    show_col_types = FALSE)

colnames(pb_data)[1] <- "sample_id"

pb_data <- pb_data %>%
  column_to_rownames("sample_id") %>%
  as.matrix() %>%
  t()

# Importing coldata of the matrices - ignoring not annotated cell types
coldat <- read_csv("./data_HCMDCM_Nature/pb_coldata.csv",
                   show_col_types = FALSE)[,-1] %>%
  column_to_rownames("colname") %>%
  dplyr::rename(ncells = "counts") %>%
  dplyr::filter(cell_type != "none")

pb_data <- pb_data[,rownames(coldat)]

# Defining cts
cts <- coldat$cell_type %>% 
  unique() %>%
  set_names()

# Pipeline for differential expression

de_res <- map(cts, function(ct) {
  print(ct)
  ct_meta_data <- coldat %>%
    mutate(test_column = ifelse(cell_type == ct, ct, "rest"))
  
  dat <- DGEList(pb_data, samples = DataFrame(ct_meta_data))
  
  keep <- filterByExpr(dat, group = ct_meta_data$test_column)
  
  dat <- dat[keep,]
  
  dat <- calcNormFactors(dat)
  
  design <- model.matrix(~factor(test_column,
                                 levels = c("rest",ct)), dat$samples)
  
  colnames(design) <- c("int", ct)
  
  dat <- estimateDisp(dat, design)
  
  fit <- glmQLFit(dat, design, robust=TRUE)
  
  res <- glmQLFTest(fit, coef=ncol(design))
  
  de_res <- topTags(res, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  return(de_res)
  
})

de_res <- de_res %>% 
  enframe() %>%
  unnest()

de_res %>%
  dplyr::filter(logFC > 0) %>%
  arrange(name, FDR, - logFC) %>%
  write_csv(file = "./results/DCM_broad/edgeR_cellmrkrs.csv")


