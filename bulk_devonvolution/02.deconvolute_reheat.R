  
## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-19
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
## run mutliple deconvolution methods 
## reference = HCA
## target bulk = chaffin and Reichart HFsc atlases

#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# 
# normalize= args[2]
# normalize_sc= args[3]

set.seed(20)

library(tidyverse)
library(Matrix)
library(Seurat)

directory = "/net/data.isilon/ag-saez/bq_jlanzer/ReHeaT2/"

source(paste0(directory, "R-scripts/utils_decon.R"))
source(paste0(directory, "src/utils1.R"))

#source("analysis/utils_decon.R")

### load data: 

# bulk target::
reheat= readRDS(paste0(directory, "data/metaheart/MetaHeart_RNA_counts.rds"))
#pb= readRDS("data/processed/HF_studiespb.rds")
#pb= readRDS(paste0(directory, "data/HF_studiespb.rds"))

# sc ref: 
seu = readRDS(paste0(directory, "output/HCA/HCA_seu_filtered_nuc.rds"))

print(colnames(seu@meta.data))
print(unique(seu@meta.data$cell_type))
print(wanted_cells)

# make sure to only use samples of interest in reference
#seu = Seurat::subset(seu, cell_type %in% wanted_cells)
seu= subset(seu, subset = (region %in% c("LV", "AX)") & 
                             cell_type %in% wanted_cells &
                             source == "Nuclei" & 
                             Used == "Yes"
                           )
            )

print(unique(seu@meta.data$cell_type))
print("input loaded")

# prepare and normalize -----------------------------------------------------------------------

## we will work with raw counts (counts-> counts)
## and with TPM -> TPM (linear scale)

## TPM single cell:
seu = TPM_normalize_SC(seu, local = "remote")
print("TPMSc")

## TPM bulk
pb2= lapply(reheat, function(x){
  bulk= Normalization(x$gex, method = "TPM", local = "remote")
  })

print("TPMpb")
print(Assays(seu))

##func for single methods deconvolution call:
run_deconvo = function( methods= c("Music", "Bisque", "SCDC"),
                        refSeu,
                        bulkGEX, 
                        assay, 
                        slot, 
                        #bulkProp,
                        ...){
  
  decon_res= list()
  
  for(i in methods){
    # run decon
    print(paste0("running ", i))
    dec_res= sc_deconvo(refSeu= refSeu, 
                        bulkB = bulkGEX,
                        clusters= "cell_type", 
                        samples= "sample", 
                        method = i,
                        slot= slot, 
                        assay= assay,
                        ...)
    
    
    message(paste0(i, " ran"))
    
    #evaluate decon
    
    #dec_res_eval = evaluate_decon_res(dec_res, bulkProp)
    decon_res[[i]]= dec_res
    # decon_res[[i]]= list(predprop= dec_res, 
    #                      eval= dec_res_eval)
  }
  
  return(decon_res)
} 

#calls: 
# 1 TPM
decon_res_tpm= map(pb2,function(x){
  #bulk= Normalization(data= bulks[["T"]], method = x)
  
  decon_res_bulk_normed= run_deconvo(refSeu = seu, 
                                     bulkGEX= as.matrix(x),
                                     assay= "TPM",
                                     slot= "data",
                                     methods = "SCDC"
                                     
  )
  
})

names(decon_res_tpm) = names(pb2)
print("tpm deconvo ran")

# 2 counts
# 
# #create lists of pb counts
# pb3 = list("Chaffin2022"= pb$Chaffin2022$counts,
#            "Reichart2022" = pb$Reichart2022$counts)
# 
# decon_res_counts= map(pb3,function(x){
#   
#   decon_res_bulk_normed= run_deconvo(refSeu = seu, 
#                                      bulkGEX= as.matrix(x),
#                                      assay= "RNA",
#                                      slot= "counts"
#                                      
#   )
#   
# })
# 
# names(decon_res_counts) = names(pb3)
# print("counts deconvo ran")

# 
# saveRDS(list("tpm"= decon_res_tpm, 
#              "counts"= decon_res_counts), 
#         paste0(directory, "output/synced/benchmark_hca_to_hf.rds")
# )
saveRDS( decon_res_tpm, paste0(directory, "output/synced/hca_reheat_scdc.rds"))
