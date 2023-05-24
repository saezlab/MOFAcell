## utils
col.set <- c("#fb8500", "#0099cc", "#f60000", "#fde800","#00b600","#00896e", "#0053c8", "#ff3399", "#b374e0", "#00cc99",  "#fbc9e2","#0000cc")


wanted_cells= c("Ventricular_Cardiomyocyte", 
                "Endothelial", 
                #"Adipocytes", 
                "Lymphoid", 
                "Fibroblast", 
                "Pericytes", 
                #"Neuronal", 
                "Smooth_muscle_cells", 
                "Myeloid")



# pseudobulk ----------------------------------------------------------------------------------
# testingS: 

# randomize cells from the GEX matrix

pseudomize= function(seu_obj,
                     counts= T,
                     permutations= 10,
                     poolsize = 1000, 
                     col_celltype= "CellType", 
                     seed= 20){
  
  set.seed(seed)
  
  meta= seu_obj[[]]
  
  all_cells= unique(meta[,col_celltype])
  
  DefaultAssay(seu_obj)= "RNA"
  
  if(counts!=T){
    message(" only works with count data so far. ")
  }else{
    GEX = GetAssayData(seu_obj, slot ="counts") #work with counts as a default
  }
  set.seed(seed)
  ## create pseudo bulk samples 
  pseudobulknames= paste0("psb_", c(1:permutations))
  
  #list of a pseudobulk and cell proportion
  pseudo_list = map(pseudobulknames, function(x){
    randomcells= sample(colnames(GEX), poolsize)
    randomcelltypes= meta[randomcells,col_celltype]
    
    
    cell_props = enframe(prop.table(table(randomcelltypes)), name = "celltype", value= x)
    
    # append for cell types that turned out to be 0 
    
    #missing.CT = cell_props$celltype[!cell_props$celltype %in% all_cells]
    missing.CT = all_cells[!all_cells %in% cell_props$celltype]
    
    if(length(missing.CT)>0){
      missing.CT = data.frame(celltype = missing.CT, val = rep(0, length(missing.CT)), stringsAsFactors = FALSE)
      colnames(missing.CT)[colnames(missing.CT)=="val"]= x
      cell_props = rbind(cell_props, missing.CT)
    }
    
    bulk_matrix  = GEX[, randomcells]
    
    bulk_vec= enframe(rowSums(as.matrix(bulk_matrix)), name= "gene", value = x)# %>% column_to_rownames("gene")
    
    list("props"= cell_props,"gex"= bulk_vec)
  })
  names(pseudo_list)= pseudobulknames
  
  
  #change format to create 
  # A) matrix of gene expression (sum of counts) where each column is a pseudobulk sample
  # B) matrix of cell proportions where each column is also a pseudobulk sample
  
  gexes= lapply(pseudo_list, function(x){
    x$gex %>% column_to_rownames("gene")
  })
  A= do.call(cbind, gexes)
  
  props= lapply(pseudo_list, function(x){
    x$props %>% column_to_rownames("celltype")
  }) 
  B= do.call(cbind, props)
  
  # return both dfs:
  return(list(GEX_pseudobulk= A, cell_proportions= B))
  
  
}


#randomize cell type ditribution and than randomize cells to fill it up
Generator <- function(sce, 
                      phenoData, 
                      Num.mixtures = 1000, 
                      pool.size = 100,
                      min.percentage = 1,
                      max.percentage = 99, 
                      seed = 24){ 
  
  CT = unique(phenoData$cellType)
  ?stopifnot(length(CT) >= 2)
  
  set.seed(seed)
  require(dplyr)
  require(gtools)
  
  cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE) 
  colnames(cell.distribution) = c("CT","max.n")
  
  Tissues = list()
  Proportions = list()
  
  for(y in 1:Num.mixtures){
    
    #Only allow feasible mixtures based on cell distribution
    while(!exists("P")){
      
      num.CT.mixture = sample(x = 2:length(CT),1)
      selected.CT = sample(CT, num.CT.mixture, replace = FALSE)
      
      P = runif(num.CT.mixture, min.percentage, max.percentage) 
      P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
      P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)
      
      missing.CT = CT[!CT %in% selected.CT]
      missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)
      
      P = rbind.data.frame(P, missing.CT)
      potential.mix = merge(P, cell.distribution)
      potential.mix$size = potential.mix$expected * pool.size
      
      if( !all(potential.mix$max.n >= potential.mix$size) | sum(P$expected) != 1){
        rm(list="P") 
      }
      
    }
    
    # Using info in P to build T simultaneously
    chosen_cells <- sapply(which(potential.mix$expected != 0), function(x){
      
      n.cells = potential.mix$expected[x] * pool.size
      chosen = sample(phenoData$cellID[phenoData$cellType == potential.mix$CT[x]],
                      n.cells)
      
      chosen
    }) %>% unlist()
    
    
    T <- Matrix::rowSums(sce[,colnames(sce) %in% chosen_cells]) %>% as.data.frame()
    colnames(T) = paste("mix",y,sep="")
    
    P = P[,c("CT","expected")]
    P$mix = paste("mix",y,sep="")
    
    Tissues[[y]] <- T
    Proportions[[y]] <- P
    
    rm(list=c("T","P","chosen_cells","missing.CT"))
    
  }
  
  P = do.call(rbind.data.frame, Proportions)
  T = do.call(cbind.data.frame, Tissues)
  
  P = data.table::dcast(P, CT ~ mix, 
                        value.var = "expected",
                        fun.aggregate = sum) %>% data.frame(.,row.names = 1) 
  
  P = P[,gtools::mixedsort(colnames(P))]
  
  return(list(T = T, P = P))
  
} 

# deconvo functions ---------------------------------------------------------------------------

## bayes prism 
#' @param refSeu, reference Seurat object. most contian raw counts and is prefiltered 
#' @param refA, reference count matrix, needs to be implemented
#' @param bulkB, bulk matrix to deconvolute, cols= samples, rows= genes
#' @param cell_types, pheno types for refA, needs to be implemented as well

runBayesPrism = function(refA= NULL, 
                         refSeu= NULL, 
                         bulkB, 
                         cell_annos= NULL){
 
   if(is.null(refA) & is.null(refSeu)){
    message("add reference matrix or reference seurat object")
  }
  
  require(TED)
  
  if(is.null(refA)){
    refA= refA= GetAssayData(refSeu, slot = "counts", assay = "RNA")
    
  }
  
  # filter genes in reference for bulk 
  genes= rownames(bulkB)
  refA = refA[rownames(refA) %in% genes, ]
  refA = t(as.matrix(refA))
  
  # prepare phenotype input
  if(is.null(cell_annos)){
    Idents(refSeu) = refSeu[[]]$cell_type
    pheno.lab = as.data.frame(Idents(refSeu))
    cell_annos= pheno.lab  %>% pull(`Idents(refSeu)`)
    
  }
  
  # prepare bulk input 
  x= t(bulkB)
  
  print("input prepared")
  
  ### run bayes prism: 
  TED.res= run.Ted(ref.dat = refA,
                   cell.type.labels = cell_annos,
                   X= x,
                   input.type = "scRNA")
  
  
  return(TED.res)
  }



## MuSiC
#' @param refSeu, reference Seurat object. most contian raw counts and is prefiltered 
#' @param refA, reference count matrix
#' @param bulkB, bulk matrix to deconvolute, cols= samples, rows= genes
#' @param pheno_data, meta data for sc referecne
#' @param clusters, column name in pheno data that contain the cell types to be clustered on,
#' @param samples, column name in pheno_data that label the individuals in the sc-reference

runMuSiC = function(refA= NULL, 
                    refSeu= NULL, 
                    bulkB, 
                    pheno_data= NULL,
                    clusters = 'cell_type',
                    samples = 'sample',
                    slot= "data",
                    assay= "SCT",
                    ...){
  require(xbioc)
  require(Biobase)
  
  if(is.null(refA) & is.null(refSeu)){
    message("add reference matrix or reference seurat object")
  }
  
  if(is.null(refA)){
    refA= GetAssayData(refSeu, slot = "counts", assay = "RNA")
  }
  
  ### prepare bulk
  bulk.eset <- ExpressionSet(assayData=as.matrix(bulkB))
  
  ### prepare sc reference
  if(is.null(pheno_data)){
    pheno_data= refSeu[[]]
  }
  
  if(all(rownames(pheno_data)== colnames(refA))){
    message("columns for sc eset agree: ")
  }
  
  phenoData <- new("AnnotatedDataFrame",
                   data=pheno_data)
  
  # filter genes in reference for bulk 
  genes= rownames(bulkB)
  refA = refA[rownames(refA) %in% genes, ]
  
  #create sc eset
  sc.eset <- ExpressionSet(assayData = as.matrix(refA), 
                           phenoData = phenoData)
  
  
  #run music
  
  Est.prop = music_prop(bulk.eset = bulk.eset, 
                        sc.eset = sc.eset, 
                        clusters = clusters,
                        samples = samples, 
                        select.ct =NULL , 
                        verbose = F)
                        #...)
  
  return(Est.prop)
}




## SCDC
#' @param refSeu, reference Seurat object. most contian raw counts and is prefiltered 
#' @param refA, reference count matrix
#' @param bulkB, bulk matrix to deconvolute, cols= samples, rows= genes
#' @param pheno_data, meta data for sc referecne
#' @param clusters, column name in pheno data that contain the cell types to be clustered on,
#' @param samples, column name in pheno_data that label the individuals in the sc-reference

runSCDC = function(refA= NULL, 
                    refSeu= NULL, 
                    bulkB, 
                    pheno_data= NULL,
                    clusters = 'cell_type',
                    samples = 'sample',
                    slot="data",
                    assay="SCT"){
  
  if(is.null(refA) & is.null(refSeu)){
    message("add reference matrix or reference seurat object")
  }
  
  require(SCDC)
  require(xbioc)
  require(Biobase)
  
  if(is.null(refA)){
    refA= GetAssayData(refSeu, slot = slot, assay = assay)
  }
  
  #prepare bulk
  bulk.eset <- ExpressionSet(assayData=as.matrix(bulkB))
  
  
  #prepare sc reference
  if(is.null(pheno_data)){
    pheno_data= refSeu[[]]
  }
  print("columns for sc eset agree: ")
  all(rownames(pheno_data)== colnames(refA))
  
  phenoData <- new("AnnotatedDataFrame",
                   data=pheno_data)
  
  
  # filter genes in reference for bulk 
  # genes= rownames(bulkB)
  # refA = refA[rownames(refA) %in% genes, ]
  # refA = t(as.matrix(refA))
  # 
  sc.eset <- ExpressionSet(assayData = as.matrix(refA), 
                           phenoData = phenoData)
  
  cells= unique(pheno_data[[clusters]])
  #run music
  
  Est.prop = SCDC_prop(bulk.eset = bulk.eset, 
                        sc.eset = sc.eset, 
                        ct.varname = clusters,
                        sample = samples, 
                       ct.sub = cells)
  
  return(Est.prop)
}

## CIBERSORT
#' @param refSeu, reference Seurat object. most contian raw counts and is prefiltered 
#' @param refA, reference count matrix
#' @param bulkB, bulk matrix to deconvolute, cols= samples, rows= genes
#' @param pheno_data, meta data for sc referecne
#' @param clusters, column name in pheno data that contain the cell types to be clustered on,
#' @param samples, column name in pheno_data that label the individuals in the sc-reference

runCibersort = function(refA= NULL, 
                   #refSeu= NULL, 
                   bulkB, 
                   #pheno_data= NULL,
                   #clusters = 'cell_type',
                   #samples = 'sample',
                   #slot="data",
                   #assay="SCT",
                   local = "remote"){
  
  if(local =="local"){
    source("src/CIBERSORT.R")
  }else if(local =="remote"){
    source("/net/data.isilon/ag-saez/bq_jlanzer/ReHeaT2/R-scripts/src/CIBERSORT.R")
  }
  if(is.null(refA) & is.null(refSeu)){
    message("add reference matrix or reference seurat object")
  }
  
  if(is.null(refA)){
    refA= GetAssayData(refSeu, slot = slot, assay = assay)
  }
  
  #prepare bulk
  bulk.eset <- ExpressionSet(assayData=as.matrix(bulkB))
  
  
  #prepare sc reference
  if(is.null(pheno_data)){
    pheno_data= refSeu[[]]
  }
  print("columns for sc eset agree: ")
  all(rownames(pheno_data)== colnames(refA))
  
  phenoData <- new("AnnotatedDataFrame",
                   data=pheno_data)
  
  
  # filter genes in reference for bulk 
  # genes= rownames(bulkB)
  # refA = refA[rownames(refA) %in% genes, ]
  # refA = t(as.matrix(refA))
  # 
  sc.eset <- ExpressionSet(assayData = as.matrix(refA), 
                           phenoData = phenoData)
  
  cells= unique(pheno_data[[clusters]])
  #run music
  
  Est.prop = SCDC_prop(bulk.eset = bulk.eset, 
                       sc.eset = sc.eset, 
                       ct.varname = clusters,
                       sample = samples, 
                       ct.sub = cells)
  
  return(Est.prop)
}

###### deconvo function wrapper
#' @param refSeu, reference Seurat object. most contian raw counts and is prefiltered 
#' @param refA, reference count matrix
#' @param bulkB, bulk matrix to deconvolute, cols= samples, rows= genes
#' @param pheno_data, meta data for sc referecne (not necessary if refSeu is provided)
#' @param clusters, column name in pheno data that contain the cell types to be deconvoluted on,
#' @param samples, column name in pheno_data that label the individuals in the sc-reference

sc_deconvo= function(refA= NULL,    #sc
                     refSeu= NULL,  #sc
                     bulkB,   #bulk
                     pheno_data= NULL, #sc
                     clusters = 'cell_type', #sc
                     samples = 'sample', #sc
                     method = c("Bisque", "Music", "SCDC"), 
                     slot,
                     assay
                     ){
  require(xbioc)
  require(Biobase)
  
  if(is.null(refA) & is.null(refSeu)){
    message("add reference matrix or reference seurat object")
  }
  print("start")
  if(is.null(refA)){
    refA= GetAssayData(refSeu, slot = slot, assay = assay)
  }
  print("refA done")
  #print(class(bulkB))
  #print(bulkB[1:10, 1:10])
  ### prepare bulk
  bulk.eset <- ExpressionSet(assayData=as.matrix(bulkB))
  
  #print("bulk eset")
  
  ### prepare sc reference
  if(is.null(pheno_data)){
    pheno_data= refSeu[[]]
  }
  
  
  if(all(rownames(pheno_data)== colnames(refA))){
    message("columns for sc eset agree ")
  }
  
  phenoData <- new("AnnotatedDataFrame",
                   data=pheno_data)
  #print("phenodata")
  # filter genes in reference for bulk 
  genes= rownames(bulkB)
  refA = refA[rownames(refA) %in% genes, ]
  #print("refA subset")
  #create sc eset
  sc.eset <- ExpressionSet(assayData = as.matrix(refA), 
                           phenoData = phenoData)
  print("sc.eset")
  
  if (method == "Music"){
    
    require(MuSiC)
    print("go Music")
    Est.prop = t(MuSiC::music_prop(bulk.eset = bulk.eset, 
                          sc.eset = sc.eset, 
                          clusters = clusters,
                          samples = samples, 
                          select.ct =NULL , 
                          verbose = F)$Est.prop.weighted) #add this for pipeline
  
  }else if (method == "Bisque"){
    
    require(BisqueRNA)
    
    Est.prop <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = bulk.eset,
                                                       sc.eset = sc.eset,
                                                       markers = NULL,
                                                       cell.types = clusters,
                                                       subject.names = samples,
                                                       use.overlap = FALSE,
                                                       verbose = FALSE,
                                                       old.cpm = FALSE)$bulk.props
  }else if (method == "SCDC"){
    
    require(SCDC)
    
    cells= unique(as.character(pheno_data[[clusters]]))
    
    Est.prop <- t(SCDC::SCDC_prop(bulk.eset = bulk.eset,
                                  sc.eset = sc.eset,
                                  ct.varname = clusters,
                                  sample = samples,
                                  ct.sub = cells, 
                                  iter.max = 200)$prop.est.mvw)
  }
  
  
  return(Est.prop)
}

##func:
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
                        assay= assay,
                        slot= slot,
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


# deconvo results -----------------------------------------------------------------------------
####### function to evaluate deconvolution results
#' @param pred_prop, predicted proportions (decon results)
#' @param real_prop, ground truth of proportions
#' 

evaluate_decon_res= function(pred_prop, real_prop, method= "decon"){
  require(pheatmap)
  require(tidyverse)
  require(ComplexHeatmap)
  require(cowplot)
  require(reshape2)
  require(reshape)
  
  pred_prop= as.matrix(pred_prop)
  real_prop= as.matrix(real_prop)
  
  #ensure same ordering: 
  pred_prop = pred_prop[sort(rownames(pred_prop)), sort(colnames(pred_prop))]
  real_prop = real_prop[sort(rownames(real_prop)), sort(colnames(real_prop))]
  
  rownames(real_prop)
  class(pred_prop)
  class(real_prop)
  
  p.heat.pred= pred_prop %>% 
    reshape2::melt(., value.name ="proportion", varnames= c("celltype", "sample")) %>% 
    ggplot(., aes(x= sample, y= celltype, fill = proportion))+
    geom_tile()+
    scale_fill_gradient(low="grey", high="darkred", limits = c(0,1))+
    labs(x = "pseudobulk samples", 
         y= "")+
    ggtitle("predicted")+
    theme(axis.text.x = element_blank())
  
  p.heat.real= real_prop %>% 
    reshape2::melt(., value.name ="proportion", varnames= c("celltype", "sample")) %>%
    ggplot(., aes(x= sample, y= celltype, fill = proportion))+
    geom_tile()+
    scale_fill_gradient(low="grey", high="darkred",limits = c(0,1))+
    labs(x = "pseudobulk samples", 
         y= "")+
    ggtitle("truth")+
    theme(axis.text.x = element_blank())
  
  #hm1= Heatmap(pred_prop, name= "pred")
  #hm2= Heatmap(real_prop, name= "truth")
  
  # bring in same order
  
  #only run if cell types agree: 
  if(all(colnames(pred_prop) == colnames(real_prop)) &
     all(rownames(pred_prop)== rownames(real_prop))){
    
    print("dfs agree")
    
    ##1 calulcate RMSE
    delta = as.matrix(pred_prop) - as.matrix(real_prop)
    #delta[1:9, 1:10]
    cell_rmse= sqrt(rowMeans(delta^2)) 
    total_rmse= sqrt(mean(delta^2))
    
    patient_rmse= sqrt(colMeans(delta^2)) 
    
    total_rmse
    
    
    
    p.rmse = as_data_frame(cell_rmse) %>% mutate(cell_type = names(cell_rmse),
                                            cells= "cells") %>%
      ggplot(., aes(x= cells, y= value))+
      geom_boxplot()+
      geom_jitter(size= 4, alpha = 0.8, aes(color = cell_type))+
      #geom_point(size= 4, aes(color = cell_type), alpha = 0.5)+
      scale_color_manual(values= col.set)+
      labs(y= "rmse", 
           x= "")
    
    ##2 calculate pearson
    #correlating per cell type (thats why transpose)
    
    #rownames(pred_prop) = paste0(rownames(pred_prop),"pred")
    #pred_prop[pred_prop == 0] = 0.000001
    cell_cors= diag(cor(t(pred_prop), t(real_prop), method= "pearson"))
    cor(t(pred_prop), t(real_prop), method= "pearson")
    
    ##3 
    # correlating per cell type
    pat_cors= diag(cor((pred_prop), (real_prop), method= "pearson"))
    
    tidyres= melt(pred_prop, varnames = c("celltype", "sample")) %>%
      dplyr::rename(pred= value) %>%
      left_join( melt(as.matrix(real_prop), varnames = c("celltype", "sample")))%>%
      dplyr::rename(truth= value)
    
    #plot direct: 
    p_prop = tidyres %>%
      ggplot(., aes(x= truth, y= pred))+
      geom_point()+
      geom_abline(size= 1, color = "darkgrey")+
      geom_smooth(method='lm', formula= y~x, fullrange = T)+
      ylim(c(0,0.5))+
      xlim(c(0,0.5))+
      facet_grid(cols= vars(celltype ))+
      theme_minimal()+
      theme(axis.text.x = element_text(angle= 60, hjust= 1))
      
     results=tidyres %>% dplyr::summarise(RMSE = sqrt(mean((pred-truth)^2)) %>% round(.,4),
                                          Pearson=cor(pred,truth) %>% round(.,4))## mean 
    
      p.cor= as_data_frame(cell_cors) %>% mutate(cell_type = names(cell_cors), 
                                          cells= "cells") %>% 
      dplyr::rename(correlation= value) %>%
      ggplot(., aes(x= cells, y= correlation))+
      geom_boxplot()+
      #geom_point(size= 4, alpha = 0.8, aes(color = cell_type))+
      geom_jitter(size= 4, alpha = 0.8, aes(color = cell_type))+
      scale_color_manual(values= col.set)+
      labs(x= "",y= "pearson.cor")
    
    p2= plot_grid(p.rmse+ theme(legend.position = "none"), p.cor, rel_widths = c(1,1.5))
    p1= plot_grid(p.heat.real, p.heat.pred)
    
    #p.sum= plot_grid(p1, p2, ncol = 1)
    p.sum2= plot_grid(p2, p_prop, ncol = 1)
    #p.sum
    
    return(list(plot = p.sum2, 
                rmse= list(rmse_cell = cell_rmse, rmse_total=  results$RMSE), 
                corr = list(cors_cell= cell_cors, cors_total=  results$Pearson),
                pat_cors= pat_cors,
                pat_rmse= patient_rmse,
                props= list("real"= real_prop, 
                            "pred"=pred_prop)
                )
           ) 
    
  }else{ message("reorder") }
  
}





# PREPROCESSING -------------------------------------------------------------------------------

## classic bulk preprocesing
Normalization <- function(data, 
                          method,
                          sc_pheno =NULL,
                          celltype= "cell_type",
                          local = "remote"){
  
  require(edgeR)
  
  if(is.null(sc_pheno)){
    groups = colnames(data)
    }else{
    groups= as.character(sc_pheno[[celltype]][rownames(sc_pheno) %in% colnames(data)])
  }
  
  dge <- edgeR::DGEList(data, group = groups)
  message(paste0("dge object created", method))
  #usuallly genes are filtered, for minimal GEX 
  #dge <- dge[keep,,keep.lib.sizes=FALSE] 
  dge <- calcNormFactors(dge, method = "TMM")
  message(paste0("dge Factors calculated", method))
  #all(data==dge$counts)
  
  # use limma voom to transform count data to log2 counts per million
  if(method == "voom"){
    normed <- voom(dge, plot=FALSE)$E
  }else if(method == "logCPM"){
    normed= edgeR::cpm(dge, log= T, prior.count = 2)
  }else if(method == "CPM"){
    normed= edgeR::cpm(dge, log= F, normalized.lib.sizes= T)
  } else if (method == "none"){
    normed= dge$counts
  } else if (method == "SCtransform"){ #not working yet
    matrix = sctransform::vst(matrix, return_corrected_umi=TRUE, show_progress = FALSE)$umi_corrected
  } else if (method == "TPM"){
    
    #load gene lengths= 
    if (local == "remote"){
      gene_info= readRDS("/net/data.isilon/ag-saez/bq_jlanzer/ReHeaT2/output/Gene_lengths.rds")
    } else if (local =="local"){
      gene_info = readRDS("~/R-projects/ReHeaT2/output/synced/Gene_lengths.rds")
    }
    #  filter genes in data and bring in order
    filter_data = data[rownames(data) %in% gene_info$hgnc_symbol, ]
    g.lenghts= gene_info$length[match(rownames(filter_data), table = gene_info$hgnc_symbol)]
    # TPM func
    tpm3 <- function(counts,len) {
      x <- counts/len
      tpm = t(t(x)*1e6/colSums(x))
      # normed= x*1e6/colSums(x)
      return(tpm)
    }
    normed= tpm3(counts= filter_data,len=  g.lenghts)
    
    
   # RPKM <- rpkm(dge, gene.length = g.lenghts)
    
     
  }
  #print(hist(normed))
  return(normed)
  
}



#function from the benchmark paper: 


marker.fc <- function(fit2, log2.threshold = 1, output_name = "markers"){
  
  topTable_RESULTS = limma::topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2.threshold)
  AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
  topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]
  
  if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }
  
  markers <- apply(topTable_RESULTS,1,function(x){
    temp = sort(x)
    ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)
    
  })
  
  topTable_RESULTS = topTable_RESULTS[markers,]
  
  markers <- cbind.data.frame(rownames(topTable_RESULTS),
                              t(apply(topTable_RESULTS, 1, function(x){
                                temp = max(x)
                                if(temp < log2.threshold){
                                  temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
                                } else {
                                  temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
                                } 
                                temp
                              })))
  
  colnames(markers) <- c("gene","log2FC","CT")
  markers$log2FC = as.numeric(as.character(markers$log2FC))
  markers <- markers %>% dplyr::arrange(CT,desc(log2FC)) 
  
  markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
  markers$gene <- as.character(markers$gene)
  markers$CT <- as.character(markers$CT)
  
  #write.table(markers, file = output_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  return(markers)
  
}


# transpose dgC matrix  -----------------------------------------------------------------------

IndexToPointer <- function(j) {
  p <- vector(mode = 'integer', length = max(j) + 1)
  index <- seq.int(from = 2, to = length(x = p))
  for (i in seq_along(along.with = index)) {
    p[index[i]] <- sum(j <= i)
  }
  return(p)
}

#' @param p A vector of sparse matrix pointers
#'
#' @return \code{PointerToIndex}: A vector of column (j) indices
#'
#' @rdname SparsePointers
#'
#' @keywords internal
#'
#' @source \code{PointerToIndex} came from
#' \href{https://stackoverflow.com/questions/20008200/r-constructing-sparse-matrix}{StackOverflow}
#' @author \code{PointerToIndex} was written by
#' \href{https://stackoverflow.com/users/980833/josh-obrien}{Josh O'Brien on StackOverflow}
#'
PointerToIndex <- function(p) {
  dp <- diff(x = p)
  j <- rep.int(x = seq_along(along.with = dp), times = dp)
  return(j)
}

Transpose.dgCMatrix <- function(x, ...) {
  i.order <- order(slot(object = x, name = 'i'))
  return(sparseMatrix(
    i = PointerToIndex(p = slot(object = x, name = 'p'))[i.order],
    p = IndexToPointer(j = slot(object = x, name = 'i') + 1),
    x = slot(object = x, name = 'x')[i.order],
    dims = rev(x = dim(x = x)),
    dimnames = rev(x = dimnames(x = x)),
    giveCsparse = TRUE
  )
  )
}



TPM_normalize_SC= function(seu,
                           local = "remote"){
  
  require(Matrix)
  print(Assays(seu))
  
  #get counts
  data= GetAssayData(seu, slot= "counts")
  
  #load gene lengths=

  if (local == "remote"){
    gene_info= readRDS("/net/data.isilon/ag-saez/bq_jlanzer/ReHeaT2/output/Gene_lengths.rds")
  } else if (local =="local"){
    gene_info = readRDS("~/R-projects/ReHeaT2/output/synced/Gene_lengths.rds")
  }
  
  #  filter genes in data and bring in order
  filter_data = data[rownames(data) %in% gene_info$hgnc_symbol, ]
  g.lenghts= gene_info$length[match(rownames(filter_data), table = gene_info$hgnc_symbol)]
  
  # TPM func
  tpm3 <- function(counts,len) {
    x <- counts/len
    tpm = Transpose.dgCMatrix(Transpose.dgCMatrix(x)*1e6/colSums(x))
    # normed= x*1e6/colSums(x)
    return(tpm)
  }
  
  normed= tpm3(counts= filter_data,len=  g.lenghts)
  
  seu[["TPM"]]= CreateAssayObject(data=normed)
  print(Assays(seu))
  return(seu)
}

