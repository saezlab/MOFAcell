## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-01-24
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## 
library(tidyverse)

pred= readRDS("output/synced/benchmark_hca_to_hf.rds")
pb= readRDS("data/processed/HF_studiespb.rds")

source("analysis/utils_decon.R")

reich_meta= read.csv("data/processed/reichart/Reichart2022_meta_data.csv")[,-1]%>% mutate(study= "Reichart2022")
chaff_meta= read.csv("data/processed/scell_chaffin/Chaffin2022_meta_data.csv")[,-1]%>% mutate(study= "Chaffin2022")

chaff= pb$Chaffin2022$props %>% column_to_rownames("donor_id")%>% t()
reich= pb$Reichart2022$props %>% column_to_rownames("donor_id")%>% t()

cell_dic= data.frame("real"= rownames(chaff), "pred"= c("Ventricular_Cardiomyocyte",
                                              "Endothelial" ,
                                              "Fibroblast"  ,
                                              "Lymphoid",
                                              "Myeloid", 
                                              "Pericytes" ,
                                              "Smooth_muscle_cells")
                     )
                      

chaff= pb$Chaffin2022$props %>% column_to_rownames("donor_id")%>% t()
reich= pb$Reichart2022$props %>% column_to_rownames("donor_id")%>% t()

#map deconvo methods
res= map(names(pred$Chaffin2022), function(x){
  print(x)
  res1= map(names(pred), function(study){
    print(study)
    pred_prop = pred[[study]][[x]]      
    real_prop = pb[[study]]$props %>% column_to_rownames("donor_id")%>% t()
    
    #translate to HCA cell ontology
    rownames(real_prop) = cell_dic[match(rownames(real_prop), cell_dic$real),"pred"]
    p1= evaluate_decon_res(pred_prop, real_prop, "decon")
    
  })
  names(res1)= names(pred)
  return(res1)
})

names(res)= names(pred$Chaffin2022)

# plot cell centric ---------------------------------------------------------------------------


rm.cor= map(res, function(method){
  rm= c(method$Chaffin2022$corr$rmse_total,  method$Reichart2022$rmse$rmse_total)
  cor= c(method$Chaffin2022$corr$cors_total,  method$Reichart2022$corr$cors_total)
  c("rmse"= rm,"corr"=  cor)
})%>% do.call( rbind, .)

plot.df= rm.cor%>%as.data.frame()%>%
  rownames_to_column("method")%>%
  pivot_longer(-method)%>% 
  mutate(study= ifelse(grepl("1", name), "Chaffin2022", "Reichart2022"),
         metric= str_replace_all(name, "1", ""),
         metric= str_replace_all(metric, "2", ""))
plot.df

ggplot(plot.df, aes(x= method, y= value, color= study))+
  geom_point()+
  facet_grid(rows= vars(metric), scales = "free")




# plot patient centric ------------------------------------------------------------------------

unique(reich_meta$disease)
unique(chaff_meta$disease)

meta= rbind(reich_meta[,c("donor_id", "disease", "study")],
      chaff_meta[,c("donor_id", "disease","study")])%>%
  mutate(HF = ifelse(grepl("NF", disease), "NF", "HF"),
         HF = ifelse(grepl("normal", disease) | HF== "NF", "NF", "HF"))


##### correlation
rm.cor= lapply(names(res), function(method){
  enframe( c(res[[method]]$Chaffin2022$pat_cors, 
             res[[method]]$Reichart2022$pat_cors), 
           name = "donor_id"
          )%>%
    mutate(method= method)
})%>%
  do.call( rbind, .)

# check summary:
rm.cor %>% group_by(method)%>% 
  summarise(med= median(value))

rm.cor %>% 
  left_join(meta)%>% 
  ggplot(., aes(x= method, y= value, fill= HF))+
  geom_boxplot()

p.corr= rm.cor %>% 
  left_join(meta)%>% 
  ggplot(., aes(x= method, y= value, fill= study))+
  geom_boxplot()+
  labs(y= "RMSE", 
       x= "", 
       fill = "study")+
  scale_fill_manual(values= c("darkred", "darkgrey", "red", "white"))+
  theme_minimal()+
  theme(axis.text = element_text(size= 11, color= "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p.corr2= rm.cor %>% 
  left_join(meta)%>%
  mutate(study.hf= paste0(study, "_",HF))%>%
  ggplot(., aes(x= method, y= value, fill= study.hf))+
  geom_boxplot()+
  labs(y= "pearson correlation", 
       x= "", 
       fill = "study_disease")+
  scale_fill_manual(values= c("darkred", "darkgrey", "red", "white"))+
  theme_minimal()+
  theme(axis.text = element_text(size= 11, color= "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


##### rmse

rm.cor= lapply(names(res), function(method){
  enframe( c(res[[method]]$Chaffin2022$pat_rmse, 
             res[[method]]$Reichart2022$pat_rmse), 
           name = "donor_id"
  )%>%
    mutate(method= method)
})%>%
  do.call( rbind, .)

# check summary:
rm.cor %>% group_by(method)%>% 
  summarise(med= median(value))

rm.cor %>% 
  left_join(meta)%>% 
  ggplot(., aes(x= method, y= value, fill= HF))+
  geom_boxplot()

p.rmse1= rm.cor %>% 
  left_join(meta)%>% 
  ggplot(., aes(x= method, y= value, fill= study))+
  geom_boxplot()+
  labs(y= "RMSE", 
       x= "", 
       fill = "study")+
  scale_fill_manual(values= c("darkred", "darkgrey", "red", "white"))+
  theme_minimal()+
  theme(axis.text = element_text(size= 11, color= "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p.rmse2= rm.cor %>% 
  left_join(meta)%>%
  mutate(study.hf= paste0(study, "_",HF))%>%
  ggplot(., aes(x= method, y= value, fill= study.hf))+
  geom_boxplot()+
  labs(y= "RMSE", 
       x= "", 
       fill = "study_disease")+
  scale_fill_manual(values= c("darkred", "darkgrey", "red", "white"))+
  theme_minimal()+
  theme(axis.text = element_text(size= 11, color= "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

p_bench= cowplot::plot_grid(p.rmse1, p.corr)
p_bench
pdf("output/figures/mofacell_paper/benchmark_deconvo.pdf",
    width= 9 ,
    height= 3)
p_bench
dev.off()

#save predicted props for contrast analysis
saveRDS(list("Reichart2022"= res$SCDC$Reichart2022$props,
             "Chaffin2022"= res$SCDC$Chaffin2022$props), 
        "output/deconvo/deconvo_proportions.rds")
