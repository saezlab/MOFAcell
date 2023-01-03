# Copyright (c) [2022] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de
# Using a list of 20 colors provided in:
# https://sashamaps.net/docs/resources/20-colors/
# To create a useful map to combine states from MOFA cell

library(tidyverse)
library(Seurat)

# # The last 3 colors are ignored 
# colors <- c(c(230, 25, 75), 
#             c(60, 180, 75), 
#             c(255, 225, 25), 
#             c(0, 130, 200), 
#             c(245, 130, 48), 
#             c(145, 30, 180), 
#             c(70, 240, 240), 
#             c(240, 50, 230), 
#             c(210, 245, 60), 
#             c(250, 190, 212), 
#             c(0, 128, 128), 
#             c(220, 190, 255), 
#             c(170, 110, 40), 
#             c(255, 250, 200), 
#             c(128, 0, 0), 
#             c(170, 255, 195), 
#             c(128, 128, 0), 
#             c(255, 215, 180), 
#             c(0, 0, 128))
# 
# color_names <- c("red", "green", "yellow", "blue", 
#                  "orange", "purple", "cyan", "magenta", 
#                  "lime", "pink", "teal", "lavender", 
#                  "brown", "beige", "maroon", "mint", 
#                  "olive", "apricot", "navy")


# The last 3 colors are ignored 
colors <- c(c(255, 0, 0),
            c(0, 255, 0),
            c(0, 0, 255))

color_names <- c("red", "green",  "blue")

RGB_mat <- matrix(colors, nrow = 3)
rownames(RGB_mat) <- c("R","G","B")
colnames(RGB_mat) <- color_names

# For a given matrix of at most 17 variables create an RGB combination
matrix_to_hex <- function(spot_matrix, colors = NULL) {
  
  if(is.null(colors)) {
  
  color_labs <- colnames(RGB_mat[,1:nrow(spot_matrix)])
  
  
  } else {
    
    color_labs <- colors[colors %in% colnames(RGB_mat)]
    
    if(length(colors) != length(color_labs)) {
      
      print("You gave wrong labels, try again")
      
      return(NULL)
    }
    
  }
  
  if(length(color_labs) > ncol(RGB_mat)) {
    
    print("we can not use more than 17 variables (rows)")
    return(NULL)
    
  } else{
    
    # Trim RGB mat
    RGB_usemat <- RGB_mat[,color_labs, drop = F]
    colnames(RGB_usemat) <- rownames(spot_matrix)
    # Get the annotation for these colors
    RGB_use_hex <- rgb(red = t(RGB_usemat), 
                       names = colnames(RGB_usemat),
                       maxColorValue = 255)
    
    # Do the linear combination
    entry_RGB_vals <- RGB_usemat %*% spot_matrix %>%
      t()
    
    entry_RGB_hex <- rgb(red = entry_RGB_vals, 
                         names = rownames(entry_RGB_vals),
                         maxColorValue = 255)
    
    
    return(list("legend_key" = enframe(RGB_use_hex,
                                       value = "hex"),
                "data_hex" = enframe(entry_RGB_hex,
                                     name = "sample_id",
                                     value = "hex")))
  }
  
}

# Generate legends for paper

legend_matrix <- expand.grid(seq(0,1,0.2), seq(0,1,0.2))

legend_matrix_hex <- matrix_to_hex(spot_matrix = legend_matrix %>%
                                     as.matrix() %>%
                                     t(),
                                   colors = c("red","blue"))

legend_matrix$color <-  legend_matrix_hex$data_hex$hex

legend_plt <- ggplot(legend_matrix, aes(x = Var1, y = Var2, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = legend_matrix$color) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text = element_text(size =12)) +
  coord_equal() +
  xlab("") +
  ylab("")

pdf(height = 2, width = 2, "./results/aes/paintR_legend_redblue.pdf")

plot(legend_plt)

dev.off()

# red green

legend_matrix <- expand.grid(seq(0,1,0.2), seq(0,1,0.2))


legend_matrix_hex <- matrix_to_hex(spot_matrix = legend_matrix %>%
                                     as.matrix() %>%
                                     t(),
                                   colors = c("red","green"))

legend_matrix$color <-  legend_matrix_hex$data_hex$hex

legend_plt <- ggplot(legend_matrix, aes(x = Var1, y = Var2, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = legend_matrix$color) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text = element_text(size =12)) +
  coord_equal() +
  xlab("") +
  ylab("")

pdf(height = 2, width = 2, "./results/aes/paintR_legend_redgreen.pdf")

plot(legend_plt)

dev.off()

# blue green

legend_matrix <- expand.grid(seq(0,1,0.2), seq(0,1,0.2))


legend_matrix_hex <- matrix_to_hex(spot_matrix = legend_matrix %>%
                                     as.matrix() %>%
                                     t(),
                                   colors = c("blue","green"))

legend_matrix$color <-  legend_matrix_hex$data_hex$hex

legend_plt <- ggplot(legend_matrix, aes(x = Var1, y = Var2, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = legend_matrix$color) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text = element_text(size =12)) +
  coord_equal() +
  xlab("") +
  ylab("")

pdf(height = 2, width = 2, "./results/aes/paintR_legend_bluegreen.pdf")

plot(legend_plt)

dev.off()


