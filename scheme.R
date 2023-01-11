library(tidyverse)

view1 <- matrix(rnorm(140),14) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)

view2 <- matrix(rnorm(120),12) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)

view3 <- matrix(rnorm(100),10) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)


pdf("./view_schemes.pdf", height = 4, width = 12)

view1 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2()

view2 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2(low = "darkgreen", high = "orange")

view3 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2(low = "purple", high = "black")

dev.off()

# Loadings

view1 <- matrix(rnorm(120),20) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)

view2 <- matrix(rnorm(72),12) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)

view3 <- matrix(rnorm(60),10) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)



pdf("./view_loadings.pdf", height = 4, width = 12)

view1 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2(low = "red", high = "red")

view2 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2(low = "red", high = "red")

view3 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2(low = "red", high = "red")

dev.off()


facts <- matrix(rnorm(60),10) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample)

pdf("./factor_scores.pdf", height = 4, width = 12)

view3 %>%
  ggplot(aes(x = sample, y = name, fill = value)) +
  geom_tile(color = "black") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")+
  coord_equal() +
  scale_fill_gradient2(low = "black", high = "black")

dev.off()






















