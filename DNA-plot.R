# read data
library(tidyverse)
RNA_bac <- readxl::read_excel(here::here("B50_陳柏豪.xlsx"), 1)
DNA_bac <- readxl::read_excel(here::here("B50_陳柏豪.xlsx"), 2,
                              guess_max= 3000) # maximum of NA values from the top

# join the tables
RNA_bac_join <-
  RNA_bac %>%
  left_join(DNA_bac, by = "Prokka_ID")

# data wrangling
select_RNA_DNA <-
  RNA_bac_join %>%
  select(Prokka_ID,
         `Log(FPKM+1)/cholesterol`,
         `Log(FPKM+1)/testosterone`,
         `Log(FPKM+1)/estrone`,
         `steroid catabolic function.y`) %>%
  mutate(function_group = ifelse(is.na(`steroid catabolic function.y`), "Y", "N")) %>%
  #mutate(`steroid catabolic function.y` = ifelse(is.na(`steroid catabolic function.y`),
                                                 #"",`steroid catabolic function.y`)) %>%
  mutate(selected_DNA =  ifelse(Prokka_ID %in% c("GMFMDNLD_02935",
                                                 "GMFMDNLD_03019"), # specify data point here
                                Prokka_ID, NA)) %>%
  mutate(`Log(FPKM+1)/estrone - Log(FPKM+1)/cholesterol` =
           `Log(FPKM+1)/estrone`-`Log(FPKM+1)/cholesterol`) %>%
  mutate(`Log(FPKM+1)/estrone - Log(FPKM+1)/testosterone` =
           `Log(FPKM+1)/estrone`-`Log(FPKM+1)/testosterone`)

# make a plot
library(ggplot2)
library(ggforce)
library(ggrepel)
select_RNA_DNA %>%
  ggplot(aes(`Log(FPKM+1)/cholesterol`, `Log(FPKM+1)/testosterone`)) +
  geom_point(aes(color = `steroid catabolic function.y`,
             alpha = function_group),
             size = 2) +
  scale_alpha_discrete(range = c(1, 0.15)) +
  viridis::scale_color_viridis(discrete = TRUE) +
  guides(alpha = FALSE) +
  scale_color_discrete(name="") + # remove legend title
  theme_minimal() +
  theme(legend.position="top")

# use viridis
plot_ct <-
select_RNA_DNA %>%
  ggplot(aes(`Log(FPKM+1)/cholesterol`, `Log(FPKM+1)/testosterone`)) +
  geom_point(color = "darkgray", alpha = 0.5) +
  geom_point(aes(color = `steroid catabolic function.y`),
             size = 2.5,
             alpha = 0.8) +
  #geom_label_repel(aes(label = selected_DNA)) +
                   #box.padding   = 0.35,
                   #point.padding = 0.5,
                   #segment.color = 'grey50') +
  viridis::scale_color_viridis(discrete = TRUE) +
  geom_point(data = subset(select_RNA_DNA, !is.na(selected_DNA)),
             col = "red", stroke = 0.8, shape = 21) +
  guides(alpha = FALSE) +
  theme_minimal() +
  theme(legend.position="top") +
  theme(legend.title=element_blank())

plot_te <-
  select_RNA_DNA %>%
  ggplot(aes(`Log(FPKM+1)/testosterone`, `Log(FPKM+1)/estrone`)) +
  geom_point(color = "darkgray", alpha = 0.5) +
  geom_point(aes(color = `steroid catabolic function.y`),
             size = 2.5,
             alpha = 0.8) +
  #geom_label_repel(aes(label = selected_DNA)) +
  #box.padding   = 0.35,
  #point.padding = 0.5,
  #segment.color = 'grey50') +
  viridis::scale_color_viridis(discrete = TRUE) +
  geom_point(data = subset(select_RNA_DNA, !is.na(selected_DNA)),
             col = "red", stroke = 0.8, shape = 21) +
  guides(alpha = FALSE) +
  theme_minimal() +
  theme(legend.position="none")

plot_ce <-
  select_RNA_DNA %>%
  ggplot(aes(`Log(FPKM+1)/cholesterol`, `Log(FPKM+1)/estrone`)) +
  geom_point(color = "darkgray", alpha = 0.5) +
  geom_point(aes(color = `steroid catabolic function.y`),
             size = 2.5,
             alpha = 0.8) +
  #geom_label_repel(aes(label = selected_DNA)) +
  #box.padding   = 0.35,
  #point.padding = 0.5,
  #segment.color = 'grey50') +
  viridis::scale_color_viridis(discrete = TRUE) +
  geom_point(data = subset(select_RNA_DNA, !is.na(selected_DNA)),
             col = "red", stroke = 0.8, shape = 21) +
  guides(alpha = FALSE) +
  theme_minimal() +
  theme(legend.position="none")

library(cowplot)
plot_grid(plot_ct,
          plot_ce,
          plot_te,
          labels = c('A', 'B', 'C'),
          ncol = 1,
          rel_heights = c(1.3, 1,1))

ggsave(here::here("gene_expression.png"),
       width = 230,
       height = 230,
       dpi = 300,
       units = "mm")

# Plot the difference of gene_expression

diff_gene <-
  select_RNA_DNA %>%
  ggplot(aes(`Log(FPKM+1)/estrone - Log(FPKM+1)/cholesterol`,
             `Log(FPKM+1)/estrone - Log(FPKM+1)/testosterone`))+
  geom_point(color = "darkgray", alpha = 0.5) +
  geom_point(aes(color = `steroid catabolic function.y`),
             size = 2.5,
             alpha = 0.9) +
  #geom_label_repel(aes(label = selected_DNA)) +
  #box.padding   = 0.35,
  #point.padding = 0.5,
  #segment.color = 'grey50') +
  viridis::scale_color_viridis(discrete = TRUE) +
  geom_point(data = subset(select_RNA_DNA, !is.na(selected_DNA)),
             col = "red", stroke = 0.8, shape = 21) +
  guides(alpha = FALSE) +
  theme_minimal() +
  theme(legend.position="top") +
  theme(legend.title=element_blank())

ggsave(here::here("gene_expression_diff.png"),
       width = 200,
       height = 200,
       dpi = 300,
       units = "mm")
