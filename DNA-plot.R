# read data
library(tidyverse)
RNA_bac <- readxl::read_excel(here::here("B50_陳柏豪.xlsx"), 1)
DNA_bac <- readxl::read_excel(here::here("B50_陳柏豪.xlsx"), 2, guess_max= 3000)

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
  mutate(alpha_group = ifelse(is.na(`steroid catabolic function.y`), "Y", "N"))

# make a plot
library(ggplot2)
library(ggforce)
select_RNA_DNA %>%
  ggplot(aes(`Log(FPKM+1)/cholesterol`, `Log(FPKM+1)/testosterone`)) +
  geom_point(aes(color = `steroid catabolic function.y`,
             alpha = alpha_group),
             size = 2) +
  scale_alpha_discrete(range = c(1, 0.15)) +
  viridis::scale_color_viridis(discrete = TRUE) +
  guides(alpha = FALSE) +
  scale_color_discrete(name="") + # remove legend title
  theme_minimal() +
  theme(legend.position="top")

