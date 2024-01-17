library(rioja)
library(readxl)
library(tidyverse)
library(readxl)
library(tidypaleo)
library(patchwork)
library(analogue)
library(ggvegan)
library(ggrepel)
library(janitor)
theme_set(theme_paleo(10))

#read data-------

chiro_counts <- read_csv("data/chiro_counts.csv")

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% chiro_counts$Depth)

sand <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_sed.csv") %>%
  filter(depth <= max(tl18_age$depth)) %>% 
  select(depth, Sand) %>% 
  rename(count = Sand) %>% 
  mutate(param = "sand", unit = "number per 1 cc")

#percent diagram-------------

chiro_long <- chiro_counts %>% 
  pivot_longer(!Depth & !volume, names_to = "taxon", values_to = "count")

chiro_sum <- chiro_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count)) %>%
  ungroup() %>% 
  mutate(volume = chiro_counts$volume)

chiro_conc <- chiro_sum %>% 
  mutate(Concentration = count_sum/volume) %>% 
  pivot_longer(!Depth, names_to = "param", values_to = "value") %>% 
  filter(param != "volume")

chiro_perc <- chiro_long %>%
  left_join(chiro_sum) %>% 
  mutate(rel_abund = count/count_sum*100)

chiro_zero <- chiro_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(chiro_counts$Depth) - 1) #select taxa with at least 2 occurences 

chiro_red <- chiro_perc %>% 
  filter(!taxon %in% chiro_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>% 
  ungroup()

chiro_coniss <- chiro_red %>%
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

chiro_order <- unique(chiro_red$taxon)

chiro_plot <- chiro_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, chiro_order)) %>% 
  ggplot(aes(x = rel_abund, y = Depth)) +
  geom_col_segsh() + 
  geom_lineh() +
  facet_abundanceh(vars(taxon)) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(chiro_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#sand plot------

chiro_sand_plot <- ggplot(sand, aes(x = count, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(chiro_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#Concentration and count sum plots-------

age_depth <- tibble(Depth = chiro_counts$Depth, Age = tl18_age$ageAD)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

chiro_adds_plot <- chiro_conc %>% 
  ggplot(aes(x = value, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_depth_age(
    tl18_adm,
    age_name = "Age (Year AD)",
    breaks = c(breaks = seq(0, 40, by = 5)),
    labels = as.character(c(breaks = seq(0, 40, by = 5))),
    expand = expansion(mult = c(0.02, 0.02)),
    age_breaks = c(2019, seq(1700, 2000, by = 20)),
    age_labels = as.character(c(2019, seq(1700, 2000, by = 20)))
  ) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_dendrogram(chiro_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(chiro_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#wrap plot-----

chiro_wrapped_plots <- wrap_plots(
  chiro_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  chiro_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  chiro_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 2)
)

ggsave(filename="figures/chiro_wrapped.svg",
       plot = chiro_wrapped_plots,
       device = svg,
       width = 12.5,
       height = 5,
       units = "in")

ggsave(filename="figures/chiro_wrapped.pdf",
       plot = chiro_wrapped_plots,
       device = pdf,
       width = 12.5,
       height = 5,
       units = "in")

