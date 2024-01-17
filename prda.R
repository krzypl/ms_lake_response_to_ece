library(rioja)
library(readxl)
library(tidyverse)
library(readxl)
library(tidypaleo)
library(patchwork)
library(analogue)
library(ggvegan)
library(ggrepel)
theme_set(theme_paleo(10))

#prepere data---------

apnap_prc_scores <- read_csv("data/apnap_prc_scores.csv")

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% apnap_prc_scores$Depth)

expl_var <- tibble(Depth = apnap_prc_scores$Depth, sand_peak = 0)
expl_var[which(expl_var$Depth %in% c(5.75, 23.75, 31.75, 35.75)),2] <- 100

expl_var <- expl_var %>% 
  mutate(apnap_prc = read_csv("data/apnap_prc_scores.csv")$value,
         Age = tl18_age$ageAD)

sand_peaks <- expl_var %>% 
  filter(sand_peak == 100) %>% 
  mutate(time_elapsed_since_event = 0)

time_elapsed_ev1 <- expl_var %>% 
  filter(sand_peak == 0 & Depth < sand_peaks$Depth[1]) %>% 
  mutate(time_elapsed_since_event = -(sand_peaks$Age[1] - Age))

time_elapsed_ev2 <- expl_var %>% 
  filter(sand_peak == 0 & Depth < sand_peaks$Depth[2] & Depth > sand_peaks$Depth[1] ) %>% 
  mutate(time_elapsed_since_event = -(sand_peaks$Age[2] - Age))

time_elapsed_ev3 <- expl_var %>% 
  filter(sand_peak == 0 & Depth < sand_peaks$Depth[3] & Depth > sand_peaks$Depth[2]) %>% 
  mutate(time_elapsed_since_event = -(sand_peaks$Age[3] - Age))

time_elapsed_ev4 <- expl_var %>% 
  filter(sand_peak == 0 & Depth < sand_peaks$Depth[4] & Depth > sand_peaks$Depth[3]) %>% 
  mutate(time_elapsed_since_event = -(sand_peaks$Age[4] - Age))

time_elapsed_ev5 <- expl_var %>% 
  filter(sand_peak == 0 & Depth > sand_peaks$Depth[4]) %>% 
  mutate(time_elapsed_since_event = -(1678 - Age))

ece_effect_prep <- sand_peaks %>% 
  add_row(time_elapsed_ev1) %>% 
  add_row(time_elapsed_ev2) %>% 
  add_row(time_elapsed_ev3) %>%
  add_row(time_elapsed_ev4) %>%
  add_row(time_elapsed_ev5)

expl_var_full <- expl_var %>% 
  left_join(ece_effect_prep) %>% 
  mutate(ece_effect = 100*exp(1)^(-0.5*time_elapsed_since_event))

ggplot(expl_var_full, aes(x = Age, y = ece_effect)) + 
  geom_point() +
  geom_line()

#diatom rda -----
diatom_wide <- read_csv("data/diatom_red.csv") %>% 
  pivot_wider(id_cols = !count & !count_sum, names_from = "taxon", values_from = "rel_abund") %>% 
  select(!depth)

diatom_rda <- rda(sqrt(diatom_wide) ~ ece_effect + Condition(Age, apnap_prc), data = expl_var_full)

set.seed(12)
(diatom_R2adj <- RsquareAdj(diatom_rda)$adj.r.squared)#6.9% 

set.seed(12)
(diatom_anova <- anova(diatom_rda, permutations = how(nperm = 999)))#significance level 0.001***

diatom_rda_all <- rda(sqrt(diatom_wide) ~ ece_effect + Age + apnap_prc, data = expl_var_full)

set.seed(12)
(diatom_R2adj_all <- RsquareAdj(diatom_rda_all)$adj.r.squared)#18.4% 

set.seed(12)
(diatom_anova_all <- anova(diatom_rda_all, permutations = how(nperm = 999)))#significance level 0.001***

diatom_rda_age <- rda(sqrt(diatom_wide) ~ Age + Condition(apnap_prc, ece_effect), data = expl_var_full)

set.seed(12)
(diatom_R2adj_age <- RsquareAdj(diatom_rda_age)$adj.r.squared)#1.8% 

set.seed(12)
(diatom_anova_age <- anova(diatom_rda_age, permutations = how(nperm = 999)))#not significant

diatom_rda_apnap <- rda(sqrt(diatom_wide) ~ apnap_prc + Condition(Age, ece_effect), data = expl_var_full)

set.seed(12)
(diatom_R2adj_apnap <- RsquareAdj(diatom_rda_apnap)$adj.r.squared)#2.8% 

set.seed(12)
(diatom_anova_apnap <- anova(diatom_rda_apnap, permutations = how(nperm = 999)))#not significant

diatom_rda_tprep <- tibble("ECE_xage_xapnap" = c(diatom_R2adj, diatom_anova$`Pr(>F)`[[1]]),
                           "apnap_xECE_xage" = c(diatom_R2adj_apnap, diatom_anova_apnap$`Pr(>F)`[[1]]),
                           "age_xECE_xapnap" = c(diatom_R2adj_age, diatom_anova_age$`Pr(>F)`[[1]]),
                           "all" = c(diatom_R2adj_all, diatom_anova_all$`Pr(>F)`[[1]]),
                           "unexplained" = c(1 - diatom_R2adj_all, NA),
                           "proxy" = c("diatom", "diatom"))

#telmalg rda -----
telmalg_wide <- read_csv("data/telmalg_red.csv") %>% 
  pivot_wider(id_cols = !count & !count_sum, names_from = "taxon", values_from = "rel_abund") %>% 
  select(!c(Depth, apnap_count_sums))

telmalg_rda <- rda(sqrt(telmalg_wide) ~ ece_effect + Condition(Age, apnap_prc), data = expl_var_full)

set.seed(12)
(telmalg_R2adj <- RsquareAdj(telmalg_rda)$adj.r.squared)#-1.1% 

set.seed(12)
(telmalg_anova <- anova(telmalg_rda, permutations = how(nperm = 999)))#not significant

telmalg_rda_all <- rda(sqrt(telmalg_wide) ~ ece_effect + Age + apnap_prc, data = expl_var_full)

set.seed(12)
(telmalg_R2adj_all <- RsquareAdj(telmalg_rda_all)$adj.r.squared)#8.9% 

set.seed(12)
(telmalg_anova_all <- anova(telmalg_rda_all, permutations = how(nperm = 999)))#significant at p 0.05

telmalg_rda_age <- rda(sqrt(telmalg_wide) ~ Age + Condition(apnap_prc, ece_effect), data = expl_var_full)

set.seed(12)
(telmalg_R2adj_age <- RsquareAdj(telmalg_rda_age)$adj.r.squared)#0% 

set.seed(12)
(telmalg_anova_age <- anova(telmalg_rda_age, permutations = how(nperm = 999)))#not significant

telmalg_rda_apnap <- rda(sqrt(telmalg_wide) ~ apnap_prc + Condition(Age, ece_effect), data = expl_var_full)

set.seed(12)
(telmalg_R2adj_apnap <- RsquareAdj(telmalg_rda_apnap)$adj.r.squared)#-1.7% 

set.seed(12)
(telmalg_anova_apnap <- anova(telmalg_rda_apnap, permutations = how(nperm = 999)))#not significant

telmalg_rda_tprep <- tibble("ECE_xage_xapnap" = c(telmalg_R2adj, telmalg_anova$`Pr(>F)`[[1]]),
                            "apnap_xECE_xage" = c(telmalg_R2adj_apnap, telmalg_anova_apnap$`Pr(>F)`[[1]]),
                            "age_xECE_xapnap" = c(telmalg_R2adj_age, telmalg_anova_age$`Pr(>F)`[[1]]),
                            "all" = c(telmalg_R2adj_all, telmalg_anova_all$`Pr(>F)`[[1]]),
                            "unexplained" = c(1 - telmalg_R2adj_all, NA),
                            "proxy" = c("telmalg", "telmalg"))

#clado rda ------
clado_wide <- read_csv("data/clado_red.csv") %>% 
  pivot_wider(id_cols = !count & !count_sum & !volume, names_from = "taxon", values_from = "rel_abund") %>% 
  select(!Depth)

clado_rda <- rda(sqrt(clado_wide) ~ ece_effect + Condition(Age, apnap_prc), data = expl_var_full)

set.seed(12)
(clado_R2adj <- RsquareAdj(clado_rda)$adj.r.squared)#-0.8% 

set.seed(12)
(clado_anova <- anova(clado_rda, permutations = how(nperm = 999)))#not significant

clado_rda_all <- rda(sqrt(clado_wide) ~ ece_effect + Age + apnap_prc, data = expl_var_full)

set.seed(12)
(clado_R2adj_all <- RsquareAdj(clado_rda_all)$adj.r.squared)#8.8% 

set.seed(12)
(clado_anova_all <- anova(clado_rda_all, permutations = how(nperm = 999)))#significance level 0.001**

clado_rda_age <- rda(sqrt(clado_wide) ~ Age + Condition(apnap_prc, ece_effect), data = expl_var_full)

set.seed(12)
(clado_R2adj_age <- RsquareAdj(clado_rda_age)$adj.r.squared)#2.7% 

set.seed(12)
(clado_anova_age <- anova(clado_rda_age, permutations = how(nperm = 999)))#significance level 0.05*

clado_rda_apnap <- rda(sqrt(clado_wide) ~ apnap_prc + Condition(Age, ece_effect), data = expl_var_full)

set.seed(12)
(clado_R2adj_apnap <- RsquareAdj(clado_rda_apnap)$adj.r.squared)#1.5% 

set.seed(12)
(clado_anova_apnap <- anova(clado_rda_apnap, permutations = how(nperm = 999)))#not significant

clado_rda_tprep <- tibble("ECE_xage_xapnap" = c(clado_R2adj, clado_anova$`Pr(>F)`[[1]]),
                          "apnap_xECE_xage" = c(clado_R2adj_apnap, clado_anova_apnap$`Pr(>F)`[[1]]),
                          "age_xECE_xapnap" = c(clado_R2adj_age, clado_anova_age$`Pr(>F)`[[1]]),
                          "all" = c(clado_R2adj_all, clado_anova_all$`Pr(>F)`[[1]]),
                          "unexplained" = c(1 - clado_R2adj_all, NA),
                          "proxy" = c("clado", "clado"))

#tabele RDA ----
rda_table <- diatom_rda_tprep[1,] %>% 
  add_row(telmalg_rda_tprep[1,]) %>% 
  add_row(clado_rda_tprep[1,]) %>% 
  pivot_longer(!proxy, 
               names_to = "category", values_to = "proportion_of_v_expl") %>% 
  mutate(proxy = factor(proxy, levels = c("telmalg","diatom", "clado")),
         category = factor(category, levels = c("age_xECE_xapnap", "apnap_xECE_xage", "ECE_xage_xapnap", "all", "unexplained")))

rda_table_sign <- diatom_rda_tprep[2,] %>% 
  add_row(telmalg_rda_tprep[2,]) %>% 
  add_row(clado_rda_tprep[2,]) %>% 
  pivot_longer(!proxy, 
               names_to = "category", values_to = "pval") %>% 
  mutate(proxy = factor(proxy, levels = c("telmalg","diatom", "clado")),
         category = factor(category, levels = c("age_xECE_xapnap", "apnap_xECE_xage", "ECE_xage_xapnap", "all", "unexplained"))) %>% 
  filter(!category == "unexplained") %>% 
  mutate(stars = cut(pval, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")))

rda_prop_var_plot <- ggplot(rda_table, aes(x = proxy, y = proportion_of_v_expl, fill = category)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(labels = c("Age", "Vegetation shifts (VS)", "Extreme coastal events (ECE)",
                                  "Age + VS + ECE", "Unexplained")) +
  labs(x = "", y = "Proportion of variance explained")#+

rda_pval_plot <- ggplot(rda_table_sign, aes(x = proxy, y = category, color = "white")) +
  geom_tile() + 
  scale_y_discrete(limits=rev) +
  geom_text(aes(label=stars), color="white", size=5)

rda_wrap_plot <- wrap_plots(
  rda_prop_var_plot,
  rda_pval_plot,
  nrow = 1
)

ggsave(filename="figures/rda_wrap_plot.svg", 
       plot = rda_wrap_plot, 
       device = svg, 
       width = 11, 
       height = 5, 
       units = "in")

ggsave(filename="figures/rda_wrap_plot.pdf", 
       plot = rda_wrap_plot, 
       device = pdf, 
       width = 11, 
       height = 5, 
       units = "in")
