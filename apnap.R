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
library(RRatepol)
theme_set(theme_paleo(10))

#read data----

apnap_counts <- read_csv("data/apnap_counts.csv")

ap <- apnap_counts %>% 
  select(`Pinus sylvestris type`:`Parthenocissus`) %>% 
  mutate(ap_counts = rowSums(.))

nap <- apnap_counts %>% 
  select(!Depth & !`Pinus sylvestris type`:`Parthenocissus`) %>% 
  mutate(nap_counts = rowSums(.))

pollen_concentration_reference <- read_csv("data/pollen_concentration_reference.csv") %>%
  mutate(ap_nap_ratio = ap$ap_counts/nap$nap_counts)

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% apnap_counts$Depth)

tl18_age_raw <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD)

sed_prep <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_sed.csv") %>%
  filter(depth <= max(tl18_age$depth))

loi <- sed_prep %>%  
  select(depth, LOI) %>% 
  rename(count = LOI) %>% 
  mutate(param = "LOI", unit = "percent")

sand <- sed_prep %>%  
  select(depth, Sand) %>% 
  rename(count = Sand) %>% 
  mutate(param = "sand", unit = "number per 1 cc") %>% 
  add_row(loi) %>% 
  mutate(param = factor(param, levels = c("sand", "LOI")))

#percent diagram-------------
conc_coef <- pollen_concentration_reference %>% 
  mutate(coef = lycopodium_ref/indicator)

apnap_long <- apnap_counts %>% 
  pivot_longer(!Depth, names_to = "taxon", values_to = "count")

apnap_sum <- apnap_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count))

apnap_conc <- apnap_sum %>% 
  left_join(conc_coef) %>% 
  mutate(Concentration = coef*count_sum) %>% 
  pivot_longer(!Depth, names_to = "param", values_to = "value") %>% 
  filter(param == "Concentration" | param == "ap_nap_ratio")

apnap_perc <- apnap_long %>%
  left_join(apnap_sum) %>% 
  mutate(rel_abund = count/count_sum*100)

apnap_zero <- apnap_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(apnap_counts$Depth) - 1) #select taxa with at least 2 occurences 

apnap_red <- apnap_perc %>% 
  filter(!taxon %in% apnap_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>% 
  ungroup()

apnap_coniss <- apnap_red %>% 
  mutate(rel_abund_trans = rel_abund/100) %>%
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

apnap_order <- unique(apnap_red$taxon)

apnap_plot <- apnap_red %>%
  mutate(taxon = fct_relevel(taxon, apnap_order)) %>% 
  arrange(taxon) %>% 
  ggplot(aes(x = rel_abund, y = Depth)) +
  geom_col_segsh(color = c(rep("darkgreen", 10*42), rep("darkorange", 8*42))) + 
  geom_lineh(color = "black") +
  facet_abundanceh(vars(taxon)) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#sand plot---------

apnap_sand_plot <- ggplot(sand, aes(x = count, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, linewidth = 0.5, colour = "black")

#Principal curve-----

age_depth <- tibble(Depth = apnap_counts$Depth, Age = tl18_age$ageAD)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

apnap_prc_prep <- apnap_red %>%
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(id_cols = Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth) 

#dalem sqrt() do prc, zeby wieksza wage nadac mniej licznym taksonom. tutaj nie ma to w zasadzie znaczenia, ale liczy sie to bardzo np przy okrzemkach, czy wioslarkach. sqrt() w niewielkim stopniu zmienia przebieg krzywych, za to znaczaco wplywa na roznice miedzy taksonami.
set.seed(12)
apnap_prc <- prcurve(
  sqrt(apnap_prc_prep),
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

apnap_prc# variation explained by prc = 50%

apnap_prc_plot_prep <- apnap_prc_prep %>% 
  mutate(Depth = apnap_counts$Depth, prc_scores = apnap_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

apnap_prcplot <- ggplot(apnap_prc_plot_prep,
                        aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3, scales = "free_y") +
  theme_bw()

apnap_prc_scores <- apnap_counts %>% 
  select(Depth) %>% 
  mutate(value = apnap_prc$lambda, param = rep("prc_score", nrow(apnap_counts)))

write_csv(apnap_prc_scores, "data/apnap_prc_scores.csv")

#rate of change-----

#below is the code necessary to reproduce RoC analysis. The computation takes a lot of time, so the output was saved in a separate file.

# apnap_source_community <- apnap_counts %>% 
#   janitor::clean_names() %>% 
#   select(!depth) %>% 
#   mutate(sample_id = as.character(seq(11, length(apnap_counts$Depth) + 10, by = 1))) %>% 
#   relocate(sample_id)
# 
# apnap_source_age <- apnap_counts %>%
#   rename(depth = Depth) %>% 
#   left_join(tl18_age) %>% 
#   select(depth, ageBP) %>% 
#   rename(age = ageBP) %>% 
#   mutate(sample_id = apnap_source_community$sample_id) %>% 
#   relocate(sample_id)
# 
# set_randomisations <- 10000
# 
# set.seed(12)
# apnap_roc <-
#   RRatepol::estimate_roc(
#     data_source_community = apnap_source_community,
#     data_source_age = apnap_source_age,
#     smooth_method = "age.w",
#     dissimilarity_coefficient = "chisq",
#     working_units = "MW",
#     bin_size = 30, # ~6*median of age difference between the samples
#     number_of_shifts = 5,
#     standardise = TRUE,
#     n_individuals = round(min(ap$ap_counts + nap$nap_counts)), #standardization adjusted to the smallest count sum of apnap
#     rand = set_randomisations, #randomizations are only for taxon standardization
#     use_parallel = FALSE
#   )
# 
# apnap_roc_peaks <-
#   RRatepol::detect_peak_points(
#     data_source = apnap_roc,
#     sel_method = "trend_non_linear"
#   ) %>% 
#   mutate(Peak = as.character(Peak))
# 
# write_csv(apnap_roc_peaks, "data/apnap_roc_peaks.csv")

apnap_roc_peaks <- read_csv("data/apnap_roc_peaks.csv") %>% 
  mutate(ageAD = 1950-Age,
         Depth = approx(tl18_age_raw$ageAD, tl18_age_raw$depth, xout = ageAD)$y)

apnap_roc_2join <- apnap_roc_peaks %>% 
  mutate(value = ROC,
         param = "RoC") %>% 
  select(Depth, param, value)

apnap_roc_peaks_2plot <- apnap_roc_peaks %>% 
  select(Depth, Peak) %>% 
  filter(Peak == TRUE)

#combined PrC, RoC, and pollen concentration plot----------
apnap_adds_plot_prep <- apnap_conc %>% 
  add_row(apnap_prc_scores) %>% 
  add_row(apnap_roc_2join)

apnap_adds_plot <- ggplot(apnap_adds_plot_prep, aes(x = value, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_point(data = filter(apnap_adds_plot_prep,
                           Depth %in%  apnap_roc_peaks_2plot$Depth), color = "green", size = 3) +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
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
  layer_dendrogram(apnap_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, linewidth = 0.5, colour = "black")

#wrap plot---------

apnap_wrapped_plots <- wrap_plots(
  apnap_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  apnap_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  apnap_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 3)
)

ggsave(filename="figures/pollen_wrapped.svg",
       plot = apnap_wrapped_plots,
       device = svg,
       width = 12.5,
       height = 5,
       units = "in")

ggsave(filename="figures/pollen_wrapped.pdf",
       plot = apnap_wrapped_plots,
       device = pdf,
       width = 12.5,
       height = 5,
       units = "in")
