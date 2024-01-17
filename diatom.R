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

apnap_counts <- read_csv("data/apnap_counts.csv")

diatom_counts <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_diatom.csv")%>% 
  dplyr::filter(depth <= max(apnap_counts$Depth))

colnames(diatom_counts)<-gsub("[(M)]","",colnames(diatom_counts))
colnames(diatom_counts)<-gsub("B ","",colnames(diatom_counts))
colnames(diatom_counts)<-gsub("E ","",colnames(diatom_counts))
colnames(diatom_counts)<-gsub("F ","",colnames(diatom_counts))
colnames(diatom_counts)<-gsub("I ","",colnames(diatom_counts))

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% diatom_counts$depth)

tl18_age_raw <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD)

sand <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_sed.csv") %>%
  filter(depth <= max(tl18_age$depth)) %>% 
  select(depth, Sand) %>% 
  rename(count = Sand) %>% 
  mutate(param = "sand", unit = "number per 1 cc")

#percent diagram-------------

diatom_long <- diatom_counts %>% 
  pivot_longer(!depth, names_to = "taxon", values_to = "count")

diatom_sum <- diatom_long %>% 
  group_by(depth) %>% 
  summarise(count_sum = sum(count))

diatom_perc <- diatom_long %>%
  left_join(diatom_sum) %>% 
  mutate(rel_abund = count/count_sum*100)

diatom_zero <- diatom_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(diatom_counts$depth) - 1) #select taxa with at least 2 occurences 

diatom_red <- diatom_perc %>% 
  filter(!taxon %in% diatom_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>% 
  ungroup()

write_csv(diatom_red, "data/diatom_red.csv")

diatom_coniss <- diatom_red %>%
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss() #single CONISS zone

diatom_order <- unique(diatom_red$taxon)

diatom_plot <- diatom_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, diatom_order)) %>% 
  ggplot(aes(x = rel_abund, y = depth)) +
  geom_col_segsh() + 
  geom_lineh() +
  facet_abundanceh(vars(taxon)) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02)))

#sand plot-------

diatom_sand_plot <- ggplot(sand, aes(x = count, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)")

#Principal curve---------

age_depth <- tibble(Depth = diatom_counts$depth, Age = tl18_age$ageAD)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

diatom_prc_prep <- diatom_red %>% 
  select(depth, taxon, rel_abund) %>% 
  pivot_wider(id_cols = depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!depth) 

set.seed(12)
diatom_prc <- prcurve(
  sqrt(diatom_prc_prep),
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

diatom_prc# variation explained by prc = 39%

diatom_prc_plot_prep <- diatom_prc_prep %>% 
  mutate(Depth = diatom_counts$depth, prc_scores = diatom_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

diatom_prcplot <- ggplot(diatom_prc_plot_prep,
                         aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3, scales = "free_y") +
  theme_bw()

diatom_prc_scores <- diatom_counts %>% 
  select(depth) %>% 
  mutate(value = diatom_prc$lambda, param = rep("prc_score", nrow(diatom_counts)))


#rate of change-----------

# diatom_source_community <- diatom_counts %>% 
#   janitor::clean_names() %>% 
#   select(!depth) %>% 
#   mutate(sample_id = as.character(seq(11, length(apnap_counts$Depth) + 10, by = 1))) %>% 
#   relocate(sample_id)
# 
# diatom_source_age <- diatom_counts %>%
#   left_join(tl18_age) %>% 
#   select(depth, ageBP) %>% 
#   rename(age = ageBP) %>% 
#   mutate(sample_id = diatom_source_community$sample_id) %>% 
#   relocate(sample_id)
# 
# set_randomisations <- 10000
# 
# set.seed(12)
# diatom_roc <-
#   RRatepol::estimate_roc(
#     data_source_community = diatom_source_community,
#     data_source_age = diatom_source_age,
#     smooth_method = "age.w",
#     dissimilarity_coefficient = "chisq",
#     working_units = "MW", # set to "MW" to apply the "moving window"
#     bin_size = 30, # ~6*median of age difference between the samples
#     number_of_shifts = 5,
#     standardise = FALSE,
#     rand = set_randomisations, # set number of randomisations
#     use_parallel = FALSE
#   )
# 
# diatom_roc_peaks <-
#   RRatepol::detect_peak_points(
#     data_source = diatom_roc,
#     sel_method = "trend_non_linear"
#   ) %>%
#   mutate(Peak = as.character(Peak))
# 
# write_csv(diatom_roc_peaks, "data/diatom_roc_peaks.csv")

diatom_roc_peaks <- read_csv("data/diatom_roc_peaks.csv") %>% 
  mutate(ageAD = 1950-Age,
         Depth = approx(tl18_age_raw$ageAD, tl18_age_raw$depth, xout = ageAD)$y)

diatom_roc_2join <- diatom_roc_peaks %>% 
  mutate(value = ROC,
         param = "RoC") %>% 
  select(Depth, param, value) %>% 
  rename(depth = Depth)

diatom_roc_peaks_2plot <- diatom_roc_peaks %>% 
  select(Depth, Peak) %>% 
  filter(Peak == TRUE)

#combined PrC, RoC, concentration plots----------

diatom_adds_plot_prep <- diatom_prc_scores %>%
  add_row(diatom_roc_2join)

diatom_adds_plot <- ggplot(diatom_adds_plot_prep, aes(x = value, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  geom_point(data = filter(diatom_adds_plot_prep,
                           depth %in%  diatom_roc_peaks_2plot$Depth), color = "green", size = 3) +
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
  layer_dendrogram(diatom_coniss, aes(y = depth), param = "CONISS")
  
#wrap plot---------

diatom_wrapped_plots <- wrap_plots(
  diatom_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  diatom_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  diatom_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 3)
)

ggsave(filename="figures/diatom_wrapped.svg",
         plot = diatom_wrapped_plots,
         device = svg,
         width = 12.5,
         height = 5,
         units = "in")
  
ggsave(filename="figures/diatom_wrapped.pdf",
         plot = diatom_wrapped_plots,
         device = pdf,
         width = 12.5,
         height = 5,
         units = "in")