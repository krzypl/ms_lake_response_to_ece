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

#read data-------

telmalg_counts <- read_csv("data/telmalg_counts.csv")

apnap_count_sums <- read_csv("data/apnap_counts.csv") %>% 
  select(!Depth) %>% 
  transmute(apnap_count_sums = rowSums(.),
            Depth = read_csv("data/apnap_counts.csv")$Depth)

pollen_concentration_reference <- read_csv("data/pollen_concentration_reference.csv")

conc_coef <- pollen_concentration_reference %>% 
  mutate(coef = lycopodium_ref/indicator)

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% telmalg_counts$Depth)

tl18_age_raw <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD)

sand <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_sed.csv") %>%
  filter(depth <= max(tl18_age$depth)) %>% 
  select(depth, Sand) %>% 
  rename(count = Sand) %>% 
  mutate(param = "sand", unit = "number per 1 cc")

#percent diagram-------

telmalg_long <- telmalg_counts %>% 
  pivot_longer(!Depth, names_to = "taxon", values_to = "count") %>% 
  mutate(taxon = as.factor(taxon))

telmalg_sum <- telmalg_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count))

telmalg_conc <- telmalg_sum %>% 
  left_join(conc_coef) %>% 
  mutate(Concentration = coef*count_sum) %>% 
  pivot_longer(!Depth, names_to = "param", values_to = "value") %>% 
  filter(param == "Concentration")

telmalg_perc <- telmalg_long %>%
  left_join(apnap_count_sums) %>% 
  left_join(telmalg_sum) %>% 
  mutate(rel_abund = count/(apnap_count_sums + count_sum) * 100)

telmalg_zero <- telmalg_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(telmalg_counts$Depth) - 1) #select taxa with at least 2 occurences 

telmalg_red <- telmalg_perc %>% 
  filter(!taxon %in% telmalg_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>%
  ungroup()

write_csv(telmalg_red, "data/telmalg_red.csv")

telmalg_coniss <- telmalg_red %>%
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

telmalg_plot <- ggplot(telmalg_red, aes(x = rel_abund, y = Depth)) +
  geom_col_segsh() + 
  geom_lineh() +
  facet_abundanceh(vars(taxon)) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(telmalg_coniss, aes(y = Depth), linetype = 2, linewidth = 0.5, colour = "black")

#sand plot-----------

telmalg_sand_plot <- ggplot(sand, aes(x = count, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(telmalg_coniss, aes(y = Depth), linetype = 2, linewidth = 0.5, colour = "black")

#Principal curve-------

age_depth <- tibble(Depth = telmalg_counts$Depth, Age = tl18_age$ageAD)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

telmalg_prc_prep <- telmalg_red %>%
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(id_cols = Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth)

set.seed(12)
telmalg_prc <- prcurve(
  sqrt(telmalg_prc_prep),
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

telmalg_prc# variation explained by prc = 72%

telmalg_prc_plot_prep <- telmalg_prc_prep %>% 
  mutate(Depth = telmalg_counts$Depth, prc_scores = telmalg_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

telmalg_prcplot <- ggplot(telmalg_prc_plot_prep,
                          aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3, scales = "free_y") +
  theme_bw()

telmalg_prc_scores <- telmalg_counts %>% 
  select(Depth) %>% 
  mutate(value = telmalg_prc$lambda, param = rep("prc_score", nrow(telmalg_counts)))

#rate of change-----

# telmalg_source_community_prep <- telmalg_counts %>%
#   janitor::clean_names() %>%
#   select(!depth)
# 
# #weights are added to telmalg counts to account for variation in pollen count sums
# telmalg_roc_weight <- telmalg_sum$count_sum/(telmalg_sum$count_sum + apnap_count_sums$apnap_count_sums)
# 
# telmalg_source_community <- as_tibble(telmalg_source_community_prep*telmalg_roc_weight) %>%
#   mutate(sample_id = as.character(seq(11, length(apnap_counts$Depth) + 10, by = 1))) %>%
#   relocate(sample_id)
# 
# telmalg_source_age <- telmalg_counts %>%
#   rename(depth = Depth) %>%
#   left_join(tl18_age) %>%
#   select(depth, ageBP) %>%
#   rename(age = ageBP) %>%
#   mutate(sample_id = telmalg_source_community$sample_id) %>%
#   relocate(sample_id)
# 
# set_randomisations <- 10000
# 
# set.seed(12)
# telmalg_roc <-
#   RRatepol::estimate_roc(
#     data_source_community = telmalg_source_community,
#     data_source_age = telmalg_source_age,
#     smooth_method = "age.w",
#     dissimilarity_coefficient = "chisq",
#     working_units = "MW", # set to "MW" to apply the "moving window"
#     bin_size = 30, #  ~6*median of age difference between the samples
#     number_of_shifts = 5,
#     standardise = TRUE,
#     n_individuals = round(min(rowSums(telmalg_source_community[,-1]))), #standardization adjusted to the smallest weighted count sum of telmalg
#     rand = set_randomisations, # set number of randomisations
#     use_parallel = FALSE
#   )
# 
# telmalg_roc_peaks <-
#   RRatepol::detect_peak_points(
#     data_source = telmalg_roc,
#     sel_method = "trend_non_linear"
#   ) %>%
#   mutate(Peak = as.character(Peak))
# 
# write_csv(telmalg_roc_peaks, "data/telmalg_roc_peaks.csv")

telmalg_roc_peaks <- read_csv("data/telmalg_roc_peaks.csv") %>% 
  mutate(ageAD = 1950-Age,
         Depth = approx(tl18_age_raw$ageAD, tl18_age_raw$depth, xout = ageAD)$y)

telmalg_roc_2join <- telmalg_roc_peaks %>% 
  mutate(value = ROC,
         param = "RoC") %>% 
  select(Depth, param, value)

telmalg_roc_peaks_2plot <- telmalg_roc_peaks %>% 
  select(Depth, Peak) %>% 
  filter(Peak == TRUE)

#combined PrC, RoC, concentration plots----------

telmalg_adds_plot_prep <- telmalg_conc %>% 
  add_row(telmalg_prc_scores) %>%
  add_row(telmalg_roc_2join)
  
telmalg_adds_plot <- ggplot(telmalg_adds_plot_prep, aes(x = value, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  geom_point(data = filter(telmalg_adds_plot_prep,
                           Depth %in%  telmalg_roc_peaks_2plot$Depth), color = "green", size = 3) +
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
  layer_dendrogram(telmalg_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(telmalg_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#wrap plot---------

telmalg_wrapped_plots <- wrap_plots(
  telmalg_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  telmalg_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  telmalg_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 4)
)

ggsave(filename="figures/telmalg_wrapped.svg",
       plot = telmalg_wrapped_plots,
       device = svg,
       width = 12.5,
       height = 5,
       units = "in")

ggsave(filename="figures/telmalg_wrapped.pdf",
       plot = telmalg_wrapped_plots,
       device = pdf,
       width = 12.5,
       height = 5,
       units = "in")

#PCA ----
telmalg_pca_prep <- telmalg_perc %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(id_cols = Depth,
              names_from = taxon,
              values_from = rel_abund) %>% 
  select(!Depth)

telmalg_pca <- rda(sqrt(telmalg_pca_prep)) 
screeplot(telmalg_pca, bstick = TRUE)

telmalg_fort <- fortify(telmalg_pca, axes = c(1,2), scaling = "sites")

telmalg_sites <- telmalg_fort[telmalg_fort$Score %in% "sites",] %>% 
  mutate(depth = telmalg_counts$Depth)
telmalg_sp <- telmalg_fort[telmalg_fort$Score %in% "species",]

telmalg_inertcomp <- inertcomp(telmalg_pca) #contribution of each taxon to the total inertia
telmalg_sel_sp <- tibble(taxon = rownames(telmalg_inertcomp), inertcomp = telmalg_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
telmalg_sp_red <- telmalg_sp[which(telmalg_sp$Label %in% telmalg_sel_sp$taxon), ] #leave 10 taxa with the highest contrubution to the total inertia

telmalg_ve_prep <- telmalg_pca$CA$eig / telmalg_pca$tot.chi * 100
(telmalg_PC1_ve <- round(((telmalg_ve_prep / sum(telmalg_ve_prep))[c(1)]) * 100, digits = 1))#25.9% expl. var.
(telmalg_PC2_ve <- round(((telmalg_ve_prep / sum(telmalg_ve_prep))[c(2)]) * 100, digits = 1))#11.1% expl. var.

telmalg_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", telmalg_PC2_ve, "%)", sep = ""), x = paste("PC1 (", telmalg_PC1_ve, "%)", sep = ""), title = "Borad Pond, TL18-2 telmalg PCA plot") +
  geom_segment(data = telmalg_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = telmalg_sites, aes(x = PC1, y = PC2, color = depth)) +
  scale_color_viridis_c() +
  geom_point(data = filter(telmalg_sites, depth %in% c(5.75, 23.75, 31.75, 35.75)), aes(x = PC1, y = PC2), color = "red") +
  ggrepel::geom_text_repel(data = telmalg_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(depth))) +
  ggrepel::geom_text_repel(data = telmalg_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = telmalg_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))

ggsave(filename="figures/telmalg_pca.svg", 
       plot = telmalg_pca_plot, 
       device = svg, 
       width = 7, 
       height = 7, 
       units = "in")

ggsave(filename="figures/telmalg_pca.pdf", 
       plot = telmalg_pca_plot, 
       device = pdf, 
       width = 7, 
       height = 7, 
       units = "in")
