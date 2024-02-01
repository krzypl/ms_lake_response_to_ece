library(rioja)
library(readxl)
library(tidyverse)
library(readxl)
library(tidypaleo)
library(patchwork)
library(analogue)
library(ggvegan)
library(vegan)
library(ggrepel)
library(janitor)
library(RRatepol)
theme_set(theme_paleo(10))

#read data-------

aqps_counts <- read_csv("data/aqps_counts.csv")

apnap_count_sums <- read_csv("data/apnap_counts.csv") %>% 
  select(!c(Depth, `Unidentified pollen - corroded`, `Unidentified pollen - varia`)) %>% 
  transmute(apnap_count_sums = rowSums(.),
            Depth = read_csv("data/apnap_counts.csv")$Depth)

pollen_concentration_reference <- read_csv("data/pollen_concentration_reference.csv")

conc_coef <- pollen_concentration_reference %>% 
  mutate(coef = lycopodium_ref/indicator/volume)

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% aqps_counts$Depth)

tl18_age_raw <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD)

sand <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_sed.csv") %>%
  filter(depth <= max(tl18_age$depth)) %>% 
  select(depth, Sand) %>% 
  rename(count = Sand) %>% 
  mutate(param = "Sand > 0.25 mm", unit = "# cm-3")

#percent diagram-------

aqps_long <- aqps_counts %>% 
  pivot_longer(!Depth, names_to = "taxon", values_to = "count") %>% 
  mutate(taxon = as.factor(taxon))

aqps_sum <- aqps_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count))

aqps_conc <- aqps_sum %>% 
  left_join(conc_coef) %>% 
  mutate(Concentration = coef*count_sum) %>% 
  pivot_longer(!Depth, names_to = "param", values_to = "value") %>% 
  filter(param == "Concentration")

aqps_perc <- aqps_long %>%
  left_join(apnap_count_sums) %>% 
  left_join(aqps_sum) %>% 
  mutate(rel_abund = count/apnap_count_sums * 100)

aqps_zero <- aqps_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(aqps_counts$Depth) - 1) #select taxa with at least 2 occurences 

aqps_aq <- aqps_perc %>% 
  filter(taxon %in% c("Utricularia", "Potamogeton subgen. Eupotamogeton", "Ruppia maritima", "Elodea", "Nuphar", "Utricularia - hair", "Myriophyllum", "Isoetes", "Scenedesmus", "Pediastrum", "Botryococcus", "Tetraedron", "Zygnemataceae", "Dinoflagellata", "HdV-128A", "HdV-128B", "HdV-127")) #select aquatics

aqps_ter <- aqps_perc %>% 
  filter(!taxon %in% unique(aqps_aq$taxon)) #select spores of terrestrial

write_csv(aqps_ter, "data/aqps_ter.csv")

aqps_red <- aqps_aq %>% 
  filter(!taxon %in% aqps_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>%
  ungroup() %>% 
  mutate(taxon = gsub("Potamogeton subgen. Eupotamogeton", "Potamogeton subgen. \n Eupotamogeton", taxon))

write_csv(aqps_red, "data/aqps_red.csv")

aqps_red_trans_prep <- aqps_red %>% 
  pivot_wider(id_cols = Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth)

aqps_red_trans <- decostand(aqps_red_trans_prep, method = "normalize", MARGIN = 2) %>% 
  mutate(Depth = aqps_counts$Depth) %>% 
  pivot_longer(!Depth, names_to = "taxon", values_to = "rel_abund_norm")
  

aqps_coniss <- aqps_red_trans %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_norm,
              trans = identity) %>% 
  nested_chclust_coniss()

aqps_plot <- ggplot(aqps_red, aes(x = rel_abund, y = Depth)) +
  geom_col_segsh() + 
  geom_lineh() +
  facet_geochem_gridh(vars(taxon), 
                      labeller = purrr::partial(label_species, species_facet = "taxon")) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(aqps_coniss, aes(y = Depth), linetype = 2, linewidth = 0.5, colour = "black") +
  rotated_facet_labels(
    angle = 90,
    direction = "x",
    remove_label_background = TRUE
  )

#sand plot-----------

aqps_sand_plot <- ggplot(sand, aes(x = count, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, linewidth = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param), units = c("Sand > 0.25 mm" = "# cm⁻³")) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(aqps_coniss, aes(y = Depth), linetype = 2, linewidth = 0.5, colour = "black") +
  rotated_facet_labels(
    angle = 90,
    direction = "x",
    remove_label_background = TRUE
  )

#Principal curve-------

age_depth <- tibble(Depth = aqps_counts$Depth, Age = tl18_age$ageAD)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

aqps_prc_prep <- aqps_red_trans %>%
  pivot_wider(id_cols = Depth, names_from = taxon, values_from = rel_abund_norm) %>% 
  select(!Depth)

set.seed(12)
aqps_prc <- prcurve(
  aqps_prc_prep,
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

aqps_prc# variation explained by prc = 70%

aqps_prc_plot_prep <- aqps_prc_prep %>%
  mutate(Depth = aqps_counts$Depth, prc_scores = aqps_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund_norm")

aqps_prcplot <- ggplot(aqps_prc_plot_prep,
                          aes(x = prc_scores, y = rel_abund_norm)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3, scales = "free_y") +
  theme_bw()

aqps_prc_scores <- aqps_counts %>% 
  select(Depth) %>% 
  mutate(value = aqps_prc$lambda, param = rep("prc_score", nrow(aqps_counts)))

#rate of change-----
# aqps_source_community_prep <- aqps_counts %>%
#   select(all_of(aqps_aq$taxon)) %>% 
#   janitor::clean_names()
# 
# # weights are added to aqps counts to account for variation in pollen count sums
# aqps_roc_weight <- 
#   rowSums(aqps_source_community_prep)/apnap_count_sums$apnap_count_sums
# 
# aqps_source_community <- as_tibble(aqps_source_community_prep*aqps_roc_weight) %>%
#   mutate(sample_id = as.character(seq(11, length(aqps_counts$Depth) + 10, by = 1))) %>%
#   relocate(sample_id)
# 
# aqps_source_age <- aqps_counts %>%
#   rename(depth = Depth) %>%
#   left_join(tl18_age) %>%
#   select(depth, ageBP) %>%
#   rename(age = ageBP) %>%
#   mutate(sample_id = aqps_source_community$sample_id) %>%
#   relocate(sample_id)
# 
# set_randomisations <- 10000
# 
# set.seed(12)
# aqps_roc <-
#   RRatepol::estimate_roc(
#     data_source_community = aqps_source_community,
#     data_source_age = aqps_source_age,
#     smooth_method = "age.w",
#     dissimilarity_coefficient = "chisq",
#     working_units = "MW", #set to "MW" to apply the "moving window"
#     bin_size = 30,   #~6*median of age difference between the samples
#     number_of_shifts = 5,
#     standardise = TRUE,
#     n_individuals = round(min(rowSums(aqps_source_community[,-1]))), #standardization adjusted to the smallest weighted count sum of aqps
#     rand = set_randomisations, #set number of randomisations
#     use_parallel = FALSE
#   )
# 
# aqps_roc_peaks <-
#   RRatepol::detect_peak_points(
#     data_source = aqps_roc,
#     sel_method = "trend_non_linear"
#   ) %>%
#   mutate(Peak = as.character(Peak))
# 
# write_csv(aqps_roc_peaks, "data/aqps_roc_peaks.csv")

aqps_roc_peaks <- read_csv("data/aqps_roc_peaks.csv") %>% 
  mutate(ageAD = 1950-Age,
         Depth = approx(tl18_age_raw$ageAD, tl18_age_raw$depth, xout = ageAD)$y)

aqps_roc_2join <- aqps_roc_peaks %>% 
  mutate(value = ROC,
         param = "RoC") %>% 
  select(Depth, param, value)

aqps_roc_peaks_2plot <- aqps_roc_peaks %>% 
  select(Depth, Peak) %>% 
  filter(Peak == TRUE)

#combined PrC, RoC, concentration plots----------

aqps_adds_plot_prep <- aqps_conc %>% 
  add_row(aqps_prc_scores) %>%
  add_row(aqps_roc_2join) %>% 
  mutate(param = gsub("prc_score", "PrC score", param))
  
aqps_adds_plot <- ggplot(aqps_adds_plot_prep, aes(x = value, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  geom_point(data = filter(aqps_adds_plot_prep,
                           Depth %in%  aqps_roc_peaks_2plot$Depth), color = "green", size = 3) +
  scale_y_depth_age(
    tl18_adm,
    age_name = "Age (Year AD)",
    breaks = c(breaks = seq(0, 40, by = 5)),
    labels = as.character(c(breaks = seq(0, 40, by = 5))),
    expand = expansion(mult = c(0.02, 0.02)),
    age_breaks = c(2019, seq(1700, 2000, by = 20)),
    age_labels = as.character(c(2019, seq(1700, 2000, by = 20)))
  ) +
  labs(x = NULL, y = "Depth (cm)") +
  facet_geochem_gridh(vars(param), rotate_axis_labels = 90,
                      units = c("Concentration" = "# cm⁻³", 
                                "PrC score" = "70%",
                                "RoC" = NA,
                                "CONISS" = "Total sum \n of squares")) +
  layer_dendrogram(aqps_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(aqps_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black") +
  rotated_facet_labels(
    angle = 90,
    direction = "x",
    remove_label_background = TRUE
  )

#wrap plot---------

aqps_wrapped_plots <- wrap_plots(
  aqps_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  aqps_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  aqps_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 4)
)

ggsave(filename="figures/aqps_wrapped.svg",
       plot = aqps_wrapped_plots,
       device = svg,
       width = 12.5,
       height = 5,
       units = "in")

ggsave(filename="figures/aqps_wrapped.pdf",
       plot = aqps_wrapped_plots,
       device = pdf,
       width = 12.5,
       height = 5,
       units = "in")

ggsave(filename="figures/aqps_wrapped.jpeg",
       plot = aqps_wrapped_plots,
       device = jpeg,
       width = 12.5,
       height = 5,
       units = "in")
#PCA ----
aqps_pca_prep <- decostand(aqps_red_trans_prep, method = "normalize", MARGIN = 2)

aqps_pca <- rda(aqps_pca_prep) 
screeplot(aqps_pca, bstick = TRUE)

aqps_fort <- fortify(aqps_pca, axes = c(1,2), scaling = "sites")

aqps_sites <- aqps_fort[aqps_fort$Score %in% "sites",] %>% 
  mutate(depth = aqps_counts$Depth)
aqps_sp <- aqps_fort[aqps_fort$Score %in% "species",]

aqps_inertcomp <- inertcomp(aqps_pca) #contribution of each taxon to the total inertia
aqps_sel_sp <- tibble(taxon = rownames(aqps_inertcomp), inertcomp = aqps_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
aqps_sp_red <- aqps_sp[which(aqps_sp$Label %in% aqps_sel_sp$taxon), ] #leave 10 taxa with the highest contrubution to the total inertia

aqps_ve_prep <- aqps_pca$CA$eig / aqps_pca$tot.chi * 100
(aqps_PC1_ve <- round(((aqps_ve_prep / sum(aqps_ve_prep))[c(1)]) * 100, digits = 1))#55.4% expl. var.
(aqps_PC2_ve <- round(((aqps_ve_prep / sum(aqps_ve_prep))[c(2)]) * 100, digits = 1))#17.7% expl. var.

aqps_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", aqps_PC2_ve, "%)", sep = ""), x = paste("PC1 (", aqps_PC1_ve, "%)", sep = ""), title = "Borad Pond, TL18-2 AqPS PCA plot") +
  geom_segment(data = aqps_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = aqps_sites, aes(x = PC1, y = PC2, color = depth)) +
  scale_color_viridis_c() +
  geom_point(data = filter(aqps_sites, depth %in% c(5.75, 23.75, 31.75, 35.75)), aes(x = PC1, y = PC2), color = "red") +
  ggrepel::geom_text_repel(data = aqps_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(depth))) +
  ggrepel::geom_text_repel(data = aqps_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = Label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))

ggsave(filename="figures/aqps_pca.svg", 
       plot = aqps_pca_plot, 
       device = svg, 
       width = 7, 
       height = 7, 
       units = "in")

ggsave(filename="figures/aqps_pca.pdf", 
       plot = aqps_pca_plot, 
       device = pdf, 
       width = 7, 
       height = 7, 
       units = "in")

ggsave(filename="figures/aqps_pca.jpeg", 
       plot = aqps_pca_plot, 
       device = jpeg, 
       width = 7, 
       height = 7, 
       units = "in")