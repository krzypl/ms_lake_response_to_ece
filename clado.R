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

clado_counts <- read_csv("data/clado_counts.csv")

tl18_age <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD) %>% 
  filter(depth %in% clado_counts$Depth)

tl18_age_raw <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_age_ad.csv") %>% 
  mutate(ageBP = 1950 - ageAD)

sand <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_sed.csv") %>%
  filter(depth <= max(tl18_age$depth)) %>% 
  select(depth, Sand) %>% 
  rename(count = Sand) %>% 
  mutate(param = "Sand > 0.25 mm", unit = "number per 1 cc")

#percent diagram-------------

clado_long <- clado_counts %>% 
  pivot_longer(!Depth & !volume, names_to = "taxon", values_to = "count") %>% 
  mutate(taxon = gsub("Chydorus sphaericus", "Chydorus biovatus/brevilabris", taxon))

clado_sum <- clado_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count)) %>%
  ungroup() %>% 
  mutate(volume = clado_counts$volume)

clado_conc <- clado_sum %>% 
  mutate(Concentration = count_sum/volume) %>% 
  pivot_longer(!Depth, names_to = "param", values_to = "value") %>% 
  filter(param != "volume")

clado_perc <- clado_long %>%
  left_join(clado_sum) %>% 
  mutate(rel_abund = count/count_sum*100)

clado_zero <- clado_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(clado_counts$Depth) - 1) #select taxa with at least 2 occurences 

clado_red <- clado_perc %>% 
  filter(!taxon %in% clado_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>% 
  ungroup()

write_csv(clado_red, "data/clado_red.csv")

clado_coniss <- clado_red %>%
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

clado_order <- unique(clado_red$taxon)

clado_plot <- clado_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, clado_order)) %>% 
  ggplot(aes(x = rel_abund, y = Depth)) +
  geom_col_segsh() + 
  geom_lineh() +
  facet_abundanceh(vars(taxon), rotate_facet_labels = 90) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(clado_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#sand plot---------

clado_sand_plot <- ggplot(sand, aes(x = count, y = depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param),
                      units = c("Sand > 0.25 mm" = "# cm⁻³")) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(clado_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")+
  rotated_facet_labels(
    angle = 90,
    direction = "x",
    remove_label_background = TRUE
  )

#Principal curve------

age_depth <- tibble(Depth = clado_counts$Depth, Age = tl18_age$ageAD)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

clado_prc_prep <- clado_red %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth) 

set.seed(12)
clado_prc <- prcurve(
  sqrt(clado_prc_prep),
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

clado_prc# variation explained by prc = 21%

clado_prc_plot_prep <- clado_prc_prep %>% 
  mutate(Depth = clado_counts$Depth, prc_scores = clado_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

clado_prcplot <- ggplot(clado_prc_plot_prep,
                        aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3, scales = "free_y") +
  theme_bw()

clado_prc_scores <- clado_counts %>% 
  select(Depth) %>% 
  mutate(value = clado_prc$lambda, param = rep("prc_score", nrow(clado_counts)))

#Rate of change-------
# 
# clado_source_community <- clado_counts %>% 
#   janitor::clean_names() %>% 
#   select(!depth) %>% 
#   mutate(sample_id = as.character(seq(11, length(apnap_counts$Depth) + 10, by = 1))) %>% 
#   relocate(sample_id)
# 
# clado_source_age <- clado_counts %>%
#   rename(depth = Depth) %>% 
#   left_join(tl18_age) %>% 
#   select(depth, ageBP) %>% 
#   rename(age = ageBP) %>% 
#   mutate(sample_id = clado_source_community$sample_id) %>% 
#   relocate(sample_id)
# 
# set_randomisations <- 10000
# 
# set.seed(12)
# clado_roc <-
#   RRatepol::estimate_roc(
#     data_source_community = clado_source_community,
#     data_source_age = clado_source_age,
#     smooth_method = "age.w",
#     dissimilarity_coefficient = "chisq",
#     working_units = "MW",
#     bin_size = 30,
#     number_of_shifts = 5,
#     standardise = TRUE,
#     n_individuals = round(min(clado_sum$count_sum)),
#     rand = set_randomisations,
#     use_parallel = FALSE
#   )
# 
# clado_roc_peaks <-
#   RRatepol::detect_peak_points(
#     data_source = clado_roc,
#     sel_method = "trend_non_linear"
#   ) %>%
#   mutate(Peak = as.character(Peak))
# 
# write_csv(clado_roc_peaks, "data/clado_roc_peaks.csv")

clado_roc_peaks <- read_csv("data/clado_roc_peaks.csv") %>% 
  mutate(ageAD = 1950-Age,
         Depth = approx(tl18_age_raw$ageAD, tl18_age_raw$depth, xout = ageAD)$y)

clado_roc_2join <- clado_roc_peaks %>% 
  mutate(value = ROC,
         param = "RoC") %>% 
  select(Depth, param, value)

clado_roc_peaks_2plot <- clado_roc_peaks %>% 
  select(Depth, Peak) %>% 
  filter(Peak == TRUE)

#combined PrC, RoC, concentration plots----------

clado_adds_plot_prep <- clado_conc %>% 
  add_row(clado_prc_scores) %>%
  filter(!param == "count_sum") %>% 
  add_row(clado_roc_2join) %>% 
  mutate(param = gsub("prc_score", "PrC score", param))

clado_adds_plot <- ggplot(clado_adds_plot_prep, aes(x = value, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_point(data = filter(clado_adds_plot_prep,
                           Depth %in%  clado_roc_peaks_2plot$Depth), color = "green", size = 3) +
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
  facet_geochem_gridh(vars(param), rotate_axis_labels = 90,
                                          units = c("Concentration" = "# cm⁻³", 
                                                    "PrC score" = "21%",
                                                    "RoC" = NA,
                                                    "CONISS" = "Total sum \n of squares")) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_dendrogram(clado_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(clado_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black") +
  rotated_facet_labels(
    angle = 90,
    direction = "x",
    remove_label_background = TRUE
  )

#wrap plot----------

clado_wrapped_plots <- wrap_plots(
  clado_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank())+
    labs(title = "(A)"),
  clado_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  clado_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 2)
)

#PCA ----
clado_pca_prep <- clado_perc %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(id_cols = Depth,
              names_from = taxon,
              values_from = rel_abund) %>% 
  select(!Depth)

clado_pca <- rda(sqrt(clado_pca_prep)) 
screeplot(clado_pca, bstick = TRUE)

clado_fort <- fortify(clado_pca, axes = c(1,2), scaling = "sites")

clado_sites <- clado_fort[clado_fort$Score %in% "sites",] %>% 
  mutate(depth = clado_counts$Depth)
clado_sp <- clado_fort[clado_fort$Score %in% "species",]

clado_inertcomp <- inertcomp(clado_pca) #contribution of each taxon to the total inertia
clado_sel_sp <- tibble(taxon = rownames(clado_inertcomp), inertcomp = clado_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
clado_sp_red <- clado_sp[which(clado_sp$Label %in% clado_sel_sp$taxon), ] #leave 10 taxa with the highest contrubution to the total inertia

clado_ve_prep <- clado_pca$CA$eig / clado_pca$tot.chi * 100
(clado_PC1_ve <- round(((clado_ve_prep / sum(clado_ve_prep))[c(1)]) * 100, digits = 1))#25.9% expl. var.
(clado_PC2_ve <- round(((clado_ve_prep / sum(clado_ve_prep))[c(2)]) * 100, digits = 1))#11.1% expl. var.

clado_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", clado_PC2_ve, "%)", sep = ""), x = paste("PC1 (", clado_PC1_ve, "%)", sep = ""), title = "(B)") +
  geom_segment(data = clado_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  geom_point(data = clado_sites, aes(x = PC1, y = PC2, color = depth)) +
  scale_color_viridis_c() +
  geom_point(data = filter(clado_sites, depth %in% c(5.75, 23.75, 31.75, 35.75)), aes(x = PC1, y = PC2), color = "red") +
  ggrepel::geom_text_repel(data = clado_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(depth))) +
  ggrepel::geom_text_repel(data = clado_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = clado_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', linewidth = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', linewidth = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2)) +
  theme(legend.position = "right")

#final plot --------

clado_final_plot <- wrap_plots(
  design = "AAAA
            #BB#",
  clado_wrapped_plots,
  clado_pca_plot,
  heights = c(5,7)
)

ggsave(filename="figures/clado_final.jpeg", 
       plot = clado_final_plot, 
       device = jpeg, 
       width = 12.5, 
       height = 12, 
       units = "in")

ggsave(filename="figures/clado_final.pdf", 
       plot = clado_final_plot, 
       device = pdf, 
       width = 12.5, 
       height = 12, 
       units = "in")

ggsave(filename="figures/clado_final.svg", 
       plot = clado_final_plot, 
       device = svg, 
       width = 12.5, 
       height = 12, 
       units = "in")
