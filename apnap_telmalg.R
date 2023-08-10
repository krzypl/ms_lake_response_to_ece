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

#terrestrial pollen
apnap_counts <- read_csv("apnap_counts.csv")

pollen_concentration_reference <- read_csv("pollen_concentration_reference.csv")

conc_coef <- pollen_concentration_reference %>% 
  mutate(coef = lycopodium_ref/indicator)

tl18_sand <- read_csv("TL18_sand.csv")

sand <- tl18_sand %>% 
  filter(Depth < 38.5)

apnap_long <- apnap_counts %>% 
  pivot_longer(!Depth & !Age, names_to = "taxon", values_to = "count")

apnap_sum <- apnap_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count))

apnap_conc <- apnap_sum %>% 
  left_join(conc_coef) %>% 
  mutate(Concentration = coef*count_sum) %>% 
  pivot_longer(!Depth & !Age, names_to = "param", values_to = "value") %>% 
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

age_depth <- tibble(Depth = apnap_counts$Depth, Age = apnap_counts$Age)

tl18_adm <- age_depth_model(
  age_depth,
  depth = Depth,
  age = Age
)

apnap_prc_prep <- apnap_red %>%
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth) 

set.seed(12)
apnap_prc <- prcurve(
  apnap_prc_prep,
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

apnap_prc# variation explained by prc = 53%

apnap_prc_plot_prep <- apnap_prc_prep %>% 
  mutate(Depth = apnap_counts$Depth, prc_scores = apnap_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

apnap_prcplot <- ggplot(apnap_prc_plot_prep,
                        aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3) +
  theme_bw()

apnap_prc_scores <- apnap_counts %>% 
  select(Depth) %>% 
  mutate(value = apnap_prc$lambda, param = rep("prc_score", nrow(apnap_counts)))

apnap_sand_plot <- ggplot(sand, aes(x = count, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

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
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

telmalg_counts <- read_csv("telmalg_counts.csv")

telmalg_long <- telmalg_counts %>% 
  pivot_longer(!Depth & !Age, names_to = "taxon", values_to = "count") %>% 
  mutate(taxon = as.factor(taxon))

telmalg_sum <- telmalg_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count))

telmalg_conc <- telmalg_sum %>% 
  left_join(conc_coef) %>% 
  mutate(Concentration = coef*count_sum) %>% 
  pivot_longer(!Depth & !Age, names_to = "param", values_to = "value") %>% 
  filter(param == "Concentration")

telmalg_perc <- telmalg_long %>%
  left_join(apnap_sum) %>% 
  rename(count_sum_apnap = count_sum) %>% 
  left_join(telmalg_sum) %>% 
  mutate(rel_abund = count/(count_sum_apnap + count_sum) * 100)

telmalg_zero <- telmalg_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(telmalg_counts$Depth) - 1) #select taxa with at least 2 occurences 

telmalg_red2 <- telmalg_perc %>% 
  filter(!taxon %in% telmalg_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>%
  ungroup()

telmalg_red1 <- telmalg_perc %>% 
  filter(!taxon %in% telmalg_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 1) %>%
  ungroup()

Sphagnum_plot <- telmalg_red1 %>% 
  filter(taxon == "Sphagnum") %>% 
  ggplot(aes(x = rel_abund, y = Depth)) +
  geom_col_segsh(color = "darkblue") + 
  geom_lineh(color = "black") +
  facet_abundanceh(vars(taxon)) +
  labs(x = "Relative abundance (%)", y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

apnap_adds_plot <- apnap_conc %>% 
  add_row(apnap_prc_scores) %>% 
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
  layer_dendrogram(apnap_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(apnap_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

wrap_plots(
  apnap_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  apnap_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  Sphagnum_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(x = NULL, y = NULL),
  apnap_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 10*(1/16), 3)
)

#telmatophytes and algae (dane wczytane wczesniej, zeby dodac sphagnum do plota)
telmalg_coniss <- telmalg_red2 %>%
  filter(!taxon == "Sphagnum") %>% 
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

# terrestrial pollen ------------------------------------------------------


telmalg_prc_prep <- telmalg_red2 %>%
  filter(!taxon == "Sphagnum") %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth)

set.seed(12)
telmalg_prc <- prcurve(
  telmalg_prc_prep,
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

telmalg_prc# variation explained by prc = 84%

telmalg_prc_plot_prep <- telmalg_prc_prep %>% 
  mutate(Depth = telmalg_counts$Depth, prc_scores = telmalg_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

telmalg_prcplot <- ggplot(telmalg_prc_plot_prep,
                          aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3) +
  theme_bw()

telmalg_prc_scores <- telmalg_counts %>% 
  select(Depth) %>% 
  mutate(value = telmalg_prc$lambda, param = rep("prc_score", nrow(telmalg_counts)))

telmalg_sand_plot <- ggplot(sand, aes(x = count, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(telmalg_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

telmalg_plot <- telmalg_red1 %>% 
  filter(!taxon == "Sphagnum") %>%
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
  layer_zone_boundaries(telmalg_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

telmalg_adds_plot <- telmalg_conc %>% 
  add_row(telmalg_prc_scores) %>% 
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
  layer_dendrogram(telmalg_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(telmalg_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

#przy wrap plocie jest inna skala dla scenedesmusa i pediastrum!!!
wrap_plots(
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

#diatoms
diatom_counts <- read_csv("https://raw.githubusercontent.com/krzypl/ms_refininig_history_of_ece/main/data/tl18_diatom.csv")

diatom_long <- diatom_counts %>% 
  pivot_longer(!Depth & !Age, names_to = "taxon", values_to = "count")

diatom_sum <- diatom_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count))

diatom_m <- diatom_counts %>% 
  select(Depth, "Nitzschia tubicola" : "Cocconeis scutelum") %>% 
  pivot_longer(!Depth, names_to = "taxa", values_to = "count") %>% 
  group_by(Depth) %>% 
  summarise(rel_abund = sum(count)/400*100) %>% 
  ungroup() %>% 
  mutate(eco_group = rep("Marine", length(diatom_counts$Depth)))

diatom_b <- diatom_counts %>% 
  select(Depth, "Parlibellus plicatus" : "Navicula veneta") %>% 
  pivot_longer(!Depth, names_to = "taxa", values_to = "count") %>% 
  group_by(Depth) %>% 
  summarise(rel_abund = sum(count)/400*100) %>% 
  ungroup() %>% 
  mutate(eco_group = rep("Brackish", length(diatom_counts$Depth)))

diatom_e <- diatom_counts %>% 
  select(Depth, "Pinnuavis elegans" : "Nitzschia dubia") %>% 
  pivot_longer(!Depth, names_to = "taxa", values_to = "count") %>% 
  group_by(Depth) %>% 
  summarise(rel_abund = sum(count)/400*100) %>% 
  ungroup() %>% 
  mutate(eco_group = rep("Euryhaline", length(diatom_counts$Depth)))

diatom_f <- diatom_counts %>% 
  select(Depth, "Ctenophora pulchella" : "Stauroneis acidoclinata") %>% 
  pivot_longer(!Depth, names_to = "taxa", values_to = "count") %>% 
  group_by(Depth) %>% 
  summarise(rel_abund = sum(count)/400*100) %>% 
  ungroup() %>% 
  mutate(eco_group = rep("Freshwater", length(diatom_counts$Depth)))

diatom_i <- diatom_counts %>% 
  select(Depth, "Stauroneis capitata" : "Pinnularia sp.") %>% 
  pivot_longer(!Depth, names_to = "taxa", values_to = "count") %>% 
  group_by(Depth) %>% 
  summarise(rel_abund = sum(count)/400*100) %>% 
  ungroup() %>% 
  mutate(eco_group = rep("Indifferent", length(diatom_counts$Depth)))

diatom_eco_groups <- diatom_m %>% 
  add_row(diatom_b) %>% 
  add_row(diatom_e) %>% 
  add_row(diatom_f) %>% 
  add_row(diatom_i)

diatom_perc <- diatom_long %>%
  left_join(diatom_sum) %>% 
  mutate(rel_abund = count/count_sum*100)

diatom_zero <- diatom_perc %>% 
  group_by(taxon) %>%
  filter(rel_abund == 0) %>% 
  summarise(nb_zero = n()) %>% 
  filter(nb_zero >= length(diatom_counts$Depth) - 2) #select taxa with at least 2 occurences 

diatom_red <- diatom_perc %>% 
  filter(!taxon %in% diatom_zero$taxon) %>% 
  group_by(taxon) %>% 
  filter(max(rel_abund) > 2) %>% 
  ungroup()

diatom_coniss <- diatom_red %>%
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

diatom_prc_prep <- diatom_red %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth) 

set.seed(12)
diatom_prc <- prcurve(
  diatom_prc_prep,
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

diatom_prc# variation explained by prc = 52%

diatom_prc_plot_prep <- diatom_prc_prep %>% 
  mutate(Depth = diatom_counts$Depth, prc_scores = diatom_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

diatom_prcplot <- ggplot(diatom_prc_plot_prep,
                         aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3) +
  theme_bw()

diatom_prc_scores <- diatom_counts %>% 
  select(Depth) %>% 
  mutate(value = diatom_prc$lambda, param = rep("prc_score", nrow(diatom_counts)))

diatom_sand_plot <- ggplot(sand, aes(x = count, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(diatom_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

diatom_order <- unique(diatom_red$taxon)

diatom_plot <- diatom_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, diatom_order)) %>% 
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
  layer_zone_boundaries(diatom_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

diatom_eco_groups_order <- unique(diatom_eco_groups$eco_group)

diatom_eco_group_plot <- diatom_eco_groups %>% 
  mutate(eco_group = as.factor(eco_group)) %>% 
  mutate(eco_group = fct_relevel(eco_group, diatom_eco_groups_order)) %>% 
  ggplot(aes(x = rel_abund, y = Depth)) +
  geom_col_segsh() + 
  geom_lineh() +
  facet_abundanceh(vars(eco_group)) +
  labs(x = NULL, y = "Depth (cm)") +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  layer_zone_boundaries(diatom_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

diatom_adds_plot <- ggplot(diatom_prc_scores, aes(x = value, y = Depth)) +
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
  layer_dendrogram(diatom_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(diatom_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

wraped_diatoms <- wrap_plots(
  diatom_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  diatom_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  diatom_eco_group_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  diatom_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 2, 2)
)

#ggsave(filename = "wraped_diatoms.svg",
#       plot = print(wraped_diatoms),
#       device = CairoSVG)

#chironomidae
chiro_counts <- read_csv("chiro_counts.csv")

chiro_long <- chiro_counts %>% 
  pivot_longer(!Depth & !Age & !volume, names_to = "taxon", values_to = "count")

chiro_sum <- chiro_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count)) %>%
  ungroup() %>% 
  mutate(volume = chiro_counts$volume)

chiro_conc <- chiro_sum %>% 
  mutate(Concentration = count_sum/volume, Age = chiro_counts$Age) %>% 
  pivot_longer(!Depth & !Age, names_to = "param", values_to = "value") %>% 
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

#do prc usuwam wszystkie obserwacje onizej najnizeszego piku piasku, bo zliaczenia sa tam max 3, wiec nie ma sensu ich uwzgldniach w obliczeniach
chiro_prc_prep <- chiro_red %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>%
  filter(Depth <= 35.75) %>% 
  select(!Depth) 

set.seed(12)
chiro_prc <- prcurve(
  chiro_prc_prep,
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

chiro_prc# variation explained by prc = 20%

chiro_bot_na <- chiro_red %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>% 
  filter(Depth > 35.75) %>% 
  select(!Depth) %>% 
  mutate(across(c(1:36), 
                ~ ifelse(Ablabesmyia == 1, 1, NA)))

chiro_prc_plot_prep <- chiro_prc_prep %>% 
  add_row(chiro_bot_na) %>% 
  mutate(Depth = chiro_counts$Depth, prc_scores = c(chiro_prc$lambda, NA, NA, NA, NA, NA)) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

chiro_prcplot <- ggplot(chiro_prc_plot_prep,
                        aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3) +
  theme_bw()

chiro_prc_scores <- chiro_counts %>% 
  select(Depth) %>% 
  mutate(value = c(chiro_prc$lambda, NA, NA, NA, NA, NA), param = rep("prc_score", nrow(chiro_counts)))

chiro_sand_plot <- ggplot(sand, aes(x = count, y = Depth)) +
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

chiro_adds_plot <- chiro_conc %>% 
  add_row(chiro_prc_scores) %>%
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

wrap_plots(
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

#cladocera
clado_counts <- read_csv("clado_counts.csv")

clado_long <- clado_counts %>% 
  pivot_longer(!Depth & !Age & !volume, names_to = "taxon", values_to = "count")

clado_sum <- clado_long %>% 
  group_by(Depth) %>% 
  summarise(count_sum = sum(count)) %>%
  ungroup() %>% 
  mutate(volume = clado_counts$volume)

clado_conc <- clado_sum %>% 
  mutate(Concentration = count_sum/volume, Age = clado_counts$Age) %>% 
  pivot_longer(!Depth & !Age, names_to = "param", values_to = "value") %>% 
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

clado_coniss <- clado_red %>%
  mutate(rel_abund_trans = rel_abund/100) %>% 
  nested_data(qualifiers = Depth, key = taxon, value = rel_abund_trans, trans = sqrt) %>% 
  nested_chclust_coniss()

clado_prc_prep <- clado_red %>% 
  select(Depth, taxon, rel_abund) %>% 
  pivot_wider(Depth, names_from = taxon, values_from = rel_abund) %>% 
  select(!Depth) 

set.seed(12)
clado_prc <- prcurve(
  clado_prc_prep,
  method = "ca",
  smoother = smoothSpline,
  trace = TRUE,
  vary = FALSE,
  penalty = 1.4
) 

clado_prc# variation explained by prc = 57%

clado_prc_plot_prep <- clado_prc_prep %>% 
  mutate(Depth = clado_counts$Depth, prc_scores = clado_prc$lambda) %>% 
  pivot_longer(!Depth & !prc_scores, names_to = "taxon", values_to = "rel_abund")

clado_prcplot <- ggplot(clado_prc_plot_prep,
                        aes(x = prc_scores, y = rel_abund)) +
  geom_point() +
  geom_smooth(span = 0.3, se = FALSE) +
  facet_wrap(~ taxon, nrow = 3) +
  theme_bw()

clado_prc_scores <- clado_counts %>% 
  select(Depth) %>% 
  mutate(value = clado_prc$lambda, param = rep("prc_score", nrow(clado_counts)))

clado_sand_plot <- ggplot(sand, aes(x = count, y = Depth)) +
  geom_lineh() +
  geom_point() +
  geom_hline(yintercept = c(5.75, 23.75, 31.75, 35.75),
             col = "darkblue", lty = 1, alpha = 0.1, size = 2) +
  scale_y_reverse(breaks = c(breaks = seq(0, 40, by = 5)), 
                  labels = as.character(c(breaks = seq(0, 40, by = 5))),
                  expand = expansion(mult = c(0.02, 0.02))) +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  layer_zone_boundaries(clado_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

clado_order <- unique(clado_red$taxon)

clado_plot <- clado_red %>% 
  mutate(taxon = as.factor(taxon)) %>% 
  mutate(taxon = fct_relevel(taxon, clado_order)) %>% 
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
  layer_zone_boundaries(clado_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

clado_adds_plot <- clado_conc %>% 
  add_row(clado_prc_scores) %>%
  filter(!param == "count_sum") %>% 
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
  layer_dendrogram(clado_coniss, aes(y = Depth), param = "CONISS") +
  layer_zone_boundaries(clado_coniss, aes(y = Depth), linetype = 2, size = 0.5, colour = "black")

wrap_plots(
  clado_sand_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  clado_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  clado_adds_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 10, 2)
)
#########rda
########probowalem z wartoscia = 100 dla wsystkich probek, ktore przekroczyly treshold w peak analysis. wynik gorszy
expl_var <- as_tibble(read.csv("expl_variables.csv", sep = ";"))#dodac tutaj jeszcze apnap_prc, ktore oryginalnie bylo juz w pliku. ja usunalem, zeby bylo zassane bezposrednio z environment.

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

emi_effect_prep <- sand_peaks %>% 
  add_row(time_elapsed_ev1) %>% 
  add_row(time_elapsed_ev2) %>% 
  add_row(time_elapsed_ev3) %>%
  add_row(time_elapsed_ev4) %>%
  add_row(time_elapsed_ev5)

expl_var_full <- expl_var %>% 
  left_join(emi_effect_prep) %>% 
  mutate(emi_effect = 100*exp(1)^(-0.3*time_elapsed_since_event))

ggplot(expl_var_full, aes(x = Age, y = emi_effect)) + 
  geom_point() +
  geom_line()

#diatom rda
diatom_wide <- diatom_red %>% 
  pivot_wider(!count & !count_sum, names_from = "taxon", values_from = "rel_abund") %>% 
  select(!Depth & !Age)

diatom_rda <- rda(sqrt(diatom_wide) ~ emi_effect + Condition(Age, apnap_prc), data = expl_var_full)

set.seed(12)
(diatom_R2adj <- RsquareAdj(diatom_rda)$adj.r.squared)#10.0% 

set.seed(12)
(diatom_anova <- anova(diatom_rda, permutations = how(nperm = 999)))#significance level 0.01**

diatom_rda_all <- rda(sqrt(diatom_wide) ~ emi_effect + Age + apnap_prc, data = expl_var_full)

set.seed(12)
(diatom_R2adj_all <- RsquareAdj(diatom_rda_all)$adj.r.squared)#25.5% 

set.seed(12)
(diatom_anova_all <- anova(diatom_rda_all, permutations = how(nperm = 999)))#significance level 0.001***

diatom_rda_age <- rda(sqrt(diatom_wide) ~ Age + Condition(apnap_prc, emi_effect), data = expl_var_full)

set.seed(12)
(diatom_R2adj_age <- RsquareAdj(diatom_rda_age)$adj.r.squared)#6.3% 

set.seed(12)
(diatom_anova_age <- anova(diatom_rda_age, permutations = how(nperm = 999)))#significance level 0.01**

diatom_rda_apnap <- rda(sqrt(diatom_wide) ~ apnap_prc + Condition(Age, emi_effect), data = expl_var_full)

set.seed(12)
(diatom_R2adj_apnap <- RsquareAdj(diatom_rda_apnap)$adj.r.squared)#7.9% 

set.seed(12)
(diatom_anova_apnap <- anova(diatom_rda_apnap, permutations = how(nperm = 999)))#significance level 0.01**

diatom_rda_tprep <- tibble("EMI_xage_xapnap" = diatom_R2adj,
                           "apnap_xEMI_xage" = diatom_R2adj_age,
                           "age_xEMI_xapnap" = diatom_R2adj_apnap,
                           "EMI_age_apnap_shared" = 
                             diatom_R2adj_all - diatom_R2adj - diatom_R2adj_age - diatom_R2adj_apnap,
                           "unexplained" = 1 - diatom_R2adj_all,
                           "proxy" = "diatom")

#telmalg rda
telmalg_wide <- telmalg_red2 %>% 
  pivot_wider(!count & !count_sum, names_from = "taxon", values_from = "rel_abund") %>% 
  select(!c(Depth, Age, count_sum_apnap, Sphagnum))

telmalg_rda <- rda(sqrt(telmalg_wide) ~ emi_effect + Condition(Age, apnap_prc), data = expl_var_full)

set.seed(12)
(telmalg_R2adj <- RsquareAdj(telmalg_rda)$adj.r.squared)#-1.2% 

set.seed(12)
(telmalg_anova <- anova(telmalg_rda, permutations = how(nperm = 999)))#not significant

telmalg_rda_all <- rda(sqrt(telmalg_wide) ~ emi_effect + Age + apnap_prc, data = expl_var_full)

set.seed(12)
(telmalg_R2adj_all <- RsquareAdj(telmalg_rda_all)$adj.r.squared)#17.6% 

set.seed(12)
(telmalg_anova_all <- anova(telmalg_rda_all, permutations = how(nperm = 999)))#significance level 0.01**

telmalg_rda_age <- rda(sqrt(telmalg_wide) ~ Age + Condition(apnap_prc, emi_effect), data = expl_var_full)

set.seed(12)
(telmalg_R2adj_age <- RsquareAdj(telmalg_rda_age)$adj.r.squared)#7.7% 

set.seed(12)
(telmalg_anova_age <- anova(telmalg_rda_age, permutations = how(nperm = 999)))#significance level 0.05*

telmalg_rda_apnap <- rda(sqrt(telmalg_wide) ~ apnap_prc + Condition(Age, emi_effect), data = expl_var_full)

set.seed(12)
(telmalg_R2adj_apnap <- RsquareAdj(telmalg_rda_apnap)$adj.r.squared)#8.2% 

set.seed(12)
(telmalg_anova_apnap <- anova(telmalg_rda_apnap, permutations = how(nperm = 999)))#significance level 0.05*

telmalg_rda_tprep <- tibble("EMI_xage_xapnap" = telmalg_R2adj,
                            "apnap_xEMI_xage" = telmalg_R2adj_age,
                            "age_xEMI_xapnap" = telmalg_R2adj_apnap,
                            "EMI_age_apnap_shared" = 
                              telmalg_R2adj_all - telmalg_R2adj - telmalg_R2adj_age - telmalg_R2adj_apnap,
                            "unexplained" = 1 - telmalg_R2adj_all,
                            "proxy" = "telmalg")

#chiro rda
chiro_wide <- chiro_prc_prep

chiro_rda <- rda(sqrt(chiro_wide) ~ emi_effect + Condition(Age, apnap_prc), data = expl_var_full[c(1:37),])

set.seed(12)
(chiro_R2adj <- RsquareAdj(chiro_rda)$adj.r.squared)#-0.1% 

set.seed(12)
(chiro_anova <- anova(chiro_rda, permutations = how(nperm = 999)))#not significant

chiro_rda_all <- rda(sqrt(chiro_wide) ~ emi_effect + Age + apnap_prc, data = expl_var_full[c(1:37),])

set.seed(12)
(chiro_R2adj_all <- RsquareAdj(chiro_rda_all)$adj.r.squared)#7.4% 

set.seed(12)
(chiro_anova_all <- anova(chiro_rda_all, permutations = how(nperm = 999)))#significance level 0.001***

chiro_rda_age <- rda(sqrt(chiro_wide) ~ Age + Condition(apnap_prc, emi_effect), data = expl_var_full[c(1:37),])

set.seed(12)
(chiro_R2adj_age <- RsquareAdj(chiro_rda_age)$adj.r.squared)#3.9% 

set.seed(12)
(chiro_anova_age <- anova(chiro_rda_age, permutations = how(nperm = 999)))#significance level 0.001***

chiro_rda_apnap <- rda(sqrt(chiro_wide) ~ apnap_prc + Condition(Age, emi_effect), data = expl_var_full[c(1:37),])

set.seed(12)
(chiro_R2adj_apnap <- RsquareAdj(chiro_rda_apnap)$adj.r.squared)#4.5% 

set.seed(12)
(chiro_anova_apnap <- anova(chiro_rda_apnap, permutations = how(nperm = 999)))#significance level 0.001**

chiro_rda_tprep <- tibble("EMI_xage_xapnap" = chiro_R2adj,
                          "apnap_xEMI_xage" = chiro_R2adj_age,
                          "age_xEMI_xapnap" = chiro_R2adj_apnap,
                          "EMI_age_apnap_shared" = 
                            chiro_R2adj_all - chiro_R2adj - chiro_R2adj_age - chiro_R2adj_apnap,
                          "unexplained" = 1 - chiro_R2adj_all,
                          "proxy" = "chiro")

#clado rda
clado_wide <- clado_red %>% 
  pivot_wider(!count & !count_sum, names_from = "taxon", values_from = "rel_abund") %>% 
  select(!Depth & !Age)

clado_rda <- rda(sqrt(clado_wide) ~ emi_effect + Condition(Age, apnap_prc), data = expl_var_full)

set.seed(12)
(clado_R2adj <- RsquareAdj(clado_rda)$adj.r.squared)#-0.8% 

set.seed(12)
(clado_anova <- anova(clado_rda, permutations = how(nperm = 999)))#not significant

clado_rda_all <- rda(sqrt(clado_wide) ~ emi_effect + Age + apnap_prc, data = expl_var_full)

set.seed(12)
(clado_R2adj_all <- RsquareAdj(clado_rda_all)$adj.r.squared)#18.6% 

set.seed(12)
(clado_anova_all <- anova(clado_rda_all, permutations = how(nperm = 999)))#significance level 0.001***

clado_rda_age <- rda(sqrt(clado_wide) ~ Age + Condition(apnap_prc, emi_effect), data = expl_var_full)

set.seed(12)
(clado_R2adj_age <- RsquareAdj(clado_rda_age)$adj.r.squared)#9.4% 

set.seed(12)
(clado_anova_age <- anova(clado_rda_age, permutations = how(nperm = 999)))#significance level 0.001***

clado_rda_apnap <- rda(sqrt(clado_wide) ~ apnap_prc + Condition(Age, emi_effect), data = expl_var_full)

set.seed(12)
(clado_R2adj_apnap <- RsquareAdj(clado_rda_apnap)$adj.r.squared)#10.5% 

set.seed(12)
(clado_anova_apnap <- anova(clado_rda_apnap, permutations = how(nperm = 999)))#significance level 0.001***

clado_rda_tprep <- tibble("EMI_xage_xapnap" = clado_R2adj,
                          "apnap_xEMI_xage" = clado_R2adj_age,
                          "age_xEMI_xapnap" = clado_R2adj_apnap,
                          "EMI_age_apnap_shared" = 
                            clado_R2adj_all - clado_R2adj - clado_R2adj_age - clado_R2adj_apnap,
                          "unexplained" = 1 - clado_R2adj_all,
                          "proxy" = "clado")

#tabele RDA
rda_table <- diatom_rda_tprep %>% 
  add_row(telmalg_rda_tprep) %>% 
  add_row(chiro_rda_tprep) %>% 
  add_row(clado_rda_tprep) %>% 
  pivot_longer(!proxy, 
               names_to = "category", values_to = "proportion_of_v_expl")

ggplot(rda_table, aes(x = proxy, y = proportion_of_v_expl, fill = category)) +
  geom_col() +
  scale_fill_viridis_d(labels = c("Age", "Vegetation shifts", "Age, VS, and EMI shared", "EMI", "unexplained"))

#pca biplots: diatoms
diatom_pca <- rda(sqrt(diatom_wide))

diatom_fort <- fortify(diatom_pca, axes = c(1,2), scaling = "sites")

diatom_sites <- diatom_fort[diatom_fort$Score %in% "sites",]
diatom_sp <- diatom_fort[diatom_fort$Score %in% "species",]

diatom_inertcomp <- inertcomp(diatom_pca) #obliczam jaki ma udzial kazdy z taksonow dla calosci inertii
diatom_sel_sp <- tibble(taxon = rownames(diatom_inertcomp), inertcomp = diatom_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
diatom_sp_red <- diatom_sp[which(diatom_sp$Label %in% diatom_sel_sp$taxon), ] #pozostawiam 10 taksonow, ktory maja najwiekszy wklad w total inertia

diatom_ve_prep <- diatom_pca$CA$eig / diatom_pca$tot.chi * 100
(diatom_PC1_ve <- round(((diatom_ve_prep / sum(diatom_ve_prep))[c(1)]) * 100, digits = 1))#35.6% expl. var.
(diatom_PC2_ve <- round(((diatom_ve_prep / sum(diatom_ve_prep))[c(2)]) * 100, digits = 1))#13.3% expl. var.

pca_color <- c(
  rep(times = length(time_elapsed_ev1$Depth), "blue"),
  "red",
  rep(times = length(time_elapsed_ev2$Depth), "darkblue"),
  "red",
  rep(times = length(time_elapsed_ev3$Depth), "gray30"),
  "red",
  rep(times = length(time_elapsed_ev4$Depth), "gray60"),
  "red",
  rep(times = length(time_elapsed_ev5$Depth), "black")
)

diatom_sites <- diatom_sites %>% 
  mutate(pca_col = pca_color)

diatom_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", diatom_PC2_ve, "%)", sep = ""), x = paste("PC1 (", diatom_PC1_ve, "%)", sep = ""), title = "Diatom PCA plot") +
  geom_point(data = diatom_sites, aes(x = PC1, y = PC2, color = pca_col), color = pca_color) +
  ggrepel::geom_text_repel(data = diatom_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(diatom_counts$Depth))) +
  labs(color = "zone") +
  geom_segment(data = diatom_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(data = diatom_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = diatom_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', size = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', size = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))

#pca biplots: clado
clado_pca <- rda(sqrt(clado_wide))

clado_fort <- fortify(clado_pca, axes = c(1,2), scaling = "sites")

clado_sites <- clado_fort[clado_fort$Score %in% "sites",]
clado_sp <- clado_fort[clado_fort$Score %in% "species",]

clado_inertcomp <- inertcomp(clado_pca) #obliczam jaki ma udzial kazdy z taksonow dla calosci inertii
clado_sel_sp <- tibble(taxon = rownames(clado_inertcomp), inertcomp = clado_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
clado_sp_red <- clado_sp[which(clado_sp$Label %in% clado_sel_sp$taxon), ] #pozostawiam 10 taksonow, ktory maja najwiekszy wklad w total inertia

clado_ve_prep <- clado_pca$CA$eig / clado_pca$tot.chi * 100
(clado_PC1_ve <- round(((clado_ve_prep / sum(clado_ve_prep))[c(1)]) * 100, digits = 1))#21.2% expl. var.
(clado_PC2_ve <- round(((clado_ve_prep / sum(clado_ve_prep))[c(2)]) * 100, digits = 1))#16.3% expl. var.

pca_color <- c(
  rep(times = length(time_elapsed_ev1$Depth), "blue"),
  "red",
  rep(times = length(time_elapsed_ev2$Depth), "darkblue"),
  "red",
  rep(times = length(time_elapsed_ev3$Depth), "gray30"),
  "red",
  rep(times = length(time_elapsed_ev4$Depth), "gray60"),
  "red",
  rep(times = length(time_elapsed_ev5$Depth), "black")
)

clado_sites <- clado_sites %>% 
  mutate(pca_col = pca_color)

clado_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", clado_PC2_ve, "%)", sep = ""), x = paste("PC1 (", clado_PC1_ve, "%)", sep = ""), title = "clado PCA plot") +
  geom_point(data = clado_sites, aes(x = PC1, y = PC2, color = pca_col), color = pca_color) +
  ggrepel::geom_text_repel(data = clado_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(clado_counts$Depth))) +
  labs(color = "zone") +
  geom_segment(data = clado_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(data = clado_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = clado_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', size = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', size = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))

#pca biplots: chiro
chiro_pca <- rda(sqrt(chiro_wide))

chiro_fort <- fortify(chiro_pca, axes = c(1,2), scaling = "sites")

chiro_sites <- chiro_fort[chiro_fort$Score %in% "sites",]
chiro_sp <- chiro_fort[chiro_fort$Score %in% "species",]

chiro_inertcomp <- inertcomp(chiro_pca) #obliczam jaki ma udzial kazdy z taksonow dla calosci inertii
chiro_sel_sp <- tibble(taxon = rownames(chiro_inertcomp), inertcomp = chiro_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
chiro_sp_red <- chiro_sp[which(chiro_sp$Label %in% chiro_sel_sp$taxon), ] #pozostawiam 10 taksonow, ktory maja najwiekszy wklad w total inertia

chiro_ve_prep <- chiro_pca$CA$eig / chiro_pca$tot.chi * 100
(chiro_PC1_ve <- round(((chiro_ve_prep / sum(chiro_ve_prep))[c(1)]) * 100, digits = 1))#15.1% expl. var.
(chiro_PC2_ve <- round(((chiro_ve_prep / sum(chiro_ve_prep))[c(2)]) * 100, digits = 1))#12.8% expl. var.

pca_color_chiro <- c(
  rep(times = length(time_elapsed_ev1$Depth), "blue"),
  "red",
  rep(times = length(time_elapsed_ev2$Depth), "darkblue"),
  "red",
  rep(times = length(time_elapsed_ev3$Depth), "gray30"),
  "red",
  rep(times = length(time_elapsed_ev4$Depth), "gray60"),
  "red"
)

chiro_sites <- chiro_sites %>% 
  mutate(pca_col = pca_color_chiro)

chiro_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", chiro_PC2_ve, "%)", sep = ""), x = paste("PC1 (", chiro_PC1_ve, "%)", sep = ""), title = "chiro PCA plot") +
  geom_point(data = chiro_sites, aes(x = PC1, y = PC2, color = pca_col), color = pca_color_chiro) +
  ggrepel::geom_text_repel(data = chiro_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = chiro_counts$Depth[c(1:37)])) +
  labs(color = "zone") +
  geom_segment(data = chiro_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(data = chiro_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = chiro_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', size = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', size = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))

#pca biplots: telmalg
telmalg_pca <- rda(sqrt(telmalg_wide))

telmalg_fort <- fortify(telmalg_pca, axes = c(1,2), scaling = "sites")

telmalg_sites <- telmalg_fort[telmalg_fort$Score %in% "sites",]
telmalg_sp <- telmalg_fort[telmalg_fort$Score %in% "species",]

telmalg_inertcomp <- inertcomp(telmalg_pca) #obliczam jaki ma udzial kazdy z taksonow dla calosci inertii
telmalg_sel_sp <- tibble(taxon = rownames(telmalg_inertcomp), inertcomp = telmalg_inertcomp[,1]) %>%
  mutate(inert_ratio = inertcomp/sum(inertcomp), inert_rank = rank(desc(inert_ratio))) %>%
  filter(inert_rank <= 10)
telmalg_sp_red <- telmalg_sp[which(telmalg_sp$Label %in% telmalg_sel_sp$taxon), ] #pozostawiam 10 taksonow, ktory maja najwiekszy wklad w total inertia

telmalg_ve_prep <- telmalg_pca$CA$eig / telmalg_pca$tot.chi * 100
(telmalg_PC1_ve <- round(((telmalg_ve_prep / sum(telmalg_ve_prep))[c(1)]) * 100, digits = 1))#66.7% expl. var.
(telmalg_PC2_ve <- round(((telmalg_ve_prep / sum(telmalg_ve_prep))[c(2)]) * 100, digits = 1))#20.0% expl. var.

pca_color <- c(
  rep(times = length(time_elapsed_ev1$Depth), "blue"),
  "red",
  rep(times = length(time_elapsed_ev2$Depth), "darkblue"),
  "red",
  rep(times = length(time_elapsed_ev3$Depth), "gray30"),
  "red",
  rep(times = length(time_elapsed_ev4$Depth), "gray60"),
  "red",
  rep(times = length(time_elapsed_ev5$Depth), "black")
)

telmalg_sites <- telmalg_sites %>% 
  mutate(pca_col = pca_color)

telmalg_pca_plot <- ggplot() +
  labs(y = paste("PC2 (", telmalg_PC2_ve, "%)", sep = ""), x = paste("PC1 (", telmalg_PC1_ve, "%)", sep = ""), title = "telmalg PCA plot") +
  geom_point(data = telmalg_sites, aes(x = PC1, y = PC2, color = pca_col), color = pca_color) +
  ggrepel::geom_text_repel(data = telmalg_sites, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = as.character(telmalg_counts$Depth))) +
  labs(color = "zone") +
  geom_segment(data = telmalg_sp_red,
               color = "black", size = 0.7,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = 
                 grid::arrow(length = grid::unit(0.25, "cm"))) +
  ggrepel::geom_text_repel(data = telmalg_sp_red, color = "black",
                           size = 2.5, segment.alpha = 0,
                           aes(x = PC1, y = PC2, 
                               label = telmalg_sp_red$Label)) +
  geom_vline(xintercept = 0, color = 'black', size = 0.6,linetype=2) + 
  geom_hline(yintercept = 0, color = 'black', size = 0.6,linetype=2) +
  theme(legend.position = "bottom", panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_line(color = "grey80", size = 0.2))
