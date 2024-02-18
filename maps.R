library(ggplot2)
library(spData)
library(maps)
library(ggspatial)
library(rnaturalearth)
library(ggmap)
library(patchwork)
library(geodata)
library(grid)
library(sf)
library(rcartocolor)
library(cowplot)
library(metR)
library(OpenStreetMap)
library(ggsn)
data(world)

nf <- gadm("GADM", country="CAN", level=0, resolution = 2)

nf2plot <- st_as_sf(nf)

nf_map <- ggplot(nf2plot) +
  geom_sf() +
  coord_sf(xlim = c(-60, -52),
           ylim = c(46, 52),
           expand = TRUE) +
  theme_bw() +
#  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering,
                         location = "tl") +
  labs(x = NULL, y = NULL, title = "(A)") +
  annotate("text", x = -56, y = 48.5, label = "Newfoundland") +
  annotate("text", x = -58.4, y = 46.8, label = "Burin Peninsula") +
  annotate("rect", xmin = -55.86, xmax = -55.83, ymin = 46.865, ymax = 46.89, color = "magenta", linewidth = 2) +
  annotate("text", x = -55.84, y = 46.68, label = "(B)") +
  annotate("segment", x = -57, xend = -55.7, y = 46.8, yend = 47, linewidth = 1,
           arrow = arrow(length = unit(0.2, "cm"))) 
  

inset_map <- ggplot(world) +
  geom_sf() +
  labs(x = NULL, y = NULL) +
  coord_sf(xlim = c(-120, -20),
           ylim = c(20, 60),
           expand = TRUE) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black")) +
  geom_rect(aes(xmin = -60 - 1, xmax = -52 + 1, ymin = 46 - 1, ymax = 52 + 1), color = "red", fill = NA, linewidth = 0.5)




bp2map_prep <- openmap(c(46.89,-55.86), c(46.865, -55.83),
                       type='bing',
                       mergeTiles = TRUE, 
                       zoom = 15)

bp2map <- openproj(bp2map_prep, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

ewbrks <- seq(-55.85, -55.83, by = 0.01)
nsbrks <- seq(46.87, 46.89, by = 0.01)


bp_map <- autoplot.OpenStreetMap(bp2map) +
  scalebar(x.min = -55.859, x.max = -55.831,
           y.min = 46.866, y.max = 46.89,
           dist = 0.4, dist_unit = "km",
           st.dist = 0.05, st.bottom = FALSE,
           st.color = "white", st.size = 3,
           location = "bottomleft", box.color = "white",
           border.size = 0.1,
           transform = TRUE, model = "WGS84") +
  annotation_north_arrow(style = north_arrow_fancy_orienteering,
                         location = "tr") +
  scale_y_latitude(breaks = nsbrks, position = "right") +
  scale_x_longitude(breaks = ewbrks) +
  labs(x = NULL, y = NULL, title = "(B)") +
  theme_bw() +
  annotate("point", x = -55.84346, y = 46.87533, color = "yellow") +
  annotate("text", x = -55.84346, y = 46.8745, label = "TL18-2", color = "yellow") +
  annotate("segment", x = -55.848, xend = -55.845, y = 46.883, yend = 46.876, color = "lightblue",
           linewidth = 1,
           arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x = -55.848, y = 46.884, label = "Broad Pond", color = "lightblue")

maps_wrapped <- wrap_plots(
  nf_map,
  bp_map)

final_map <- ggdraw(maps_wrapped) +
  draw_plot(inset_map, x = 0.225, y = 0.77, width = 0.25, height = 0.15)

ggsave(filename="figures/final_map.pdf",
       plot = final_map,
       device = pdf,
       width = 11,
       height = 5,
       units = "in")

ggsave(filename="figures/final_map.jpeg",
       plot = final_map,
       device = jpeg,
       width = 11,
       height = 5,
       units = "in")
