# #install.packages("weathercan", 
#                  repos = c("https://ropensci.r-universe.dev", 
#                            "https://cloud.r-project.org"))
library(weathercan)
library(tidyverse)

lamaline_closests <- stations_search(coords = c(46.8680556, -55.805277), dist = 100, interval = "month") %>% 
  filter(normals_1981_2010 == TRUE)

w_data <- normals_dl(lamaline_closests$climate_id[[1]])$normals[[1]]