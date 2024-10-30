# create polygon object defining study area

library(sf)
library(ggplot2)

data_raw_path = here::here("R", "data_raw")
data_path = here::here("R", "data")
samplers <- readRDS(here::here(data_path, "samplers.RDS"))

bnd.csv = read.csv(here::here(data_raw_path, "OpenClosedForestBoundary.csv"))
bnd.xy = data.frame(x = bnd.csv$EASTING / 1000,
                    y = bnd.csv$NORTHING / 1000)

# convert to sf polygon
bnd.poly <- st_polygon(list(as.matrix(bnd.xy)))
bnd.sfc <- st_sfc(bnd.poly)
bnd.sf <- st_sf(geometry = bnd.sfc)

# plot to check
ggplot() +
  geom_sf(data = bnd.sf, alpha = 0.2) +
  geom_sf(data = samplers)

saveRDS(bnd.sf,
        file = here::here(data_path,
                          "study_area.RDS"))
