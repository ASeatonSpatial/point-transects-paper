# create polygon defining study area
# use simple convex hull and combindation with open forest polygon for now

library(sp)
library(rgeos)
library(inlabru)

data_path = here::here("analyses", "data")
samplers <- readRDS(here::here(data_path, "samplers_extended_no_crs.RDS"))
samplersnew = samplers[samplers$Stratum == "CF",]
# old_study_area = readRDS("study_area_no_crs.RDS")
# chull = gConvexHull(samplersnew)
# chull = gBuffer(chull, width = 500, capStyle = "SQUARE")
#
# ggplot() +
#   gg(old_study_area) +
#   gg(chull) +
#   gg(samplers) +
#   gg(samplersnew, colour = "red") +
#   coord_equal()
#
# study_area_new = gUnion(old_study_area, chull)

# bnd.csv = read.csv("bnd_forest2.csv")
bnd.csv = read.csv(here::here(data_path, "OpenClosedForestBoundary.csv"))
bnd.xy = data.frame(x = bnd.csv$EASTING / 1000,
                    y = bnd.csv$NORTHING / 1000)
study_area_new <- SpatialPolygons(list(Polygons(list(Polygon(bnd.xy)), ID=1)))

ggplot() +
  gg(study_area_new) +
  gg(samplers) +
  gg(samplersnew, colour = "red") +
  coord_equal()

saveRDS(study_area_new, file = here::here(data_path, "study_area_extended_no_crs.RDS"))
