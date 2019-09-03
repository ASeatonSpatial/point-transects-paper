# create polygon defining study area
# use simple convex hull and combindation with open forest polygon for now

library(sp)
library(rgeos)
library(inlabru)

samplers <- readRDS("samplers_extended_no_crs.RDS")
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

bnd.csv = read.csv("bnd_forest2.csv")
study_area_new <- SpatialPolygons(list(Polygons(list(Polygon(bnd.csv)), ID=1)))

ggplot() +
  gg(study_area_new) +
  gg(samplers) +
  gg(samplersnew, colour = "red") +
  coord_equal()

saveRDS(study_area_new, file = "study_area_extended_no_crs.RDS")
