# Evaluate fitted model
library(INLA)
library(inlabru)
library(sp)

# set ggplot theme
theme_set(theme_minimal())

set.seed(9701071)

#### Fitted model  ####
file_name = "akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

study_area = readRDS("data/study_area_extended_no_crs.RDS")
samplers = readRDS("data/samplers_extended_no_crs.RDS")
mesh = readRDS("data/mesh_extended_no_crs.RDS")

# intensity

# Rick's grid already clipped to study region
# pred.grid <- expand.grid(EASTING=seq(from=255400, to=261200, by=200),
#                          NORTHING=seq(from=2189000, to=2200800, by=200))
# coordinates(pred.grid) = ~ EASTING + NORTHING
# keep = over(pred.grid, study_area, byid = TRUE)
# keep[is.na(keep)] = 0
# keep = as.logical(keep)

# using mgcv inSide() to match Rick's code
# bnd_points = study_area@polygons[[1]]@Polygons[[1]]@coords
# bnd_points = data.frame(EASTING = bnd_points[,1], NORTHING = bnd_points[,2])
# EASTING = pred.grid@coords[,1]
# NORTHING = pred.grid@coords[,2]
# keep = mgcv::inSide(bnd_points, EASTING, NORTHING)
# 
# pred.grid.studyarea = pred.grid[keep,]

# Rick's grid already clipped to study region
pred.grid.studyarea = read.csv("data/rick-pred-grid.csv")
coordinates(pred.grid.studyarea) = ~ EASTING + NORTHING

# to get back as matrix for dsm do 
# pred.grid.matrix = pred.grid.studyarea@coords

ggplot() + gg(pred.grid.studyarea) + gg(study_area) + coord_equal()

pred.pxl = SpatialPixels(pred.grid.studyarea)
# pr.int <- predict(fit, pred.pxl, ~ exp(grf + Intercept))
pr.int <- generate(fit, pred.pxl, ~ exp(grf + Intercept), n.samples = 1000)
pr.int = do.call(cbind, pr.int)

# p1 = ggplot() + 
#   gg(study_area) +
#   gg(pr.int) +
#   scale_fill_viridis_c() +
#   ggtitle("Posterior mean intensity") +  
#   xlab("Easting") +
#   ylab("Northing")
# 
# p2 = ggplot() +
#   gg(study_area)  +
#   gg(pr.int["sd"]) +
#   gg(samplers, colour = "red") +
#   scale_fill_viridis_c() +
#   ggtitle("Posterior SD of intensity") +  
#   xlab("Easting") +
#   ylab("Northing")
# 
# multiplot(p1, p2, cols = 2)

# pr.int.df = cbind(pr.int@data, pr.int@coords)
write.csv(pr.int, file = "gam_compare_intensity.csv")

# # log intensity
# pred.pxl = SpatialPixels(pred.grid.studyarea)
# pr.int <- predict(fit, pred.pxl, ~ grf + Intercept)
# 
# p1 = ggplot() + 
#   gg(study_area) +
#   gg(pr.int) +
#   scale_fill_viridis_c() +
#   ggtitle("Posterior mean intensity") +  
#   xlab("Easting") +
#   ylab("Northing")
# 
# p2 = ggplot() +
#   gg(study_area)  +
#   gg(pr.int["cv"]) +
#   gg(samplers, colour = "red") +
#   scale_fill_viridis_c() +
#   ggtitle("Posterior CV of intensity") +  
#   xlab("Easting") +
#   ylab("Northing")
# 
# multiplot(p1, p2, cols = 2)
# 
# pr.int.df = cbind(pr.int@data, pr.int@coords)
# write.csv(pr.int.df, file = "gam_compare_log_intensity.csv")
# 
# 
# # convert intensities to abundance - cell size is 200x200
# pr.ab <- predict(fit, pred.pxl, ~ 200*200*exp(grf + Intercept))
# 
# p1 = ggplot() + 
#   gg(study_area) +
#   gg(pr.ab) +
#   scale_fill_viridis_c() +
#   ggtitle("Posterior mean intensity") +  
#   xlab("Easting") +
#   ylab("Northing")
# 
# p2 = ggplot() +
#   gg(study_area)  +
#   gg(pr.ab["cv"]) +
#   gg(samplers, colour = "red") +
#   scale_fill_viridis_c() +
#   ggtitle("Posterior CV of intensity") +  
#   xlab("Easting") +
#   ylab("Northing")
# 
# multiplot(p1, p2, cols = 2)
# 
# pr.ab.df = cbind(pr.ab@data, pr.ab@coords)
# write.csv(pr.ab.df, file = "gam_compare_abundance.csv")


