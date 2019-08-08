# Point Transect Distance Sampling With Inlabru
## A Case Study with Hawaiian Akepa Distance Sampling Data

A brief summary of contents of this repo:


### Data files

obs_2002_no_crs.RDS    -  point transects with counts
samplers_no_crs.RDS    -  point transect locations
study_area_no_crs.RDS  -  study area polygon

create_mesh.R          -  creates INLA mesh and saves as mesh_no_crs.RDS 

### Analyses

Once mesh is created:

fit_model.R            -  fit the model, saved as akepa_20020_fitted_model.RDS

eval_model.R           -  model summary figures for the paper

excursions.R           -  create excursions plots for the paper
