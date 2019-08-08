# Point Transect Distance Sampling With Inlabru
## A Case Study with Hawaiian Akepa Distance Sampling Data

A brief tour of this repository:

### Data files

**obs_2002_no_crs.RDS**    -  point transects locations with observed distances

**samplers_no_crs.RDS**    -  point transect locations

**study_area_no_crs.RDS**  -  study area polygon

**create_mesh.R**          -  creates INLA mesh and saves as mesh_no_crs.RDS 

### Analyses

Once mesh is created:

**fit_model.R**            -  fit the model, saved as akepa_20020_fitted_model.RDS

**eval_model.R**           -  model summary figures for the paper

**excursions.R**           -  create excursions plots for the paper

### Figures

All figures for the manuscript should be committed in the /figures directory  
Please ensure figure size is not too large before committing (some default resolutions for plots can be enormous)  

### Paper

Bibliography items in bibtex format are in **paper.bib**    
Please try to keep in alphabetical order and keep similar naming convention for the labels   
Only commit the **paper.bib** and **paper.tex** files.  Do not commit auxiliary tex files or the produced paper.pdf    

### General Git Reminders

Always pull from remote before pushing changes   
Before you pull check you have no uncommitted staged changes (this will stop automatic merge)   
Only commit files that are needed for the paper (figures, R scripts, .tex and .bib files should be about it)   
