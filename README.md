# Point Transect Distance Sampling With Inlabru
## A Case Study with Hawaiian Akepa Distance Sampling Data

### Data files

**obs.RDS**    -  point transects locations with observed distances

**samplers.RDS**    -  point transect locations

**study_area.RDS**  -  study area polygon

**mesh.RDS**  - INLA mesh 

**create_mesh.R** -  creates INLA mesh

### Analyses

Once mesh is created:

**fit_model.R**            -  fit the model

**eval_model.R**           -  model summary figures for the paper

**eval_spde.R** 					 -  compare observed pairwise distances with simulated point patterns

**excursions.R**           -  create excursions plots for the paper

**posterior_N.R**					 -  create posterior plot and investigate various approximation strategies for estimating abundance

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
