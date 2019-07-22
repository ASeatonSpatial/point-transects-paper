# Evaluate fitted model
library(INLA)
library(inlabru)

#### Fitted model  ####
file_name = "../akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

source("AS.Nposterior.R") 
predpts = ipoints(study_area, mesh)

# system.time(AS.Nposterior(fit, predpts, n.samples = 4000))
# abnc = AS.Nposterior(fit, predpts, n.samples = 1000)   # n.samples = 1000 default

n.samples.list = c(1000, 2000, 4000, 8000, 16000)

for (n in n.samples.list){
  abnc = AS.Nposterior(fit, predpts, n.samples = n)   # n.samples = 1000 default
  g = ggplot() + 
    gg(abnc$Npost) +
    ggtitle(paste("n = ", n)) +
    xlim(3000, 9000)
  print(g)
}


  