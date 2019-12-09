# Evaluate spde by posterior sampling

library(INLA)
library(inlabru)
library(rgeos)
source("helper_functions.R")

# set ggplot theme
theme_set(theme_minimal())
set.seed(9701071)

#### WARNING - as it stands this script will overwrite figures for the paper

#### Fitted model  ####
file_name = "akepa_2002_fitted_model.RDS"
fit = readRDS(file_name)

study_area = readRDS("data/study_area_extended_no_crs.RDS")
samplers = readRDS("data/samplers_extended_no_crs.RDS")
mesh = readRDS("data/mesh_extended_no_crs.RDS")
obs = readRDS("data/obs_extended_no_crs.RDS")

### Generate posterior point patterns
n.pp = 25  

# simulate intensities
llam = do.call(cbind,
               generate(fit, mesh, ~ grf + Intercept, n.samples = n.pp))

# create bins and matrix to store counts
max.ds = max(as.numeric(gDistance(samplers, byid = TRUE)))
breaks = seq(0, ceiling(max.ds), length.out = 20)


# plot one point pattern
# a.pp = subset(post.pp, sample == 1)
# ggplot() +
#   gg(study_area) +
#   gg(a.pp) +
#   coord_equal()

# detection functiosn for thinning
# Log half-normal
log_hn = function(distance, lsig){ 
  -0.5*(distance/exp(lsig))^2
}

# Half-normal
hn <- function(distance, lsig) exp(log_hn(distance, lsig))

# thin.pp needs a plotID columns
samplers$plotID = samplers$SampleLabel

lsig.post = generate(fit, data.frame(x=1), ~ lsig, n.samples = n.pp)

# memory leakin somewhere in this loop - think rgeos function
# maybe in thin.pp() ?
# which time?

# please don't judge me for the code below.
# I just cant find out what's happening with the memory

counts = matrix(NA, nrow = n.pp, ncol = length(breaks)-1)

# okay run this a few times to get desired sample size.
# memory problems in the loops.
for (i in 1:(n.pp)){
  
  a.pp = sample.lgcp(mesh = mesh, loglambda = llam[,i])
  a.pp = a.pp[study_area,]
  a.lsig = lsig.post[[i]]
  
  a.obs = thin.pp(pp = a.pp, 
                  W = 60, 
                  samplers = samplers, 
                  snap = TRUE,
                  detfn = hn,
                  lsig = a.lsig)
  
  # calculate point to point distances of observations at samplers
  
  # freq for posterior sample
  a.ds = gDistance(a.obs, byid = TRUE)
  a.ds[upper.tri(a.ds, diag = TRUE)] = NA
  a.ds.vec = as.numeric(a.ds)
  counts[i,] = hist(a.ds.vec, breaks = breaks, plot = FALSE)$counts
  
  print(paste("Finished sample ", i))
  
  # memory?
  rm(a.ds, a.ds.vec, a.obs, a.lsig, a.pp)
  gc()
}

# RUN THIS IF ITS FIRST TIME RUNNING THE SCRIPT
# counts.all = counts
# saveRDS(counts.all, "data/post.pp.counts.RDS")

# Other times run this
counts.all = readRDS("data/post.pp.counts.RDS")
counts.all = rbind(counts.all, counts)
saveRDS(counts.all, "data/post.pp.counts.RDS")

# freq for observations
obs.ds = gDistance(obs, byid = TRUE)
obs.ds[upper.tri(obs.ds, diag = TRUE)] = NA
obs.ds.vec = as.numeric(obs.ds)
count.obs = hist(obs.ds.vec, breaks = breaks, plot = FALSE)$counts

# plot
png(filename = "../figures/post_pp_distances.png",
    width = 7, height = 4, units = "in", res = 100)

boxplot(counts, xlab = "distance", ylab = "count", xaxt = "n", range = 0)
points(1:19, count.obs, col = "red")
# xticks = floor(breaks[-1] - (breaks[2] - breaks[1])/2)
xticks = round(breaks[-1])
axis(1, at = 1:19, labels = xticks, las = 2)

dev.off()
