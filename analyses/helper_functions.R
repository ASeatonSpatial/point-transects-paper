#' thin.pp()
#'
#' Thin a point pattern using a known detection function and survey design
#' Point transects only for now - add line transects in future
#'
#' Assumes each observation falls within only 1 transect
#'
#' TODO:
#' Change hardcoded "PlotID"
#' Avoid for loop
#'
#' @param pp point pattern SpatialPointsDataFrame
#' @param samplers SpatialPointsDataFrame of point transect locations   - MUST HAVE column called PlotID, unique identifier of each transect
#' @param W radius or half strip width
#' @param detfn  detection function.  First argument must be the vector of distances, second argument list of params
#' @param snap (default = TRUE) snap sampled points to centre of transects (or to the centre line).  To snap need to have plotID column
#' @param ... arguments passed to detfn
#' @return spatial points dataframe
#' @import sp
#' @export


library(dplyr)
library(sf)

thin.pp <- function(pp, samplers, W, detfn, snap = TRUE, ...){
  
  if (!("plotID" %in% colnames(samplers@data))) stop("samplers needs a column called plotID")
  samplers_area <- rgeos::gBuffer(samplers, byid = TRUE, width = W)
  pp_det <- rgeos::gIntersection(pp, samplers_area, byid = TRUE)  # detectable points
  
  # Associate each nest with a sampler
  IDs <- sp::over(pp_det, samplers_area)
  pp_det$plotID <- IDs$plotID
  
  # Calculate distances:
  distances <- rep(NA, nrow(pp_det))
  
  for (i in 1:nrow(pp_det)) {
    
    pt <- pp_det[i,]
    transect <- subset(samplers, plotID == pt$plotID)
    distances[i] <- rgeos::gDistance(pt, transect)
  }
  
  # add distances column and create prob detected column
  pdet <- detfn(distances, ...)
  pp_det$distance <- distances
  pp_det$pdet <- pdet
  
  # subsample detectable to get observed
  detected <- rbinom(nrow(pp_det), 1, pdet)
  pp_obs <- subset(pp_det, as.logical(detected))
  
  # snap to centre?
  if (snap) {samplers_sf <- sf::st_as_sf(samplers)  # convert sp to sf
  test <- dplyr::inner_join(samplers_sf, pp_obs@data, by = "plotID")
  pp_obs <- sf::as_Spatial(test)  # convert back to sp
  } else NULL
  
  return(pp_obs)
}