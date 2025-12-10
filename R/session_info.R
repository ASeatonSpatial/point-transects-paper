# session info for results in paper

library(ggplot2)
library(sp)
library(inlabru)
library(INLA)
library(scales)
library(cowplot)
library(sf)
library(patchwork)
library(dplyr)
library(excursions)

writeLines(
  capture.output(devtools::session_info()),
  here::here("R", "session_info.txt")
)
