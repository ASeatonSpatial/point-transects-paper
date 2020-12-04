# library(ggpolypath)

# update gg.SpatialPixelsDataFrame to play well with factor
# this could be improved by allowing pixels to have multiple columns, only do default
# behaviour if length(names(data)) == 1, otherwise make user specify

# the long lat default name behaviour comes from the fortify.sp function
# by default this labels coords long lat, regardless of crs
# I assume this is why Fabian just kept long lat as default for sp objects
# That this behaviour comes from fortify means that we cannot use the aes(...) argument
# to provide the names of the coordinates.  The created dataframe will have columns
# long and lat.
# So seems like easiest thing is to just just xlab and ylab to get labels you want

gg.SpatialPixelsDataFrame <- function(data,
                                      mapping = NULL,
                                      alpha = NULL,
                                      crs = NULL,
                                      mask = NULL, ...) {
  if (!is.null(crs)) {
    data <- spTransform(data, crs)
  }
  if (!is.null(mask)) {
    data <- data[as.vector(!is.na(over(data, mask))), ]
  }

  df <- as.data.frame(data)

  if (!is.factor(data[,c(names(data)[[1]])])){
  dmap <- aes_string(
    x = coordnames(data)[1],
    y = coordnames(data)[2],
    fill = names(data)[[1]]
  )} else
    dmap <- aes_string(
      x = coordnames(data)[1],
      y = coordnames(data)[2],
      group = names(data)[[1]]
    )

  if (!is.null(mapping)) {
    dmap <- modifyList(dmap, mapping)
  }
  if (!is.null(alpha)) dmap <- modifyList(dmap, aes_string(alpha = alpha))
  gm <- geom_tile(data = df, mapping = dmap, ...)

  # If data is not discrete (factor), add default color palette
  # if (!inherits(data[[deparse(dmap$fill)]], "factor")) {
  #  gm = c(gm, scale_fill_gradientn(colours = bru.pal()))
  # }
  gm
}


#' # Update this so plotting identifies holes in polygons
#' # replace geom_polygon with ggpolypath::geom_polypath
#'
#' #' Geom for SpatialPolygons objects
#' #'
#' #' Uses the fortitfy() function to turn the SpatialPolygons objects into a data.frame. Then
#' #' calls [geom_polygon] to plot the polygons.
#' #'
#' #'
#' #' @aliases gg.SpatialPolygons
#' #' @name gg.SpatialPolygons
#' #' @export
#' #' @import ggplot2
#' #' @param data A SpatialPolygons object.
#' #' @param mapping Aesthetic mappings created by [aes] or [aes_] used to update the default
#' #'                mapping. The default mapping is `aes_string(x = "long", y = "lat", group = "group")`.
#' #' @param crs A [CRS] object defining the coordinate system to project the data to before plotting.
#' #' @param color Filling color for the polygons.
#' #' @param alpha Alpha level for polygon filling.
#' #' @param ... Arguments passed on to [geom_polygon].
#' #' @return A [geom_polygon] return value.
#' #' @family geomes for spatial data
#' #' @example inst/examples/gg.spatial.R
#'
#' gg.SpatialPolygons <- function(data, mapping = NULL, crs = NULL, color = "black", alpha = NULL, ...) {
#'   if (!is.null(crs)) {
#'     data <- spTransform(data, crs)
#'   }
#'   df <- fortify(data)
#'   dmap <- aes_string(x = "long", y = "lat", group = "group")
#'
#'   if (!("alpha" %in% names(dmap)) & is.null(alpha)) {
#'     alpha <- 0.1
#'   }
#'   if (!("color" %in% names(dmap)) & is.null(color)) {
#'     color <- "black"
#'   }
#'
#'   if (!is.null(mapping)) {
#'     dmap <- modifyList(dmap, mapping)
#'   }
#'
#'   geom_polypath(data = df, mapping = dmap, alpha = alpha, color = color, ...)
#' }

