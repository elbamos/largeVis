#' Visualize an embedding by plotting with images
#'
#' Makes a plot of \code{n} images sampled from \code{images}, positions at coordinates given by \code{x}.
#'
#' @param x A \code{largeVis} object or [N,D] matrix of coordinates.
#' @param n The number of images to sample.
#' @param images The images. Either a list of image data objects, or an array.
#' @param scale Proportion to scale the images to.
#' @param transparency Whether to add an alpha channel to greyscale images.
#' @param ... Addiitional parameters passed to \code{plot}.
#'
#' @details The images can be passed in either as a list or a 3- or 4-dimensional array. The first dimension is \code{n}.
#'
#' If the objects in the list are \code{matrix} objects, or the array is 3-dimensional, the images will be treated as
#' greyscale. If there is an additional dimension, it must have a length of 3 and be RGB color layers.
#'
#' @importFrom grDevices as.raster
#' @importFrom graphics rasterImage
#' @export
#' @examples \dontrun{
#' load("mnist.Rda")
#' load("mnistcoordinates.Rda")
#'
#' flip <- function(x) apply(x,2,rev)
#' rotate <- function(x) t(flip(x))
#'
#' mnistimages <- apply(mnist$images,
#'    MARGIN=1,
#'    FUN = function(x) as.array(rotate(flip(x))))
#' mnistimages <- t(mnistimages)
#' dim(mnistimages) <- c(42000, 28, 28)
#'
#' manifoldMap(coords,
#'    1000,
#'    scale = 0.07,
#'    mnistimages)
#' }

manifoldMap <- function(x,
                      n,
                      images,
                      scale = 1,
                      transparency = FALSE,
                      ...) {
  if (class(x) == 'largeVis') x <- t(x$coords)
  if (ncol(x) != 2) stop("Can only visualize in 2-D.")
  N <- nrow(x)
  if (class(images) == 'list' && N != length(images)) stop("Number of images doesn't equal number of points.")
  if (N != nrow(images)) stop("Number of images doesn't equal number of points.")

  D <- length(dim(images)) -1

  if (! (D == 2 || D == 3)) stop("Wrong number of dimensions.")
  if (D == 3 && (dim(x)[3] < 2 || dim(x)[3] > 4)) stop("Wrong number of color layers.")

  selections <- sample(N, n, replace = F)
  lowerscale <- min(images)
  upperscale <- max(images)
  plot(x[selections,], pch = NA, ...)

  for (i in selections) {
    if (class(images) == 'list') imageData <- images[[i]]
    else if (D == 2) {
      imageData <- images[i,,]
    } else {
      imageData <- images[i,,,]
    }
    imageData <- 1 - ((imageData - lowerscale) / upperscale)
    if (transparency) {
      if (length(dim(imageData)) == 2)
        imageData <- abind::abind(imageData,
                                  imageData,
                                  imageData,
                                  imageData,
                                  along = 3)
      else if (length(dim(imageData)) == 3) {
        alpha <- apply(imageData, MARGIN, FUN = sum)
        alpha <- alpha / max(alpha)
        imageData <- abind::abind(imageData, alpha, along = 3)
      }
    }
    image <- grDevices::as.raster(imageData)
    offsetx <- (nrow(image) * scale) / 2
    offsety <- (ncol(image) * scale) / 2
    graphics::rasterImage(image,
                          x[i,1] - offsetx,
                          x[i,2] - offsety,
                          x[i,1] + offsetx,
                          x[i,2] + offsety,
                          interpolate = TRUE
                          )
  }
}

rotate <- function(x) t(apply(x, 2, rev))
