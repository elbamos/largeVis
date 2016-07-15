#' Visualize an embedding by plotting with images
#'
#' Makes a plot of \code{n} images sampled from \code{images}, positions at coordinates given by \code{x}.
#'
#' @param x A \code{largeVis} object or [N,D] matrix of coordinates.
#' @param n The number of images to sample.
#' @param images The images. A 3-D or 4-D array.
#' @param scale Proportion to scale the images to.
#' @param ... Addiitional parameters passed to \code{plot}.
#'
#' @details The images can be passed in either as a list or a 3- or 4-dimensional array. The first dimension is \code{n}.
#'
#' If the objects in the list are \code{matrix} objects, or the array is 3-dimensional, the images will be treated as
#' greyscale. If there is an additional dimension, it must have a length of 3 and be RGB color layers.
#'
#' @references Andrej Karpapthy. \href{http://cs.stanford.edu/people/karpathy/cnnembed/}{t-SNE Visualization of CNN Codes.}
#'
#' @importFrom grDevices as.raster
#' @importFrom graphics rasterImage
#' @seealso \code{\link{ggManifoldMap}}
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
                        n = nrow(x),
                        images,
                        scale = 1,
                        ...) { #nocov start
  if (class(x) == "largeVis") x <- t(x$coords)
  if (ncol(x) != 2) stop("Can only visualize in 2-D.")
  N <- nrow(x)
  if (class(images) == "list" &&
      N != length(images))
    stop("Number of images doesn't equal number of points.")
  if (N != nrow(images))
    stop("Number of images doesn't equal number of points.")

  D <- length(dim(images)) - 1

  if (! (D == 2 || D == 3)) stop("Wrong number of dimensions.")
  if (D == 3 &&
      (dim(images)[4] < 2 ||
       dim(images)[4] > 4)) stop("Wrong number of color layers.")

  selections <- sample(N, n, replace = F)
  lowerscale <- min(images)
  upperscale <- max(images)
  graphics::plot(x * 1.1, pch = NA, type = 'n', ...)

  for (i in selections) {
    if (D == 2) {
      image_data <- images[i, , ]
    } else {
      image_data <- images[i, , , ]
    }
    image_data <- 1 - ( (image_data - lowerscale) / upperscale)
    image <- grDevices::as.raster(image_data)
    offsetx <- (nrow(image) * scale) / 2
    offsety <- (ncol(image) * scale) / 2
    graphics::rasterImage(image,
                          x[i, 1] - offsetx,
                          x[i, 2] - offsety,
                          x[i, 1] + offsetx,
                          x[i, 2] + offsety,
                          interpolate = TRUE
    )
  }
} # nocov end

#' Visualize an embedding by ggplotting with images
#'
#' Identical to \link{manifoldMap}, but adds images to an existing \code{ggplot2} object or creates one.
#'
#' @param ggObject a \code{\link[ggplot2]{ggplot}} object.  If not provided, a new \code{ggplot}
#' object with \code{\link[ggplot2]{geom_blank}} will be created.
#' @param x A \code{largeVis} object or [N,D] matrix of coordinates.
#' @param n The number of images to sample.
#' @param images The images. A 3-D or 4-D array.
#' @param scale Proportion to scale the images to.
#' @param ... Addiitional parameters passed to \code{plot}.
#' @return A \code{ggplot} object.
#'
#' @details See \code{\link{manifoldMap}}.  Note that this function can be considerably slower to display than \code{manifoldMap}.
#' It therefore should only be used if other features of \code{ggplot2} are required.
#'
#' If the objects in the list are \code{matrix} objects, or the array is 3-dimensional, the images will be treated as
#' greyscale. If there is an additional dimension, it must have a length of 3 and be RGB color layers.
#'
#' @importFrom grDevices as.raster
#' @importFrom graphics rasterImage
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 annotation_raster
#' @importFrom ggplot2 aes
#' @export
ggManifoldMap <- function(
                          ggObject = NULL,
                          x,
                          n = nrow(x),
                          images,
                          scale = 1,
                          ...) { #nocov start
  if (class(x) == "largeVis") x <- t(x$coords)
  if (ncol(x) != 2) stop("Can only visualize in 2-D.")
  N <- nrow(x)
  if (class(images) == "list" &&
      N != length(images))
    stop("Number of images doesn't equal number of points.")
  if (N != nrow(images))
    stop("Number of images doesn't equal number of points.")

  D <- length(dim(images)) - 1

  if (! (D == 2 || D == 3)) stop("Wrong number of dimensions.")
  if (D == 3 &&
      (dim(images)[4] < 2 ||
       dim(images)[4] > 4)) stop("Wrong number of color layers.")

  selections <- sample(N, n, replace = F)
  lowerscale <- min(images)
  upperscale <- max(images)
  if (is.null(ggObject)) {
    x <- data.frame(x)
    colnames(x) <- c("x", "y")
    ggObject = ggplot2::`%+%`(ggplot2::ggplot(x,
                                              ggplot2::aes_(x = quote(x),
                                                            y = quote(y))),
                              geom_blank())
  }

  for (i in selections) {
    if (D == 2) {
      image_data <- images[i, , ]
    } else {
      image_data <- images[i, , , ]
    }
    image_data <- 1 - ( (image_data - lowerscale) / upperscale)
    image <- grDevices::as.raster(image_data)
    offsetx <- (nrow(image) * scale) / 2
    offsety <- (ncol(image) * scale) / 2

    ggObject <- ggplot2::`%+%`(ggObject, ggplot2::annotation_raster(
      image,
      xmin = x[i, 1] - offsetx,
      ymin = x[i, 2] - offsety,
      xmax = x[i, 1] + offsetx,
      ymax = x[i, 2] + offsety,
      interpolate = TRUE
    ))
  }
  return(ggObject)
} # nocov end

#' manifoldMapStretch
#'
#' A manifold map that fills the full extent of the plot.
#'
#' Ported from \url{http://cs.stanford.edu/people/karpathy/cnnembed/}.  Each position is filled with its nearest neighbor.
#'
#' @param x A [N,D] matrix of coordinates.
#' @param f A function that, called with the index number of a row of \code{x}, returns an R object representing
#' an image. See the example.
#' @param size_x The width of the requested plot, in pixels.
#' @param size_y The height of the requested plot, in pixels.
#' @param image_size The size to plot each image; each is plotted as a square.
#' @param ... Additional parameters passed to \code{plot}.
#'
#' @note This function is experimental.
#'
#' @examples
#' \dontrun{
#' # Demonstration of f
#' load(system.file("extdata", "faces.Rda", package="largeVis"))
#'
#' imagepaths <- paste("pathtoimages",
#'    faceLabels[,1], sub("png", "jpg", faceLabels[,2]), sep = "/")
#'
#' manifoldMapStretch(as.matrix(faceCoords[,1:2]),
#'    f = function(x) jpeg::readJPEG(imagePaths[x]),
#'    size_x = 5000, size_y = 5000, image_size = 100)
#' }
#'
#' @export
manifoldMapStretch <- function(x,
                               f,
                               size_x = 500,
                               size_y = 500,
                               image_size = 50,
                               ...) { #nocov start

  xnum <- size_x / image_size
  ynum <- size_y / image_size

  coordsadj <-   x  - c(min(x[,1]), min(x[,2]))
  coordsadj <- coordsadj * c(size_x / max(coordsadj[,1]),
                             size_y / max(coordsadj[,2]))

  graphics::plot(matrix(c(0, size_x,
                          size_y, 0),
                        ncol = 2),
                 pch = NA,
                 type = 'n', ...)

  abes <- matrix(c(rep(1:xnum, ynum),
                   rep(1:ynum, each = xnum)),
                 ncol = 2)

  for (i in 1:nrow(abes)) {
    img_x <- abes[i, 1]
    img_y <- abes[i, 2]
    xf <- (img_x * image_size) - (image_size/2)
    yf <- (img_y * image_size) - (image_size/2)
    dd <- apply((coordsadj - c(xf, yf))^2,
                MARGIN = 1,
                FUN = sum)
    selection <- which(dd == min(dd))
    coordsadj[selection,] <- Inf
    image <- f(selection)
    rasterImage( image,
                 xleft = xf - (image_size / 2),
                 ybottom = yf - (image_size / 2),
                 xright = xf + (image_size / 2),
                 ytop = yf + (image_size / 2),
                 interpolate = TRUE
    )
  }
} #nocov end
