#' Embedding vectors for faces in the Labelled Faces in the Wild dataset
#'
#' A dataset of OpenFace embeddings for the "Labelled Faces in the Wild" dataset, see \url{http://vis-www.cs.umass.edu/lfw/}.
#'
#' OpenFace is a facial recognition library. The similarity between two OpenFace vectors should correlate with the
#' likelihood that the vectors were generated from images of the same person. For details and discussion,
#' see \url{https://cmusatyalab.github.io/openface/}.  The images may be obtained from \url{http://vis-www.cs.umass.edu/lfw/}.
#'
#' @format A data.frame where each row represents an image.  The first column is the name of the person in the image, the second column
#' is the name of the image file, and the remaining columns are the columns of the embedding vector for each image as calculated with
#' the OpenFace `batch-represent` function.
#'
#' @source \url{http://openface-models.storage.cmusatyalab.org/lfw.nn4.small2.v1/labels.csv}
#' @source \url{http://openface-models.storage.cmusatyalab.org/lfw.nn4.small2.v1/reps.csv}
"facevectors"
