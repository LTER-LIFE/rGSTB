#' Build INLA mesh for spatial models
#'
#' Creates a 2D mesh based on plot and grid coordinates,
#' with a fixed minimal distance and offsets.
#'
#' @param xy.plots data.frame or matrix with columns X,Y (plot locations)
#' @param xy.grid data.frame or matrix with columns X,Y (grid locations)
#' @param mindist numeric, minimal allowed distance between points (default: 100)
#' @param max.edge numeric vector length 2 for inla.mesh.2d max.edge argument.
#'   If NULL, default c(mindist, 6 * mindist) is used.
#' @param offset numeric vector length 2 for inla.mesh.2d offset argument
#'   (inner and outer offset). Default: c(mindist, 6 * mindist).
#'
#' @return an \code{inla.mesh} object
#' @export
build_mesh <- function(xy.plots,
                       xy.grid,
                       mindist = 100,
                       max.edge = NULL,
                       offset = NULL) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required but not installed.")
  }

  if (is.null(max.edge)) {
    max.edge <- c(mindist, 6 * mindist)
  }
  if (is.null(offset)) {
    offset <- c(mindist, 6 * mindist)
  }

  xy.all <- rbind(as.matrix(xy.plots), as.matrix(xy.grid))

  mesh <- INLA::inla.mesh.2d(
    loc = xy.all,
    max.edge = max.edge,
    cutoff   = mindist,
    offset   = offset
  )

  mesh
}


#' Create SPDE model and projection matrices
#'
#' Given a mesh and coordinates for plots and grid, create
#' - SPDE object (Matern model) and
#' - A-matrices (projectors) for plots and grid.
#'
#' @param mesh an \code{inla.mesh} object as returned by \code{build_mesh()}
#' @param xy.plots data.frame or matrix with columns X,Y (plot locations)
#' @param xy.grid data.frame or matrix with columns X,Y (grid locations)
#' @param alpha numeric, smoothness parameter for SPDE (default: 2)
#'
#' @return list with elements:
#'   \item{spde}{SPDE model object}
#'   \item{A.plots}{projector matrix for plots}
#'   \item{A.grid}{projector matrix for grid}
#' @export
build_spde_and_projectors <- function(mesh,
                                      xy.plots,
                                      xy.grid,
                                      alpha = 2) {

  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Package 'INLA' is required but not installed.")
  }

  # SPDE model (Matern)
  spde <- INLA::inla.spde2.matern(
    mesh = mesh,
    alpha = alpha
  )

  # Projector matrices
  A.plots <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(xy.plots)
  )

  A.grid <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(xy.grid)
  )

  list(
    spde    = spde,
    A.plots = A.plots,
    A.grid  = A.grid
  )
}