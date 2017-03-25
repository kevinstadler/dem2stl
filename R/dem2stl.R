#' Map latitude/longitude coordinates onto an ellipsoid. Vectorised.
#'
#' @param lat latitude coordinates (scalar or vector, in radians)
#' @param long longitude coordinates (same length as \code{lat}, in radians)
#' @param a semi-principal axis length
#' @param b semi-principal axis length
#' @param c semi-principal axis length
#' @return a \code{length(lat)x3} matrix giving the lat/long coordinates points
#'   on an ellipsoid with the given axes lengths in 3d space
#' @export
latlong2ellipsoid <- function(lat, long, a=1, b=1, c=1)
  matrix(c(a*cos(lat)*cos(long),
           b*cos(lat)*sin(long),
           c*sin(lat)), ncol=3)

#' @export
deg2rad <- function(deg)
  (deg * pi) / 180.0

#' @export
rad2deg <- function(rad)
  (rad * 180.0) / pi

#' Calculate vector crossproduct.
#'
#' @param ... \code{n} vectors, each of length \code{n+1}
#' @return a vector of the same length as the \code{...} arguments
#' @export
xprod <- function(...) {
  m <- do.call(rbind, list(...))
  sapply(seq(ncol(m)), function(i) det(m[, -i, drop=FALSE]) * (-1)^(i+1))
}

#' Calculate the surface normal vector at a specific point on an ellipsoid.
#'
#' @param lat latitude coordinates (in radians)
#' @param long longitude coordinates (in radians)
#' @param a semi-principal axis length of the ellipsoid
#' @param b semi-principal axis length of the ellipsoid
#' @param c semi-principal axis length of the ellipsoid
#' @export
ellipsoidnormal <- function(lat, long, a=1, b=1, c=1) {
  # https://en.wikipedia.org/wiki/Ellipsoid#Parametric_representation
  n <- c(b*c*cos(lat)*cos(long),
         a*c*cos(lat)*sin(long),
         a*b*sin(lat))
  n / sqrt(sum(n^2))
}

# returns a logical matrix the same size as data
partofbiggestclump <- function(data) {
  # find biggest 4-connected clump
  clump <- raster::clump(data, directions=4, gaps=FALSE)
  # check if cells are part of biggest clump
  clump[] == (which.max(raster::freq(clump, useNA="no")[ ,2])[1])
}

# apply a 2x2 convolution filter on the given (matrix) data that erases all
# except the biggest connected clump
edgefilter <- function(data) {
  complete <- !is.na(diff(t(diff(data))))
  complete[!partofbiggestclump(raster::raster(t(complete)))] <- FALSE
  data[!as.logical(t(cbind(0, rbind(0, complete)) +
                     cbind(0, rbind(complete, 0)) +
                     cbind(rbind(0, complete), 0) +
                     cbind(rbind(complete, 0), 0)))] <- NA
  data
}

boundaryoutline <- function(data) {
  # temporarily extend by a buffer row of NAs to get *all* boundaries
  boundary <- which(!(is.na(raster::crop(
    raster::boundaries(raster::extend(data, 1), asNA=TRUE), data)[])))
  # incremental trace - start in top left corner, go clockwise
  outline <- boundary[1]
  cell <- raster::rowColFromCell(data, outline)
  direction <- 2 # 0 = west, 1 = north, 2 = east, 3 = south
  for (i in seq(boundary)) {
    # clockwise -> always try to turn left first
    for (d in (direction + -1:1) %% 4) {
      neighbour <- cell + c(sign(d-2) * d%%2, sign(d-1) * (d-1)%%2)
      if (!is.na(data[neighbour[1], neighbour[2]])) {
        outline[i+1] <- raster::cellFromRowCol(data, neighbour[1], neighbour[2])
        direction <- d
        cell <- neighbour
        break
      }
    }
  }
  return(outline)
}

#' Create a 3d mesh from raster data.
#' 
#' Creates a 3d mesh showing elevations and earth curvature based on raster
#' data.
#' 
#' This function can render arbitrarily shaped landmasses by setting some
#' values of the raster data to \code{NA}, which are consequently not rendered.
#' Only the biggest continuous clump of non-\code{NA} values in the data set
#' is rendered. The center points of each raster cell are mapped to vertices
#' in 3d space on a spheroid with the given principal axes. Protrusions out of
#' the biggest clump which are only one raster cell wide are shaved off
#' (because a single protruding vertex cannot form part of a polygon if it
#' doesn't have neighbours on two adjacent sides).
#'
#' The function does not support rendering of holes: holes (\code{NA} regions
#' within the biggest clump's interior) are filled with 0s.
#' 
#' @param data a \code{raster} object (such as of class \code{RasterLayer})
#' @param thicknessratio desired thickness of the base as a proportion of
#'   the overall size of the spheroid (i.e. between 0 and 1)
#' @param thickness desired thickness of the base (in absolute 3d coordinates)
#' @param size desired maximum extent (width/height) of the resultin mesh in mm
#' @param scale scale at which mesh should be rendered. If scale is specified,
#    the \code{size} argument is ignored.
#' @param exaggeration factor by which elevation differences are exaggerated
#' @param up desired direction of the front face of the resulting mesh
#' @param flatbase logical
#' @param a semi-principal axis length of the Earth spheroid
#' @param b semi-principal axis length of the Earth spheroid
#' @export
dem2mesh <- function(data, thicknessratio=0.01, thickness=NULL, size=NULL, scale=NULL, exaggeration=1, up=c(0, 0, 1), flatbase=FALSE, a=6378137.0, b=6356752.314245) {
  # clean up data (take only biggest 4-connect clump, shave protrusions off)
  data[] <- edgefilter(raster::as.matrix(data))
  data <- raster::trim(data)

  # fill holes (otherwise we'd need more sidewalls)
  holes <- raster::extend(data, 1)
  holes[] <- as.numeric(is.na(holes[]))
  holes <- raster::clump(holes, useNA=FALSE)
  # don't fill in surrounding
  holes[holes[] == holes[1]] <- NA
  data[which(!is.na(raster::crop(holes, data)[]))] <- 0


  # create top vertices
  coords <- deg2rad(raster::xyFromCell(data, 1:length(data)))
  # take spheroid and calculate position of surface in 3d space
  surface <- latlong2ellipsoid(coords[,"y"], coords[,"x"], a, a, b)

  validdata <- matrix(!(is.na(data[])), nrow=ncol(data))
  # offset coordinates by elevation differences
  # increment by normalize(latlong unit vector) * values * exaggeration
  surface[validdata,] <- surface[validdata,] *
    (1 + raster::values(data)[which(validdata)] * exaggeration / sqrt(rowSums(surface[validdata,]^2)))

  # calculate quad polygon corner indices (referring to surface matrix columns)
  topleftcorners <- validdata & rbind(validdata[-1, ], F) &
                                cbind(validdata[, -1], F) & 
                                cbind(rbind(validdata[-1, -1], F), F)
  corners <- sapply(which(topleftcorners),
    function(i) c(i+1, i, i+ncol(data), i+ncol(data)+1))

  # add side+bottom to make whole object
  if (is.null(thickness)) {
    # calculate thickness (proportion of ellipsoid crust included) from
    # thickness ratio (in relation to maximum extent of the object)
    thickness <- thicknessratio
  }

  center <- deg2rad(apply(matrix(raster::extent(data), nrow=2), 2, mean))
  tangentpoint <- as.vector(latlong2ellipsoid(center[2], center[1], a, a, b))
  normal <- ellipsoidnormal(center[2], center[1], a, a, b)

  # convert to indices
  validdata <- which(validdata)

  # create base plane (or curve)
  if (flatbase) {
    # project shape outline onto a base plane which is parallel to the tangent
    # plane on the center of the shape
    # determine distance from the lowest (backest) point of our figure to the
    # tangent plane of its frontest point
    # http://mathinsight.org/distance_point_plane
    depth <- max(crossprod(tangentpoint - t(surface[validdata, ]), normal))
    baseplanecenter <- tangentpoint - normal*depth*7
    basepoints <- surface - tcrossprod(crossprod(t(surface) - baseplanecenter, normal), normal)
  } else {
    # take lat/long-equivalent surface points of an ellipsoid proportionally
    # smaller than the surface one
    thickness <- 1 - thickness
    basepoints <- matrix(latlong2ellipsoid(coords[,"y"], coords[,"x"], a*thickness, a*thickness, b*thickness), ncol=3)
    basepoints[-validdata,] <- NA
  }

  x <- !complete.cases(surface)
  surface <- rbind(surface, basepoints)
  corners <- cbind(corners, length(data) + corners[4:1,])

  # connect top and base vertices
  outline <- boundaryoutline(data)
  corners <- cbind(corners,
    sapply(seq(length(outline)-1),
      function(i) c(outline[i], outline[i+1], length(data)+outline[i+1], length(data)+outline[i])))

  # remove NA vertices, adjust corner indices accordingly
  # should rather remove all vertices not referred to by indices?
  realvertices <- complete.cases(surface)
  surface <- surface[realvertices,]
  corners <- cumsum(realvertices)[corners]

  # rotate the mesh so that the normal at its center points up
  theta <- acos( (normal %*% up)[1] / (sqrt(sum(normal^2)) * sqrt(sum(up^2))) )
  rotaxis <- xprod(normal, up)
  rotmtrx <- rgl::rotationMatrix(theta, rotaxis[1], rotaxis[2], rotaxis[3])[1:3, 1:3]
  surface <- apply(surface, 1, function(vertex) rotmtrx %*% vertex)

  if (is.null(scale)) {
    maxextent <- max(diff(apply(surface, 1, range)))
    # set maximum extent to about 10cm in diameter
    # TODO use size argument
    # if (is.null(size)) {
    #} else {
    scale <- maxextent/size
    #}
  } else if (!is.null(size)) {
    warn("Both size and scale are specified, ignoring size")
  }
  # apply scale by making homogeneous coords with w=scale
  # (input is in meters, output in cm)
  message("Rendering at scale 1:", 100*scale)
  rgl::qmesh3d(rbind(surface, scale), corners)
}

writemesh <- function(writefun, filename, x, ...) {
  rgl::open3d(useNULL=TRUE)
  rgl::plot3d(x, box=FALSE, axes=FALSE, xlab=NULL, ylab=NULL, zlab=NULL)
  writefun(filename, ...)
}

#' Write mesh objects to various file formats.
#'
#' @param filename the file to write to
#' @param x an rgl \code{mesh3d} object
#' @param ... extra arguments passed to the corresponding
#'   \code{rgl::write...()} function
#' @return Invisibly returns the name of the connection to which the data was
#'   written.
#' @export
mesh2stl <- function(filename, x, ...)
  writemesh(rgl::writeSTL, filename, x, ...)

#' @rdname mesh2stl
mesh2obj <- function(filename, x, ...)
  writemesh(rgl::writeOBJ, filename, x, ...)

#' @rdname mesh2stl
mesh2ply <- function(filename, x, ...)
  writemesh(rgl::writePLY, filename, x, ...)

#' @rdname mesh2stl
mesh2webgl <- function(filename, x, ...)
  writemesh(rgl::writeWebGL, filename, x, ...)
