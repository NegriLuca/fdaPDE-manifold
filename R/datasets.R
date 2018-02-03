#' Meuse river data set
#'
#' This data set gives locations, top soil heavy metal concentrations (ppm) and other variables, collected in a flood plain of the river Meuse. 
#' For details on this dataset, see \code{meuse.all} in \url{https://cran.r-project.org/web/packages/gstat/gstat.pdf}. This version of the dataset
#' includes 155 observations. Moreover, a definition of the domain's boundary is provided through the file {\link{MeuseBorder}}, 
#' as it is used in many examples to illustrate fdaPDE features. 
#' \itemize{
#'   \item sample: original sample number. In this version of the dataset some observations has been left out because not indicative, or outliers.
#'   \item x: numeric vector indicating the x-coordinate (m) in RDM (Dutch topographical map coordinates).
#'   \item y: numeric vector indicating the y-coordinate (m) in RDM (Dutch topographical map coordinates).
#'   \item cadmium: topsoil cadmium concentration (ppm).
#'   \item copper: topsoil copper concentration (ppm).
#'   \item lead: topsoil lead concentration (ppm).
#'   \item zinc: topsoil zinc concentration (ppm).
#'   \item elev: relative elevation.
#'   \item dist.log(m): logarithm of the distance to river Meuse (metres), as obtained during the field survey.
#'   \item om: organic matter, as percentage.
#'   \item ffreq: flooding frequency class.
#'   \item soil: soil type.
#' }
#'
#' @format A data frame with 155 rows and 12 variables.
#' @source \url{https://cran.r-project.org/web/packages/gstat/gstat.pdf}
#' \url{http://spatial-analyst.net/book/meusegrids}
#' @name MeuseData
NULL

#' Boundary of the Meuse River data set
#'
#' This file provides the boundary of the domain of the {\link{MeuseData}}.
#'
#' \itemize{
#'   \item vertex_start. A vector having as entries the row's indices of the location in {\link{MeuseData}} where a boundary segment starts from.
#'   \item vertex_end. A vector having as entries the row's indices of the location in {\link{MeuseData}} where a boundary segment ends to.
#' }
#'
#' @format A data frame with 52 rows and 2 variables.
#' @name MeuseBorder
NULL

#' Simple mesh
#'
#' A simple mesh. This is a MESH2D object created with \code{create.MESH.2D}.
#'
#' @name mesh.2D.simple
NULL

#' Simple Rectangular mesh
#'
#' A simple rectangular mesh. This is a MESH2D object created with \code{create.MESH.2D}.
#'
#' @name mesh.2D.rectangular
NULL


#' Hub mesh
#'
#'A triangular mesh of a manifold representing a hub geometry.This is a \code{MESH.2.5D} object created with \code{create.MESH.2.5D}.
#'
#'@name hub
NULL

#' Candy mesh
#'
#'A triangular mesh of a manifold representing a candy geometry.This is a \code{MESH.2.5D} object created with \code{create.MESH.2.5D}.
#'
#'@name candy
NULL

#' SphereData
#'
#' A dataset with information about the connectivity matrix and the nodes locations of a sphere geometry. It containes:
#' \itemize{
#' 	\item nodes. A #nodes-by-3 matrix specifying the locations of each node.
#' 	\item triangles. A #triangles-by-3*order matrix specifying the indices of the nodes in each triangle.
#'	   }
#' This dataset is used to create a \code{MESH.2.5D} object with the function \code{create.MESH.2.5D} 
#'
#'@name sphereData
NULL

#' Sphere3Ddata
#'
#' A dataset with information about the connectivity matrix and the nodes locations of a sphere geometry. It containes:
#' \itemize{
#' 	\item nodes. A #nodes-by-3 matrix specifying the locations of each node.
#' 	\item tetrahedrons. A #tetrahedrons-by-4 matrix specifying the indices of the nodes in each tetrahedron.
#'	   }
#' This dataset is used to create a \code{MESH.3D} object with the function \code{create.MESH.3D} 
#'
#'@name sphere3Ddata
NULL

#' Sphere3D
#'
#'A tetrahedral mesh of a volume representing a sphere geometry.This is a \code{MESH.3D} object created with \code{create.MESH.3D}.
#'
#'@name sphere3D
NULL
