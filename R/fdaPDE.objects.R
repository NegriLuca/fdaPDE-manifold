#' Create a FEM basis
#' 
#' @param mesh A \code{MESH2D} or \code{MESH.2.5D} object representing the domain triangulation. See \link{create.MESH.2D}, \link{create.MESH.2.5D}.
#' @return A  \code{FEMbasis} object. This contains the \code{mesh}, along with some additional quantities:
#'
#' if \code{class(mesh) == MESH2D}
#' 	\item{\code{order}}{Either "1" or "2". Order of the Finite Element basis.} 
#' 	\item{\code{nbasis}}{Scalar. The number of basis.} 
#' 	\item{\code{transf_coord}}{An object containing 4 vectors of length #triangles. The for of them encode the tranformation matrix [diff1x diff2x; diff1y diff2y] that transforms the nodes of the reference triangle to the nodes of the i-th triangle.}
#' 	\item{\code{detJ}}{A vector of length #triangles. The ith element contains the determinant of the transformation from the reference triangle to the nodes of the i-th triangle. It's values is also the double of the area of each triangle of the basis.}
#' if \code{class(mesh) == MESH.2.5D}
#' 	\item{\code{order}}{Either "1" or "2". Order of the Finite Element basis.}
#' 	\item{\code{nbasis}}{Scalar. The number of basis.}
#' @description Sets up a Finite Element basis. It requires a triangular mesh, a \code{MESH2D} or a \code{MESH.2.5D} object, as input. 
#' The basis' functions are globally continuos surfaces, that are polynomials once restricted to a triangle in the mesh. 
#' Linear if (\code{order = 1}) in the input \code{mesh} and quadratic if (\code{order = 2}) in the input \code{mesh}
#' Finite Element are currently implemented.
#' @usage create.FEM.basis(mesh)
#' @seealso \code{\link{create.MESH.2D}}, \code{\link{create.MESH.2.5D}}
#' @examples 
#' ## Creates a simple triangulated domain with a concavity; this is a MESH2D object  
#' mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
#' segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order=1)
#' ## Plot it
#' plot(mesh)                   
#' ## Creates the basis
#' FEMbasis = create.FEM.basis(mesh)

create.FEM.basis = function(mesh)
{
  if (class(mesh)=="MESH2D"){

	  #  The number of basis functions corresponds to the number of vertices
	  #  for order = 1, and to vertices plus edge midpoints for order = 2

	  nbasis = dim(mesh$nodes)[[1]]
	  eleProp = R_elementProperties(mesh)

	  #eleProp = NULL
	  #if(CPP_CODE == FALSE)
	  #{
	  #  eleProp = R_elementProperties(mesh)
	  #}
  
  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order), nbasis = nbasis, detJ=eleProp$detJ, transf_coord = eleProp$transf_coord)
  class(FEMbasis) = "FEMbasis"
  
  FEMbasis
  } else if (class(mesh) == "MESH.2.5D"){

  	  FEMbasis = list(mesh = mesh, order = as.integer(mesh$order),nbasis = mesh$nnodes)
  	  class(FEMbasis) = "FEMbasis"
  	  FEMbasis
  }
}


#' Define a surface or spatial field by a Finite Element basis expansion
#' 
#' @param coeff A vector or a matrix containing the coefficients for the Finite Element basis expansion. The number of rows 
#' (or the vector's length) corresponds to the number of basis in \code{FEMbasis}. 
#' The number of columns corresponds to the number of functional replicates. 
#' @param FEMbasis A \code{FEMbasis} object defining the Finite Element basis, created by \link{create.FEM.basis}.
#' @description This function defines a FEM object. This is not usualled called directly by users.
#' @usage FEM(coeff,FEMbasis)
#' @return An \code{FEM} object. This contains a list with components \code{coeff} and \code{FEMbasis}.
#' @examples 
#' ## Upload a triangular mesh and plot it
#' data("mesh.2D.rectangular")
#' plot(mesh.2D.rectangular)
#' ## Create a linear Finite Element basis
#' FEMbasis = create.FEM.basis(mesh.2D.rectangular)
#' ## Define a sinusoidal function as expansion of this basis and plot it
#' coeff <- sin(mesh.2D.rectangular$nodes[,1])*cos(mesh.2D.rectangular$nodes[,2])
#' FEM_object<- FEM(coeff, FEMbasis)
#' plot(FEM_object)

FEM<-function(coeff,FEMbasis)
{
  if (is.null(coeff)) 
    stop("coeff required;  is NULL.")
  if (is.null(FEMbasis)) 
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis) != "FEMbasis")
    stop("FEMbasis not of class 'FEMbasis'")
  coeff = as.matrix(coeff)
  if(nrow(coeff) != FEMbasis$nbasis)
    stop("Number of row of 'coeff' different from number of basis")
  
  fclass = NULL
  fclass = list(coeff=coeff, FEMbasis=FEMbasis)
  class(fclass)<-"FEM"
  return(fclass)
}
