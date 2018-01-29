#' Create a \code{MESH.3D} object from the connectivty matrix and nodes locations
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node
#' @param tetrahedrons A #tetrahedrons-by-4*order matrix specifying the indices of the nodes in each tetrahedrons
#' @param order Order of the Finite Element basis default is order = 1. Only order = 1 is possible for the 3D Finite Elements.
#' @return An object of the class \code{MESH.3D} with the following output:
#' \item{\code{nnodes}}{The #nodes contained in the mesh}
#' \item{\code{ntetrahedrons}}{The #tetrahedrons contained in the mesh}
#' \item{\code{nodes}}{A vector of length #nodes*3 containing the x,y and z coordinate for each point of the mesh}
#' \item{\code{tetrahedrons}}{A vector of length #tetrahedrons*4*order specifying the indices of the nodes in each triangle of the mesh}
#' \item{\code{order}}{It specifies the order of the Finite Element basis. When order = 1, each mesh tetrahedron is represented by 4 nodes (the tetrahedron vertices). 
#' These are respectively used for linear (order = 1) Finite Elements. Default is \code{order} = 1.}
#' @examples
#' #Load the matrix nodes and tetrahedrons
#'
#' library(fdaPDE)
#' data(sphere3Ddata)
#'
#' nodes=sphere3Ddata$nodes
#' tetrahedrons=sphere3Ddata$tetrahedrons
#'
#' #Create the triangulated mesh from the connectivity matrix and nodes locations
#' mesh=create.MESH.3D(nodes,tetrahedrons)
#' 

create.MESH.3D<- function(nodes, tetrahedrons, order = 1)
{
  nnodes = dim(nodes)[1]

  ntetrahedrons = dim(tetrahedrons)[1]

  if(dim(tetrahedrons)[2]!= 4*order){
    if (order==1)
      stop("The matrix 'tetrahedrons' has the wrong number of columns. See second.order.mesh(...)")
  	stop("The matrix 'tetrahedrons' has wrong number of columns. Should be 4*order \n")
  	}
  out = list(nnodes=nnodes, ntetrahedrons=ntetrahedrons, nodes=c(t(nodes)), tetrahedrons = c(t(tetrahedrons)), order=as.integer(order))

  class(out)<-"MESH.3D"

  return(out)
}
