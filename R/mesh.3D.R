#' Create a \code{MESH.3D} object from the connectivty matrix and nodes locations
#'
#' @param nodes A #nodes-by-3 matrix specifying the locations of each node
#' @param tetrahedrons A #tetrahedrons-by-4*order matrix specifying the indices of the nodes in each tetrahedrons
#' @param order Either "1" or "2". Order of the Finite Element basis default is order = 1
#' @return An object of the class \code{MESH.3D} with the following output:
#' \item{\code{nnodes}}{The #nodes contained in the mesh}
#' \item{\code{ntetrahedrons}}{The #tetrahedrons contained in the mesh}
#' \item{\code{nodes}}{A vector of length #nodes*3 containing the x,y and z coordinate for each point of the mesh}
#' \item{\code{tetrahedrons}}{A vector of length #tetrahedrons*4*order specifying the indices of the nodes in each triangle of the mesh}
#' \item{\code{order}}{Either '1' or '2'. It specifies wether each mesh triangle should be represented by 3 nodes (the triangle' vertices) or by 6 nodes (the triangle's vertices and midpoints). 
#' These are respectively used for linear (order = 1) and quadratic (order = 2) Finite Elements. Default is \code{order} = 1.}
#' @examples
#' #Load the matrix nodes and tetrahedrons
#'
#' library(fdaPDE)
#' data(sphereData)
#'
#' nodes=sphere$nodes
#' triangles=sphere$triangles
#'
#' #Create the triangulated mesh from the connectivity matrix and nodes locations
#' mesh=create.MESH.2.5D(nodes,triangles)
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
