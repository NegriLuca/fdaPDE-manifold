#' Evaluate a FEM object at a set of point locations
#' 
#' @param FEM A \code{FEM} object to be evaluated.
#' @param locations A 2-colums(in case of planar mesh) or 3-columns(in case of 2D manifold in a 3D space or a 3D volume) matrix with the spatial locations where the FEM object should be evaluated.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of a Visibility Walk Algorithm (Devillers et al. 2001). This usually ensures a fast computation. In case of the 2D manifold in a 3D space or a 3D volume, only the C++ method is available.
#' @return 
#' A matrix of numeric evaluations of the \code{FEM} object. Each row indicates the location where the evaluation has been taken, the column indicates the 
#' function evaluated.
#' @description It evaluates a FEM object the specified set of locations.  
#' @usage eval.FEM(FEM, locations, CPP_CODE = TRUE)
#' @references 
#'  Devillers, O. et al. 2001. Walking in a Triangulation, Proceedings of the Seventeenth Annual Symposium on Computational Geometry

eval.FEM <- function(FEM, locations, CPP_CODE = TRUE)
{
  if (is.null(FEM)) 
    stop("FEM required;  is NULL.")
  if(class(FEM) != "FEM")
    stop("'FEM' is not of class 'FEM'")
  if (is.null(locations)) 
    stop("locations required;  is NULL.")
  if (is.null(CPP_CODE)) 
    stop("CPP_CODE required;  is NULL.")
  if(!is.logical(CPP_CODE))
    stop("'CPP_CODE' is not logical")
  
  locations = as.matrix(locations)
  
  res = NULL
  
  if(class(FEM$FEMbasis$mesh)=='MESH2D'){
    ndim = 2
    mydim = 2
	  if(CPP_CODE == FALSE)
	  {
	    res = R_eval.FEM(FEM, locations)
	  }else{ 
	    res = CPP_eval.FEM(FEM, locations, TRUE, ndim, mydim)
	  }
  }else if(class(FEM$FEMbasis$mesh)=='MESH.2.5D'){
      ndim = 3
      mydim = 2
  	    res = CPP_eval.manifold.FEM(FEM, locations, TRUE, ndim, mydim)
  }else if(class(FEM$FEMbasis$mesh)=='MESH.3D'){
      ndim = 3
      mydim = 3
  	    res = CPP_eval.volume.FEM(FEM, locations, TRUE, ndim, mydim)
  	  }
  	  
  
  return(as.matrix(res))
}

#' Evaluate a FEM object at a grid of point locations
#' 
#' @param FEM A \code{FEM} object to be evaluated.
#' @param length_interval_x. An integer specifying the length of the interval on the x axis of the grid of points where the FEM object should be evaluated. Default value is 25.
#' @param length_interval_y. An integer specifying the length of the interval on the y axis of the grid of points where the FEM object should be evaluated. If no value is specified, the length is assumed to be equal to the \code{length_interval_x} parameter.
#' @param length_interval_z. An integer specifying the length of the interval on the z axis of the grid of points where the FEM object should be evaluated. If no value is specified, the length is assumed to be equal to the \code{length_interval_x} parameter.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of a Visibility Walk Algorithm (Devillers et al. 2001). This usually ensures a fast computation. In case of the 2D manifold in a 3D space or a 3D volume, only the C++ method is available.
#' @return 
#' An object of class \code{FEMevaluated} with the following output:
#' \item{\code{evaluated_values}}{A matrix of numeric evaluations of the \code{FEM} object. The element of the matrix corresponds to the locations of the elements of the meshgrid generated.}
#' \item{\code{meshgrid}}{The meshgrid generated using the specified intervals on the 3 axis.}
#' \item{\code{x}}{A vector of length \code{length_interval_x} containing the x coordinate for each point of the meshgrid generated.}
#' \item{\code{y}}{A vector of length \code{length_interval_y} containing the y coordinate for each point of the meshgrid generated.}
#' \item{\code{z}}{A vector of length \code{length_interval_z} containing the z coordinate for each point of the meshgrid generated.}
#' @description It evaluates a FEM object at a specified grid of point locations. It is needed for doing 3D plots like isosurface plot, slice plot and slicecontour plot.
#' @usage eval.FEM(FEM, length_interval_x=25, length_interval_y=NULL,length_interval_z=NULL, CPP_CODE=TRUE)
#' @seealso \code{\link{isosurfaces.FEMevaluated}},\code{\link{slices.FEMevaluated}} and \code{\link{slicecontours.FEMevaluated}}
#' @references 
#'  Devillers, O. et al. 2001. Walking in a Triangulation, Proceedings of the Seventeenth Annual Symposium on Computational Geometry
evaluate.grid.FEM <- function(FEM, length_interval_x=25, length_interval_y=NULL,length_interval_z=NULL, CPP_CODE=TRUE)
{
  if (is.null(FEM)) 
    stop("FEM required;  is NULL.")
  if(class(FEM) != "FEM")
    stop("'FEM' is not of class 'FEM'")
    if(is.null(length_interval_x))
    stop("length_interval_x required;is NULL.")
    if(is.null(length_interval_y))
    length_interval_y=length_interval_x
    if(is.null(length_interval_z))
     length_interval_z=length_interval_x
    if (is.null(CPP_CODE)) 
    stop("CPP_CODE required;  is NULL.")
  if(!is.logical(CPP_CODE))
    stop("'CPP_CODE' is not logical")
    
   nodes=matrix(FEM$FEMbasis$mesh$nodes,nrow=FEM$FEMbasis$mesh$nnodes,ncol=3,byrow=TRUE)
   x=seq(min(nodes[,1]),max(nodes[,1]),length.out=length_interval_x)
   y=seq(min(nodes[,2]),max(nodes[,2]),length.out=length_interval_y)
   z=seq(min(nodes[,3]),max(nodes[,3]),length.out=length_interval_z)
   
   M=mesh(x,y,z)
   
   loc=NULL
   for(i in 1:length_interval_x)
   {
   	for(j in 1:length_interval_y)
   	{
   		for(k in 1:length_interval_z)
   		{
   			loc=rbind(loc,c(M$x[i,j,k],M$y[i,j,k],M$z[i,j,k]))
   		}
   	}
   }
   
   
   eval=eval.FEM(FEM=FEM,locations=loc,CPP_CODE=TRUE)
   ndata=dim(eval)[2]
   
   evaluated_values=NULL
   
   for(idata in 1:ndata)
   {	name=paste("values",idata,sep="")
   	values=array(data=NA,dim=c(length_interval_x,length_interval_y,length_interval_z))
   	for(i in 1:length_interval_x)
   	{
   		for(j in 1:length_interval_y)
   		{
   			for(k in 1:length_interval_z)
   			{
   				values[i,j,k]=eval[length_interval_x^2*(i-1)+length_interval_y*(j-1)+k,idata]
   			}
   		}
   	}
   	evaluated_values[[idata]]=values
   	names(evaluated_values)[idata]=name
   }
   
   
   reslist=list(evaluated_values=evaluated_values,meshgrid=M, x=x, y=y, z=z)
   class(reslist)="FEMevaluated"
   return (reslist)
}
   	
