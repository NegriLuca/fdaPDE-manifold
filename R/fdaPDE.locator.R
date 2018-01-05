#' Evaluate a FEM object at a set of point locations
#' 
#' @param FEM A \code{FEM} object to be evaluated.
#' @param locations A 2-colums(in case of planar mesh) or 3-columns(in case of 2D manifold in a 3D space) matrix with the spatial locations where the FEM object should be evaluated.
#' @param CPP_CODE Boolean. If \code{TRUE} the computation relies on the C++ implementation of a Visibility Walk Algorithm (Devillers et al. 2001). This usually ensures a fast computation. In case of the 2D manifold in a 3D space, only the C++ method is available.
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
   	
