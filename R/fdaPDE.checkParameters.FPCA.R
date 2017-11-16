checkSmoothingParametersFPCA<-function(locations = NULL, datamatrix, FEMbasis, lambda, GCV = FALSE)
{
  #################### Parameter Check #########################
  if(!is.null(locations))
  {
    if(any(is.na(locations)))
      stop("Missing values not admitted in 'locations'.")
    if(any(is.na(datamatrix)))
      stop("Missing values not admitted in 'observations' when 'locations' are specified.")
  }
  if (is.null(datamatrix)) 
    stop("observations required;  is NULL.")
  if (is.null(FEMbasis)) 
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'")
  if (is.null(lambda)) 
    stop("lambda required;  is NULL.")

  if (is.null(GCV)) 
    stop("GCV required;  is NULL.")
  if(!is.logical(GCV))
    stop("'GCV' is not logical")

}

checkSmoothingParametersSizeFPCA<-function(locations = NULL, datamatrix, FEMbasis, lambda, GCV = FALSE, ndim, mydim)
{
  #################### Parameter Check #########################
  if(nrow(datamatrix) < 1)
    stop("'observations' must contain at least one element")
  if(is.null(locations))
  {
    if(class(FEMbasis$mesh) == "MESH2D"){
    	if(nrow(datamatrix) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }else if(class(FEMbasis$mesh) == "MESH.2.5D"){
    	if(nrow(datamatrix) > FEMbasis$mesh$nnodes)
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }
  }
  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' and the mesh points have incompatible size;")
    if(nrow(locations) != nrow(datamatrix))
      stop("'locations' and 'observations' have incompatible size;")
  }
  if(ncol(lambda) != 1)
    stop("'lambda' must be a column vector")
  if(nrow(lambda) < 1)
    stop("'lambda' must contain at least one element")
  
}
