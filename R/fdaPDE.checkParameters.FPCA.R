checkSmoothingParametersFPCA<-function(locations = NULL, datamatrix, FEMbasis, lambda,nPC, validation, NFolds,GCVmethod = 2,nrealizations = 100,RNGstate = "",solver = "EigenLU",nprocessors = 1)
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
  if(is.null(nPC))
    stop("nPC required; is NULL.")
  if(!is.null(validation)){
    if(validation!="GCV" && validation!="KFold")
   	stop("'validation' needs to be 'GCV' or 'KFold'")
    if(validation=="KFold" && is.null(NFolds))
   	stop("NFolds is required if 'validation' is 'KFold'")
   }
   if (GCVmethod != 1 && GCVmethod != 2)
    stop("GCVmethod must be either 1(exact calculation) or 2(stochastic estimation)")

  if( !is.numeric(nrealizations) || nrealizations < 1)
    stop("nrealizations must be a positive integer")

  if ( !is.character(RNGstate))
    stop("Invalid RNGstate: it needs to be a character previously returned by a smoothing function")
  
  if (solver != "MumpsSparse" && solver != "EigenSparseLU")
    stop("Solver must be either \"MumpsSparse\" or \"EigenSparseLU\"")

  if (!is.numeric(nprocessors) || nprocessors < 1)
    stop("nprocessors must be a positive integer")

}

checkSmoothingParametersSizeFPCA<-function(locations = NULL, datamatrix, FEMbasis, lambda, ndim, mydim, validation, NFolds)
{
  #################### Parameter Check #########################
  if(nrow(datamatrix) < 1)
    stop("'observations' must contain at least one element")
  if(is.null(locations))
  {
    if(class(FEMbasis$mesh) == "MESH2D"){
    	if(ncol(datamatrix) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }else if(class(FEMbasis$mesh) == "MESH.2.5D"){
    	if(ncol(datamatrix) > FEMbasis$mesh$nnodes)
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }
  }
  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' and the mesh points have incompatible size;")
    if(nrow(locations) != ncol(datamatrix))
      stop("'locations' and 'observations' have incompatible size;")
  }
  if(ncol(lambda) != 1)
    stop("'lambda' must be a column vector")
  if(nrow(lambda) < 1)
    stop("'lambda' must contain at least one element")
  if(nrow(lambda)>1 && is.null(validation))
    stop("If 'lambda' contains more than one element, 'validation' needs to be specified as 'GCV' or 'KFold'")
}
