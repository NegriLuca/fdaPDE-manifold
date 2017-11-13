#observation is a matrix with nrow=number of location points (or mesh nodes) and n cols=n of observations

smooth.FEM.FPCA<-function(locations = NULL, datamatrix, FEMbasis, lambda, GCV = FALSE, CPP_CODE = TRUE, nPC=1)
{
 if(class(FEMbasis$mesh) == "MESH2D"){
 	ndim = 2
 	mydim = 2
 }else if(class(FEMbasis$mesh) == "MESH.2.5D"){
 	ndim = 3
 	mydim = 2
 }else{
 	stop('Unknown mesh class')
 }
 
##################### Checking parameters, sizes and conversion ################################

  checkSmoothingParametersFPCA(locations, datamatrix, FEMbasis, lambda,  GCV, CPP_CODE) 
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  datamatrix = as.matrix(datamatrix)
  lambda = as.matrix(lambda)
  
  checkSmoothingParametersSizeFPCA(locations, datamatrix, FEMbasis, lambda, GCV, CPP_CODE, ndim, mydim)
	  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'MESH2D'){	  
  	bigsol = NULL
	print('C++ Code Execution')
	bigsol = CPP_smooth.FEM.FPCA(locations, datamatrix, FEMbasis, lambda,
	ndim, mydim, GCV,nPC)
	numnodes = nrow(FEMbasis$mesh$nodes)
	  
  } else if(class(FEMbasis$mesh) == 'MESH.2.5D'){

	  bigsol = NULL  
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.manifold.FEM.FPCA(locations, datamatrix, FEMbasis$mesh,
	  lambda, ndim, mydim, GCV,nPC)
	  numnodes = FEMbasis$mesh$nnodes
  }
  
  loadings=bigsol[[1]]
  loadings.FEM=FEM(loadings,FEMbasis)
  
  scores=bigsol[[2]]
  
  lambda=bigsol[[3]]
  
  variance_explained=bigsol[[4]]
  
  cumsum_percentage=bigsol[[5]]
  
  reslist=list(loadings.FEM=loadings.FEM, scores=scores, lambda=lambda, variance_explained=variance_explained, cumsum_percentage=cumsum_percentage)
  return(reslist)
  }
