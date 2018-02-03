#' Smooth Functional Principal Component Analysis
#' 
#' @param datamatrix A matrix of dimensions #samples-by-#points(or #observations) with the observed data values over the domain for each sample. 
#' The locations of the observations can be specified with the \code{locations} argument.
#' The datamatrix needs to have zero mean.
#' Otherwise if only the datamatrix is given, these are considered to be located in the corresponding node in the table
#' \code{nodes} of the mesh. In this last case, a column with \code{NA} value in the \code{datamatrix} indicates that there is no observation associated to the corresponding
#'  node.
#' @param locations A #observations-by-ndim matrix where each row specifies the spatial coordinates \code{x} and \code{y} (and \code{z} if ndim=3) of the corresponding column of observations in the matrix \code{datamatrix}.
#' This parameter can be \code{NULL}. In this case the spatial coordinates of the corresponding observations are assigned as specified in \code{datamatrix}.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis, as created by \code{\link{create.FEM.basis}}.
#' @param lambda A scalar or vector of smoothing parameters.
#' @param nPC An integer specifying the number of Principal Components to compute.
#' @param validation String. If \code{lambda} is a vector, it has to be specified as \code{"GCV"} or \code{"KFold"}. This parameter specify which method of cross-validation is used to select the best parameter \code{lambda} among those values of the smoothing parameter specified in \code{lambda} for each Principal Component.
#' @param NFolds An integer specifying the number of folds to use if the KFold cross-validation method for the selection of the best parameter \code{lambda} is chosen. Default value is 5.
#' @param GCVmethod If set to 1 perform an exact (but possibly slow) computation of the GCV index. If set to 2 approximate the GCV with a stochastic algorithm.
#' @param nrealizations The number of realizations to be used in the stochastic algorithm for the estimation of GCV.
#' @return A list with the following variables:
#' \item{\code{loadings.FEM}}{A \code{FEM} object that represents the normalized functional loadings for each Principal Component computed.}
#' \item{\code{scores}}{A #samples-by-#PrincipalComponents matrix that represents the unnormalized scores or PC vectors.}
#' \item{\code{lambda}}{A vector of length #PrincipalComponents where each value is \code{lambda} chosen for that Principal Component.}
#' \item{\code{variance_explained}}{A vector of length #PrincipalComponents where each value represent the variance explained by that component.}
#' \item{\code{ cumsum_percentage}}{A vector of length #PrincipalComponents containing the cumulative percentage of the variance explained by the first components.}
#' \item{\code{var}}{If GCV is \code{TRUE} and GCVmethod = 2, a scalar or vector with the sample variance of the realizations of the stochastic edf estimator for each value of the smoothing parameter specified in \code{lambda}.}
#' @description This function implements a smooth functional principal component analysis over a planar mesh, a smooth manifold or a volume. In order to perform some regularization, the calculation involves the Laplacian of the spatial field. The computation relies only on the C++ implementation of the algorithm.
#' @usage smooth.FEM.FPCA(locations = NULL, datamatrix, FEMbasis, lambda, nPC=1, validation=NULL,NFolds=5,GCVmethod = 2, nrealizations = 100)
#' @references Sangalli, L.M., Ramsay, J.O. & Ramsay, T.O., 2013. Spatial spline regression models. Journal of the Royal Statistical Society. Series B: Statistical Methodology, 75(4), pp. 681-703.
#' @examples
#' library(fdaPDE)
#' ## Load the hub data
#' data(hub)
#' ## Plot the mesh
#' plot(hub)
#' ## Create the Finite Element basis 
#' FEMbasis = create.FEM.basis(hub)
#' ## Create a datamatrix
#'datamatrix<-NULL;
#'for(ii in 1:50){
#'a1 = rnorm(1,mean = 1, sd = 1)
#'a2 = rnorm(1,mean = 1, sd = 1)
#'a3 = rnorm(1,mean = 1, sd = 1)
#'
#'func_evaluation = numeric(nnodes)

#'for (i in 0:(nnodes-1)){
#'  func_evaluation[i+1] = a1* sin(2*pi*hub$nodes[3*i+1]) +  a2* sin(2*pi*hub$nodes[3*i+2]) +  a3*sin(2*pi*hub$nodes[3*i+3]) +1
#'}

#'data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
#'datamatrix<-rbind(datamatrix,data)
#'}
#' ## Compute the mean of the datamatrix and subtract it
#data_bar=colMeans(datamatrix)
#data_demean=matrix(rep(data_bar,50),nrow=50,byrow=TRUE)
#'
#'datamatrix_demeaned=datamatrix-data_demean
#' ##Fix the parameter lambda
#'lambda=c(0.00375)
#' ## Estimate the first 2 Principal Components
#'output_CPP =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
#'                            FEMbasis = FEMbasis, lambda = lambda,nPC=2)
#'
#' ## Plot the functional loadings of the estimated Principal Components                           
#'plot(output_CPP$loadings.FEM)


#observation is a matrix with nrow=number of location points (or mesh nodes) and n cols=n of observations

smooth.FEM.FPCA<-function(locations = NULL, datamatrix, FEMbasis, lambda, nPC=1, validation=NULL,NFolds=5,GCVmethod = 2, nrealizations = 100)
{
 if(class(FEMbasis$mesh) == "MESH2D"){
 	ndim = 2
 	mydim = 2
 }else if(class(FEMbasis$mesh) == "MESH.2.5D"){
 	ndim = 3
 	mydim = 2
 }else if(class(FEMbasis$mesh) == "MESH.3D"){
 	ndim = 3
 	mydim = 3
 }else{
 	stop('Unknown mesh class')
 }
 
##################### Checking parameters, sizes and conversion ################################

  checkSmoothingParametersFPCA(locations, datamatrix, FEMbasis, lambda,nPC, validation,NFolds,GCVmethod , nrealizations) 
  ## Coverting to format for internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  datamatrix = as.matrix(datamatrix)
  lambda = as.matrix(lambda)
  
  checkSmoothingParametersSizeFPCA(locations, datamatrix, FEMbasis, lambda, ndim, mydim, validation, NFolds)
  
	  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'MESH2D'){	  
  	bigsol = NULL
	print('C++ Code Execution')
	bigsol = CPP_smooth.FEM.FPCA(locations, datamatrix, FEMbasis, lambda,
	ndim, mydim,nPC, validation, NFolds,GCVmethod, nrealizations)
	numnodes = nrow(FEMbasis$mesh$nodes)
	  
  } else if(class(FEMbasis$mesh) == 'MESH.2.5D'){

	  bigsol = NULL  
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.manifold.FEM.FPCA(locations, datamatrix, FEMbasis$mesh,
	  lambda, ndim, mydim,nPC, validation, NFolds,GCVmethod, nrealizations)
	  numnodes = FEMbasis$mesh$nnodes
  } else if(class(FEMbasis$mesh) == 'MESH.3D'){

	  bigsol = NULL  
	  print('C++ Code Execution')
	  bigsol = CPP_smooth.volume.FEM.FPCA(locations, datamatrix, FEMbasis$mesh,
	  lambda, ndim, mydim,nPC, validation, NFolds,GCVmethod, nrealizations)
	  numnodes = FEMbasis$mesh$nnodes
  }
  
  loadings=bigsol[[1]]
  loadings.FEM=FEM(loadings,FEMbasis)
  
  scores=bigsol[[2]]
  
  lambda=bigsol[[3]]
  
  variance_explained=bigsol[[4]]
  
  cumsum_percentage=bigsol[[5]]
  
  var=bigsol[[6]]
  
  reslist=list(loadings.FEM=loadings.FEM, scores=scores, lambda=lambda, variance_explained=variance_explained, cumsum_percentage=cumsum_percentage, var=var)
  return(reslist)
  }
