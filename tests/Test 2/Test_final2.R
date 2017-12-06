rm(list= ls())
library(Rmpi)
library(fdaPDE2)

source("mesh.R")
################################################################################
# Edit here

# A vector containing the Ns of the grids to be used
N = c(20)
# The number of observations to be generated (same length as N)
n_observations = c(20)
# The "true" coefficients of the covariates
beta = rbind(0.2, -0.4, 0.7, -0.05)
# Functions to be used to generate the covariates
f = vector("list", length(beta-1))
# Specify a function for each covariate -1, which is random
#f[[1]] <- function (x,y){x}
#f[[2]] <- function (x,y){y}
#f[[3]] <- function (x,y){x*y}
f[[1]] <- function (x,y){sin(2*pi*x*y)}
f[[2]] <- function (x,y){sin(2*pi*x)*sin(2*pi*y)}
f[[3]] <- function (x,y){sin(3*pi*x)*cos(4*pi*y)}
# The lambda to be used
lambda = c(1)
# The order of FEM
order = 1
# Numbers of realizations
nreal = c(100, 200, 400, 800, 1600)
################################################################################

n_meshes = length(N)
n_covariates = length(beta)

FEMbasis = vector("list", n_meshes)
locations = vector("list", n_meshes)
covariates = vector("list", n_meshes)
covariates_on_nodes = vector("list", n_meshes)
observations = vector("list", n_meshes)
observations_on_nodes = vector("list", n_meshes)
output = vector("list", n_meshes)
set.seed(7)
for (i in 1:n_meshes) {
    grid = mesh_quadratounitario(N[i])
    mesh = create.MESH.2D(nodes=grid$nodes, order = order)
    FEMbasis[[i]] = create.FEM.basis(mesh)
    locations[[i]] = cbind(cbind(runif(n_observations[i],0,1)),
                                cbind(runif(n_observations[i],0,1)))
    fun = runif(n_observations[i],-1,1)
    for (j in 1:(length(beta)-1)){
        fun = cbind(fun,f[[j]](locations[[i]][,1],locations[[i]][,2]))
    }
    covariates[[i]] = matrix(fun,
                             nrow = n_observations[i],
                             ncol = n_covariates)
    fun_on_nodes = runif(nrow(mesh$nodes),-1,1)
    for (j in 1:(length(beta)-1)){
        fun_on_nodes = cbind(fun_on_nodes,f[[j]](mesh$nodes[,1],mesh$nodes[,2]))
    }
    covariates_on_nodes[[i]] =
        matrix(fun_on_nodes,
               nrow = nrow(mesh$nodes),
               ncol = n_covariates)
    observations[[i]] = locations[[i]][,1]^2 * locations[[i]][,2]^2
                        + covariates[[i]] %*% beta
                        + rnorm(n = nrow(locations[[i]]), sd = 0.1)
    observations_on_nodes[[i]] = locations[[i]][,1]^2 * locations[[i]][,2]^2
                                 + covariates_on_nodes[[i]] %*% beta
                                 + rnorm(n = nrow(mesh$nodes), sd = 0.1)
    indeces_to_cut= sample(1:length(observations_on_nodes[[i]]), N[i] - n_observations[i], replace=F)
    observations_on_nodes[[i]][indeces_to_cut[]]=NaN
}
var_vector = vector("list", length(nreal) )
edf_vector_n = vector("integer", 100 )


######### COVARIATES, LOC NOT ON NODES 

cat("\nCOVARIATES, LOC NOT ON NODES\n\n")
for ( j in 1:length(nreal))
{
	output_CPP_stochastic1 = smooth.FEM.basis(  
	                                observations = observations[[i]],
	                                locations=locations[[i]],
	                                FEMbasis = FEMbasis[[i]],
	                                lambda = lambda,
	                                covariates=covariates[[i]],
	                                GCV = TRUE,
	                                CPP_CODE = TRUE,
	                                GCVmethod=2,
	                                nrealizations = nreal[j] ,
	                                RNGstate = "")
	RNG_state =  output_CPP_stochastic1$RNGstate
	for ( k in 1:200){
		output_CPP_stochastic1 = smooth.FEM.basis(  
                                observations = observations[[i]],
                                locations=locations[[i]],
                                FEMbasis = FEMbasis[[i]],
                                lambda = lambda,
                                covariates=covariates[[i]],
                                GCV = TRUE,
                                CPP_CODE = TRUE,
                                GCVmethod=2,
                                nrealizations = nreal[j],
                                RNGstate = RNG_state )
		edf_vector_n[k] = output_CPP_stochastic1$edf
		RNG_state =  output_CPP_stochastic1$RNGstate
	}
	var_vector[j]=var(edf_vector_n)
}
logvar <-NULL
for (i in 1:length(var_vector)){
	logvar=c(logvar,log(var_vector[[i]]))
}

pdf("2_testlog")
plot(log(nreal), logvar ,ylim = range(c(logvar,-log(nreal))),col="blue",lwd = 2,xlab= "log(nrealization)", ylab="log(var)") 
points (log(nreal), -log(nreal), lwd = 2,col="red")
lines(log(nreal), logvar ,col="blue",lwd = 2, )
lines(log(nreal), -log(nreal),lwd = 2,col="red")
legend ("topright",legend = c ("Computed","1/nreal"),col = c("blue","red"), lty=c(1,1),lwd=c(3,3))
dev.off()
