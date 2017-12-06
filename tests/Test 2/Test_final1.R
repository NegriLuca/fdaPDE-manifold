rm(list= ls())
#library(Rmpi)
library(fdaPDE2)

source("mesh.R")

################################################################################
# Edit here

# A vector containing the Ns of the grids to be used
N = c(30)
# The number of observations to be generated (same length as N)
n_observations = c(30)
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
lambdat = seq(-3,1,0.25) 
lambda = 10^lambdat
# The order of FEM
order = 1
# Numbers of realizations
#nreal = 30 
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

######### COVARIATES, LOC NOT ON NODES 
# cat("\nCOVARIATES, LOC NOT ON NODES\n\n")
# output_CPP_exact1 = smooth.FEM.basis(  
#                                 observations = observations[[i]],
#                                 locations=locations[[i]],
#                                 FEMbasis = FEMbasis[[i]],
#                                 lambda = lambda,
#                                 covariates=covariates[[i]],
#                                 GCV = TRUE,
#                                 CPP_CODE = TRUE,
#                                 GCVmethod=1 )
# edf_exact1 = output_CPP_exact1$edf
# gcv_exact1 = output_CPP_exact1$GCV

# output_CPP_stochastic1 = smooth.FEM.basis(  
#                                 observations = observations[[i]],
#                                 locations=locations[[i]],
#                                 FEMbasis = FEMbasis[[i]],
#                                 lambda = lambda,
#                                 covariates=covariates[[i]],
#                                 GCV = TRUE,
#                                 CPP_CODE = TRUE,
#                                 GCVmethod=2,
#                                 nrealizations = 1000 )
# edf_stochastic1 = output_CPP_stochastic1$edf
# gcv_stochastic1 = output_CPP_stochastic1$GCV

# pdf('test11_edf')
# plot(lambda, edf_exact1,ylim = range(c(edf_exact1,edf_stochastic1)),lwd = 2,type="s",col="blue",ylab="edf") 
# points(lambda, edf_stochastic1, lwd = 2, type="s",col="red")
# legend ("topright",legend = c ("Exact","Stochastic"),col = c("blue","red"), lty=c(1,1),lwd=c(3,3))
# dev.off()

# pdf('test11_gcv')
# plot(lambda, gcv_exact1,ylim = range(c(gcv_exact1,gcv_stochastic1)),lwd = 2, type="s",col="blue",ylab="gcv")
# points(lambda, gcv_stochastic1,lwd = 2, type="s",col="red")
# legend ("bottomright",legend = c ("Exact","Stochastic"),col = c("blue","red"), lty=c(1,1),lwd=c(3,3))
# dev.off()


############ NO COVARIATES, LOC ON NODES
cat("\nNO COVARIATES, LOC ON NODES\n\n")
output_CPP_exact2 =smooth.FEM.basis(observations = observations[[i]],
                             FEMbasis = FEMbasis[[i]],
                             lambda = lambda,
                             GCV = TRUE,
                             CPP_CODE = TRUE,
                             GCVmethod = 1)
edf_exact2 = output_CPP_exact2$edf
gcv_exact2 = output_CPP_exact2$GCV

output_CPP_stochastic2 =smooth.FEM.basis(observations = observations[[i]],
                             FEMbasis = FEMbasis[[i]],
                             lambda = lambda,
                             GCV = TRUE,
                             CPP_CODE = TRUE,
                             nrealizations = 1000,
                             GCVmethod = 2)
edf_stochastic2 = output_CPP_stochastic2$edf
gcv_stochastic2 = output_CPP_stochastic2$GCV

pdf('test12_edf')
plot(lambda, edf_exact2,ylim = range(c(edf_exact2,edf_stochastic2)),lwd = 2,type="s",col="blue",ylab="edf") 
points(lambda, edf_stochastic2,lwd = 2,type="s",col="red")
legend ("topright",legend = c ("Exact","Stochastic"),col = c("blue","red"), lty=c(1,1),lwd=c(3,3))
dev.off()

pdf('test12_gcv')
plot(lambda, gcv_exact2,ylim = range(c(gcv_exact2,gcv_stochastic2)),lwd = 2, type="s",col="blue", ylab="gcv")
points(lambda, gcv_stochastic2,lwd = 2,type="s",col="red")
legend ("topright",legend = c ("Exact","Stochastic"),col = c("blue","red"), lty=c(1,1),lwd=c(3,3))
dev.off()