library(fdaPDE)

data(hub)

cat('Plotting the mesh \n')
plot(hub)

## Generate some random data ##

nnodes = hub$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*hub$nodes[3*i+1]) +  a2* sin(2*pi*hub$nodes[3*i+2]) +  a3*sin(2*pi*hub$nodes[3*i+3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
FEMbasis <- create.FEM.basis(hub)

lambda=c(0.00375)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda,
                             CPP_CODE = TRUE)

cat("Showing result")
plot(output_CPP$fit.FEM)
