library(fdaPDE)

data(candy)

cat('Plotting the mesh \n')
plot(candy)

## Generate some random data ##

nnodes = candy$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*candy$nodes[3*i+1]) +  a2* sin(2*pi*candy$nodes[3*i+2]) +  a3*sin(2*pi*candy$nodes[3*i+3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
FEMbasis <- create.FEM.basis(candy)

lambda=c(0.1,0.2,0.3,0.4)
output_CPP =smooth.FEM.basis(observations = data,
                             FEMbasis = FEMbasis, lambda = lambda,
                             CPP_CODE = TRUE, GCV=TRUE)

cat("Showing result for different values of lambda")
plot(output_CPP$fit.FEM)


lambda_opt = lambda[which(output_CPP$GCV == min(output_CPP$GCV))]
cat("From the GCV analysis, the best lambda is ", lambda_opt, "\n")
