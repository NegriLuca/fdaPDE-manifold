library(fdaPDE)

data(hub)

#cat('Plotting the mesh \n')
#plot.MESH2.5D(hub)

## Generate some random data ##

nnodes = hub$nnodes
cat("Nnodes: \n")
nnodes
cat("\n")

datamatrix<-NULL;
for(ii in 1:50){
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*hub$nodes[3*i+1]) +  a2* sin(2*pi*hub$nodes[3*i+2]) +  a3*sin(2*pi*hub$nodes[3*i+3]) +1
}

data=func_evaluation+rnorm(nnodes,mean=0,sd=0.5)
datamatrix<-rbind(datamatrix,data)
}

#datamatrix

data_bar=colMeans(datamatrix)
data_demean=matrix(rep(data_bar,50),nrow=50,byrow=TRUE)

#data_demean

datamatrix_demeaned=datamatrix-data_demean



FEMbasis <- create.FEM.basis(hub)

lambda=c(0.00375)
output_CPP =smooth.FEM.FPCA(datamatrix = datamatrix_demeaned,
                             FEMbasis = FEMbasis, lambda = lambda,nPC=2)
                             
str(output_CPP)


plot(output_CPP$loadings.FEM)


#cat("Showing result")
#plot.MESH2.5D(hub,output_CPP$fit.FEM$coeff)

message("Press Return to continue")
invisible(readLines("stdin",n=1))
