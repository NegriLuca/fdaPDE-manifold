## fdaPDE Library

fdaPDE implements a class of spatial regression models in between statistics and numerical analysis. These models are particularly well-suited for situations where a prior knowledge of the underlying phenomenon is known and can be described in terms of a Partial Differential Equation (PDE).

## Installation
Subfolder structure
src contains all C++ code and a special file named Makevars necessary to build and install the R package
R contains the R functions that wrap the C++ calls
data contains all .rda and .RData files useful for testings
tests contains basic R script to run tests
To install the package, please make sure that you have the package devtools already intalled. If using a Linux machine, it is also advisable to install rgl, plot3D and plot3Drgl before fdaPDE

From the root folder then type

R -e "library(devtools); install()" --silent
R -e "library(devtools); document()" --silent

To install the package from the Github repository, use the command `install_github` from the 
R package `devtools` as follows:
`install_github("NegriLuca/fdaPDE-manifold")`

## Examples

Some example can be found visualizing the help of the three main smoothing functions i.e.
`?smooth.FEM.basis` or `?smooth.FEM.PDE.basis` or `?smooth.FEM.PDE.SV.basis`

