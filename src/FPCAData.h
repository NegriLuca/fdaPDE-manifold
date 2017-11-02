#ifndef __FPCADATA_HPP__
#define __FPCADATA_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"
#include "inputData.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  FPCAData: public InputData{
	private:
		//Design matrix
		MatrixXr datamatrix_;
		std::vector<UInt> observations_indices_;
		UInt n_;
		UInt p_;
		
		//Number of Principal Components
		UInt nPC_;
		
		#ifdef R_VERSION_
		void setDatamatrix(SEXP Rdatamatrix);
		#endif

	public:
		//! A basic version of the constructor.

		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Robservations an R-vector containing the values of the observations.
			\param Rdesmat an R-matrix containing the design matrix for the regression.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rbindex an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param Rbvalues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
		*/


		//! A complete version of the constructor.
		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Robservations an R-vector containing the values of the observations.
			\param Rdesmat an R-matrix containing the design matrix for the regression.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			\param Rbindex an R-integer vector containing the indexes of the nodes the user want to apply a Dirichlet Condition,
					the other are automatically considered in Neumann Condition.
			\param Rbvalues an R-double vector containing the value to impose for the Dirichlet condition, on the indexes specified in Rbindex
			\param Rc an R-double that contains the coefficient of the REACTION term
			\param Rbeta an R-double 2-dim vector that contains the coefficients for the TRANSPORT coefficients.
			\param RK an R-double 2X2 matrix containing the coefficients for a anisotropic DIFFUSION term.
			\param (UNSUPPORTED put it zero) Ru an R-double vector of length #triangles that contaiins the forcing term integrals.
		*/
		FPCAData(){};

		#ifdef R_VERSION_
		explicit FPCAData(SEXP Rlocations, SEXP Rdatamatrix, SEXP Rorder,
		SEXP Rlambda, SEXP RnPC, SEXP DOF);
		#endif

				
		explicit FPCAData(std::vector<Point>& locations, MatrixXr& datamatrix,
		UInt order, std::vector<Real> lambda, UInt nPC, bool DOF);


		void printDatamatrix(std::ostream & out) const;
		
		//void newDatamatrix(const VectorXr& scores_,const VectorXr& loadings_);
		
		//! A method returning a reference to the observations vector
		inline MatrixXr const & getDatamatrix() const {return datamatrix_;}
		
		//! A method returning the number of observations
		inline UInt const getNumberofObservations() const {return datamatrix_.cols();}
		inline std::vector<UInt> const & getObservationsIndices() const {return observations_indices_;}

		//! A method returning the number of Principal Components to compute
		inline UInt const getNPC() const {return nPC_;}
};

#include "FPCAData_imp.h"

#endif
