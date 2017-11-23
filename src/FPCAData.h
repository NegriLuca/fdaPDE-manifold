#ifndef __FPCADATA_HPP__
#define __FPCADATA_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
class  FPCAData{
	private:
	
		std::vector<Point> locations_;

		bool locations_by_nodes_;
	
		//Design matrix
		MatrixXr datamatrix_;
		std::vector<UInt> observations_indices_;
		UInt n_;
		UInt p_;
		
		//Other parameters
		UInt order_;
		std::vector<Real> lambda_;

		//bool inputType;
		bool DOF_;
		
		//bool KF
		bool KFold_;
		
		//Number of Principal Components
		UInt nPC_;
		
		//Number of Folds for KFold
		UInt nFolds_;
		
		#ifdef R_VERSION_
		void setDatamatrix(SEXP Rdatamatrix);
		void setLocations(SEXP Rlocations);
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
		SEXP Rlambda, SEXP RnPC, SEXP Rvalidation, SEXP RnFolds);
		#endif

				
		explicit FPCAData(std::vector<Point>& locations, MatrixXr& datamatrix,
		UInt order, std::vector<Real> lambda, UInt nPC, std::string validation, UInt nFolds);


		void printDatamatrix(std::ostream & out) const;
		void printLocations(std::ostream & out) const;

		//! A method returning the locations of the observations
		inline std::vector<Point> const & getLocations() const {return locations_;}
		//! A method returning TRUE if the observations are located in the nodes of the mesh or FALSE otherwise
		inline bool isLocationsByNodes() const {return locations_by_nodes_;}

		//void newDatamatrix(const VectorXr& scores_,const VectorXr& loadings_);
		
		//! A method returning a reference to the observations vector
		inline MatrixXr const & getDatamatrix() const {return datamatrix_;}
		
		//! A method returning the number of observations
		inline UInt const getNumberofObservations() const {return datamatrix_.cols();}
		//! A method returning the observations indices
		inline std::vector<UInt> const & getObservationsIndices() const {return observations_indices_;}

		//! A method returning the number of Principal Components to compute
		inline UInt const getNPC() const {return nPC_;}
		
		//! A method returning a boolean value specifying if the Degrees of Freedom needs to be computed
		inline bool computeDOF() const {return DOF_;}
		//! A method returning the the penalization term
		inline std::vector<Real> const & getLambda() const {return lambda_;}
		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		
};

#include "FPCAData_imp.h"

#endif
