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
		
		//Number of Principal Components
		UInt nPC_;
		
		//Number of Folds for KFold
		UInt nFolds_;
		
		//Parameters for better GCV timings
		UInt GCVmethod_;
		UInt nrealizations_;      // Number of relizations for the stochastic estimation of GCV
		
		#ifdef R_VERSION_
		void setDatamatrix(SEXP Rdatamatrix);
		void setLocations(SEXP Rlocations);
		void setNrealizations(SEXP Rnrealizations);
		#endif

	public:
		//! A basic version of the constructor.

		/*!
			It initializes the object storing the R given objects. This is the simplest of the two possible interfaces with R
			\param Rdatamatrix an R-matrix containing the datamatrix.
			
			\param Rlocations an R-matrix containing the location of the observations.
			
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
			
			\param RnPC an R-integer specifying the number of principal components to compute.
			
			\param RnFolds an R-integer specifying the number of folds to use if K-Fold cross validation method is chosen.
			
			\param RGCVmethod an R-integer specifying if the GCV computation has to be exact(if = 1) or stochastic (if = 2).
			
			\param Rnrealizations an R-integer specifying the number of realizations to use when computing the GCV stochastically.
		
		*/

		FPCAData(){};

		#ifdef R_VERSION_
		explicit FPCAData(SEXP Rlocations, SEXP Rdatamatrix, SEXP Rorder,
		SEXP Rlambda, SEXP RnPC, SEXP RnFolds,SEXP RGCVmethod, SEXP Rnrealizations);
		#endif

				
		explicit FPCAData(std::vector<Point>& locations, MatrixXr& datamatrix,
		UInt order, std::vector<Real> lambda, UInt nPC, UInt nFolds);


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
		
		//! A method returning the the penalization term
		inline std::vector<Real> const & getLambda() const {return lambda_;}
		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		//! A method returning the input order
		inline UInt const getNFolds() const {return nFolds_;}
		//! A method returning the method that should be used to compute the GCV:
		//! 1: exact calculation
		//! 2: stochastic estimation
		inline UInt const & getGCVmethod() const {return GCVmethod_;}
		//! A method returning the number of vectors to use to stochastically estimate the edf
		inline UInt const & getNrealizations() const {return nrealizations_;}
		
		
};

#include "FPCAData_imp.h"

#endif
