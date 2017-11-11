#ifndef __INPUTDATA_HPP__
#define __INPUTDATA_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"

//!  An IO handler class for objects passed from R
/*!
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/

class InputData{
	private:
		std::vector<Point> locations_;

		bool locations_by_nodes_;


		//Other parameters
		UInt order_;
		std::vector<Real> lambda_;

		//bool inputType;
		bool DOF_;

		#ifdef R_VERSION_
		void setLocations(SEXP Rlocations);
		#endif

	public:
	InputData(){};

		#ifdef R_VERSION_
		explicit InputData(SEXP Rlocations, SEXP Rorder, SEXP Rlambda, SEXP DOF);
		#endif

		explicit InputData(std::vector<Point>& locations, UInt order, std::vector<Real> lambda, bool DOF);


		void printLocations(std::ostream & out) const;

		//! A method returning the locations of the observations
		inline std::vector<Point> const & getLocations() const {return locations_;}
		inline bool isLocationsByNodes() const {return locations_by_nodes_;}
		inline bool computeDOF() const {return DOF_;}
		//! A method returning the the penalization term
		inline std::vector<Real> const & getLambda() const {return lambda_;}
		//! A method returning the input order
		inline UInt const getOrder() const {return order_;}
		//inline MatrixXr const & getCovariates() const {return {};}
		//inline std::vector<UInt> const & getDirichletIndices() const {return {};}
		//inline std::vector<Real> const & getDirichletValues() const {return {};}
		//inline VectorXr const & getObservations() const {return {};}
};

#include "inputData_imp.h"

#endif
