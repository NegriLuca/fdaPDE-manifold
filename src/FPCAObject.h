#ifndef __FPCAOBJECT_HPP__
#define __FPCAOBJECT_HPP__

#include "fdaPDE.h"
#include "mesh_objects.h"
#include "param_functors.h"

class  FPCAObject{
	private:
		
		
		//Loadings and scores estimation
		VectorXr scores_;
		VectorXr loadings_;
		
		//Regression data
		VectorXr ObservationData_;

	public:
		
		FPCAObject(){};
				
		explicit FPCAObject(const MatrixXr& datamatrix_);
		
		void setScores(const MatrixXr& datamatrix_);
		void setObservationData(const MatrixXr& datamatrix_);
		//void setObservationData(const MatrixXr& datamatrix_, const SpMat& psi_);
		void setLoadingsPsi(UInt nnodes, const VectorXr& f_sol,const SpMat& psi);
		void setLoadings(UInt nnodes, const VectorXr& f_sol, const std::vector<UInt>& obs_indices);
		void finalizeLoadings(const std::vector<UInt>& obs_indices, UInt nlocations);
		
		

		//!A method for printing the Scores vector on a specified output
		void printScores(std::ostream & out) const;
		//!A method for printing the Loadings vector on a specified output
		void printLoadings(std::ostream & out) const;
		//!A method for printing the ObservationData on a specified output
		void printObservationData(std::ostream & out) const;
		//void printDatamatrix(std::ostream & out) const;
		
			
		//! A method returning a reference to the scores vector
		inline VectorXr const & getScores() const {return scores_;}
		//! A method returning a reference to the loadings vector
		inline VectorXr const & getLoadings() const {return loadings_;}
		//! A method returning a reference to the observation data vector
		inline VectorXr const & getObservationData() const {return ObservationData_;}
};

#include "FPCAObject_imp.h"

#endif
