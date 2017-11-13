#ifndef __MIXEDFEFPCA_HPP__
#define __MIXEDFEFPCA_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "solver.h"
#include "FPCAData.h"
#include "mixedFE.h"
#include "FPCAObject.h"

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCABase: public MixedFE<InputHandler, Integrator, ORDER, mydim, ndim>
{
protected:
	std::vector<VectorXr> scores_mat_;
	std::vector<VectorXr> loadings_mat_;
	std::vector<Real> lambda_PC_;
	std::vector<Real> variance_explained_;
	std::vector<Real> cumsum_percentage_;

	void computeDataMatrix(SpMat& DMat);
	void computeDataMatrixByIndices(SpMat& DMat);
	void computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput);
	void computeDegreesOfFreedom(UInt output_index);
	void computeVarianceExplained();
	void computeCumulativePercentageExplained();
	
	

public:
	//!A Constructor.
	MixedFEFPCABase(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData):MixedFE<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData){};

	template<typename A>
	void apply(EOExpr<A> oper);
	
	//! A method returning a reference to the scores matrix
	inline std::vector<VectorXr> const & getScoresMat() const {return scores_mat_;}
	//! A method returning a reference to the loadings matrix
	inline std::vector<VectorXr> const & getLoadingsMat() const {return loadings_mat_;}
	//! A method returning a reference to the vector of lambda taken for each PC
	inline std::vector<Real> const & getLambdaPC() const {return lambda_PC_;}
	//! A method returning a reference to the vector of variance explained for each PC
	inline std::vector<Real> const & getVarianceExplained() const {return variance_explained_;}
	//! A method returning a reference to the vector of the percentage explained cumulatively by the first N PC
	inline std::vector<Real> const & getCumulativePercentage() const {return cumsum_percentage_;}
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCA : public MixedFEFPCABase<InputHandler, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFEFPCA(const MeshHandler<ORDER, ndim, mydim>& mesh, const InputHandler& inputData):MixedFEFPCABase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData){};

	void apply()
	{
		std::cout << "Option not implemented! \n";
	}
};



#include "mixedFEFPCA_imp.h"

#endif
