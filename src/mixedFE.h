#ifndef __MIXEDFE_HPP__
#define __MIXEDFE_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFE{
protected:
	static constexpr Real pruning_coeff = 2.2204e-013;
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const InputHandler& inputData_;
	std::vector<coeff> tripletsData_;

	//SpMat NWblock_;
	SpMat DMat_;
	SpMat AMat_;
	SpMat MMat_;

	SpMat Psi_;

	SpMat coeffmatrix_;       //!A Eigen::VectorXr: Stores the system right hand side.
	VectorXr b_;			  //!A Eigen::VectorXr : Stores the system solution
	std::vector<VectorXr> solution_;
	std::vector<Real> dof_;

	void computeBasisEvaluations();

	void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
	
	template<typename P>
	void solve(UInt output_index);


public:
	//!A Constructor.
	MixedFE(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData):mesh_(mesh), inputData_(inputData){};


	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return solution_;};
	inline std::vector<Real> const & getDOF() const{return dof_;};
};

#include "mixedFE_imp.h"

#endif
