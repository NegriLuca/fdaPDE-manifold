#ifndef __MIXEDFEFPCA_HPP__
#define __MIXEDFEFPCA_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEfpcaBase{
protected:
	static constexpr Real pruning_coeff = 2.2204e-013;
	static constexpr Real dirichlet_penalization = 10e18;
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const FPCAData& fPCAData_;
	std::vector<coeff> tripletsData_;

	//SpMat NWblock_;
	SpMat DMat_;
	SpMat AMat_;
	SpMat MMat_;

	SpMat Psi_;
	MatrixXr P_;

	SpMat coeffmatrix_;       //!A Eigen::VectorXr: Stores the system right hand side.
	VectorXr b_;			  //!A Eigen::VectorXr : Stores the system solution
	std::vector<VectorXr> solution_;
	std::vector<Real> dof_;

	void computeBasisEvaluations();

	//	|DMat | AMat^T  |
	//	|AMat | MMat	|
	void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
	void computeDataMatrix(SpMat& DMat);
	void computeDataMatrixByIndices(SpMat& DMat);
	void computeRightHandData(VectorXr& rightHandData);
	void computeDegreesOfFreedom(UInt output_index);
	 //! A template for the system resolution: SpLu, SpQR, SpCholesky,SpConjGrad
	template<typename P>
	void solve(UInt output_index);

public:
	//!A Constructor.
	MixedFEfpcaBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const FPCAData& fPCAData):mesh_(mesh), fPCAData_(fPCAData){};

	template<typename A>
	void apply(EOExpr<A> oper);

	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return solution_;};
	inline std::vector<Real> const & getDOF() const{return dof_;};
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEfpca : public MixedFEfpcaBase<Integrator, ORDER, mydim, ndim>
{
public:
	MixedFEfpca(const MeshHandler<ORDER, ndim, mydim>& mesh, const FPCAData& fPCAData):MixedFEfpcaBase<Integrator, ORDER, mydim, ndim>(mesh, fPCAData){};

	void apply()
	{
		std::cout << "Option not implemented! \n";
	}
};



#include "mixedFEfpca_imp.h"

#endif
