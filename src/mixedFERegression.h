#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include "mixedFE.h"

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegressionBase: public MixedFE<InputHandler, Integrator, ORDER, mydim, ndim>
{
protected:
	static constexpr Real dirichlet_penalization = 10e18;

	MatrixXr P_;

	void computeProjOnCovMatrix();

	//! A normal member taking two arguments: Dirichlet Boundary condition
	/*!
	 * This member applies Dirichlet boundary conditions on the linear system with the penalization method.
	  \param bcindex is a const reference to vector<int> : the global indexes of the nodes to which the boundary condition has to be applied.
	  \param bcvalues is a const reference to vector<double> : the values of the boundary conditions relative to bcindex.
	*/
	void applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues);
	void computeDataMatrix(SpMat& DMat);
	void computeDataMatrixByIndices(SpMat& DMat);
	void computeRightHandData(VectorXr& rightHandData);
	void computeDegreesOfFreedom(UInt output_index);

public:
	//!A Constructor.
	MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData):MixedFE<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData){};

	template<typename A>
	void apply(EOExpr<A> oper);
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, ndim, mydim>& mesh, const InputHandler& inputData):MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData){};

	void apply()
	{
		std::cout << "Option not implemented! \n";
	}
};



#include "mixedFERegression_imp.h"

#endif
