#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include <memory>

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegressionBase
{
protected:
	static constexpr Real dirichlet_penalization = 10e18;
	static constexpr Real pruning_coeff = 2.2204e-013;
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const InputHandler& regressionData_;
	std::vector<coeff> tripletsData_;
	
	SpMat A_;	// System matrix with psi^T*psi in north-west block
	SpMat R1_;	// North-east block of system matrix A_
	SpMat R0_;	// South-east block of system matrix A_
	SpMat psi_;
	MatrixXr U_;	// psi^T*W padded with zeros
		
		
	Eigen::SparseLU<SpMat> Adec_; // Stores the factorization of A_
	Eigen::PartialPivLU<MatrixXr> Gdec_;	// Stores factorization of G =  C + [V * A^-1 * U]
	Eigen::PartialPivLU<MatrixXr> WTWinv_;	// Stores the factorization of W^T * W
	bool isWTWfactorized_;
	bool isRcomputed_;
	MatrixXr R_; //R1 ^T * R0^-1 * R1

	//SpMat NWblock_;
	SpMat DMat_;
	SpMat AMat_;
	SpMat MMat_;
	/*
	SpMat Psi_;
	MatrixXr P_;
	*/
	
	MatrixXr Q_;
 	MatrixXr H_;
	
	

	SpMat _coeffmatrix;        //!A Eigen::VectorXr: Stores the system right hand side.
	VectorXr _b;                     //!A Eigen::VectorXr: Stores the system right hand side.
	std::vector<VectorXr> _solution; //!A Eigen::VectorXr : Stores the system solution
	std::vector<Real> _dof;

	//void computeBasisEvaluations();
	//void computeProjOnCovMatrix();
	
	void setPsi();
	void buildA(const SpMat& Psi,  const SpMat& R1,  const SpMat& R0);
	MatrixXr LeftMultiplybyQ(const MatrixXr& u);
	void addDirichletBC();
	void setQ();
 	void setH();

	//	|DMat | AMat^T  |
	//	|AMat | MMat	|
	void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
	void getRightHandData(VectorXr& rightHandData);
	void computeDegreesOfFreedom(UInt output_index, Real lambda);
	void computeDegreesOfFreedomExact(UInt output_index, Real lambda);
	void computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);

	void system_factorize();
	template<typename Derived>
	MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);
	
	void getDataMatrix(SpMat& DMat);
 	void getDataMatrixByIndices(SpMat& DMat);
	
	//! A normal member taking two arguments: Dirichlet Boundary condition
	/*!
	 * This member applies Dirichlet boundary conditions on the linear system with the penalization method.
	  \param bcindex is a const reference to vector<int> : the global indexes of the nodes to which the boundary condition has to be applied.
	  \param bcvalues is a const reference to vector<double> : the values of the boundary conditions relative to bcindex.
	*/
	
	/*void applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues);
	void computeDataMatrix(SpMat& DMat);
	void computeDataMatrixByIndices(SpMat& DMat);
	void computeRightHandData(VectorXr& rightHandData);
	void computeDegreesOfFreedom(UInt output_index);
	
	//! A template for the system resolution: SpLu, SpQR, SpCholesky,SpConjGrad
	template<typename P>
	void solve(UInt output_index);
*/
public:
	//!A Constructor.
	MixedFERegressionBase(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& regressionData): mesh_(mesh),regressionData_(regressionData),isRcomputed_(false),isWTWfactorized_(false) {};

	template<typename A>
	void apply(EOExpr<A> oper);
	
	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return _solution;};
	inline std::vector<Real> const & getDOF() const{return _dof;};
};

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression : public MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, ndim, mydim>& mesh, const InputHandler& regressionData):MixedFERegressionBase<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		std::cout << "Option not implemented! \n";
	}
};



#include "mixedFERegression_imp.h"

#endif
