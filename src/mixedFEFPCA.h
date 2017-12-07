#ifndef __MIXEDFEFPCA_HPP__
#define __MIXEDFEFPCA_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "solver.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include <memory>

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCABase
{
protected:
	const MeshHandler<ORDER, mydim, ndim> &mesh_;
	const FPCAData& fpcaData_;
	std::vector<coeff> tripletsData_;
	
	SpMat R1_;	// North-east block of system matrix A_
	SpMat R0_;	// South-east block of system matrix A_
	SpMat psi_;

	bool isRcomputed_;
	MatrixXr R_; //R1 ^T * R0^-1 * R1

	//SpMat NWblock_;
	SpMat DMat_;
	SpMat AMat_;
	SpMat MMat_;

	SpMat Psi_;

	SpMat coeffmatrix_;       //!A Eigen::VectorXr: Stores the system right hand side.
	VectorXr b_;			  //!A Eigen::VectorXr : Stores the system solution
	std::vector<VectorXr> solution_;
	
	Sparse_LU sparseSolver_;
	std::string _finalRNGstate;
	std::vector<Real> var_;
	
	std::vector<VectorXr> scores_mat_;
	std::vector<VectorXr> loadings_mat_;
	std::vector<Real> lambda_PC_;
	std::vector<Real> variance_explained_;
	std::vector<Real> cumsum_percentage_;
	//std::vector<VectorXr> loadings_lambda_;
	//std::vector<VectorXr> scores_lambda_;
	MatrixXr datamatrixResiduals_ ;

	void computeBasisEvaluations();
	void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
	void computeDataMatrix(SpMat& DMat);
	void computeDataMatrixByIndices(SpMat& DMat);
	void computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput);
	void computeVarianceExplained();
	void computeCumulativePercentageExplained();
	void computeIterations(MatrixXr & datamatrixResiduals_,FPCAObject & FPCAinput, UInt lambda_index, UInt nnodes);
	void SetAndFixParameters();
	
	
	

public:
	//!A Constructor.
	MixedFEFPCABase(const MeshHandler<ORDER,mydim,ndim>& mesh, const FPCAData& fpcaData): mesh_(mesh),fpcaData_(fpcaData),isRcomputed_(false) {};

	
	virtual ~MixedFEFPCABase(){};	
	
	virtual  void apply()=0;
	
	//! A inline member that returns a VectorXr, returns the whole solution_.
	inline std::vector<VectorXr> const & getSolution() const{return solution_;};
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
	inline std::string const & getFinalRNGstate() const{return _finalRNGstate;}
	inline std::vector<Real> const & getVar() const{return var_;};
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCA : public MixedFEFPCABase<Integrator, ORDER, mydim, ndim>
{
public:
	MixedFEFPCA(const MeshHandler<ORDER, mydim, ndim>& mesh, const FPCAData& fpcaData):MixedFEFPCABase<Integrator, ORDER, mydim, ndim>(mesh, fpcaData){};
	
	virtual ~MixedFEFPCA(){};

	void apply();
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCAGCV: public MixedFEFPCABase<Integrator, ORDER, mydim, ndim>
{
protected:
	std::vector<VectorXr> loadings_lambda_;
	std::vector<VectorXr> scores_lambda_;
	
	std::vector<Real> dof_;
	std::vector<Real> GCV_;
	
	void computeDegreesOfFreedom(UInt output_index, Real lambda);
	void computeDegreesOfFreedomExact(UInt output_index, Real lambda);
	void computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);
	
	void computeGCV(FPCAObject& FPCAinput,UInt output_index);
	void computeDegreesOfFreedom(UInt output_index);
public:
	MixedFEFPCAGCV(const MeshHandler<ORDER, mydim, ndim>& mesh, const FPCAData& fpcaData):MixedFEFPCABase<Integrator, ORDER, mydim, ndim>(mesh, fpcaData){};
	
	virtual ~MixedFEFPCAGCV(){};

	void apply();
	
	inline std::vector<Real> const & getDOF() const{return dof_;};

};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCAKFold : public MixedFEFPCABase<Integrator, ORDER, mydim, ndim>
{
protected:
	std::vector<Real> KFold_;
	
	UInt nFolds;
	
	void computeKFolds(MatrixXr & datamatrixResiduals_, UInt lambda_index, UInt nnodes,UInt nFolds);
public:
	MixedFEFPCAKFold(const MeshHandler<ORDER, mydim, ndim>& mesh, const FPCAData& fpcaData):MixedFEFPCABase<Integrator, ORDER, mydim, ndim>(mesh, fpcaData){};
	
	virtual ~MixedFEFPCAKFold(){};

	void apply();
};


#include "mixedFEFPCA_imp.h"

#endif
