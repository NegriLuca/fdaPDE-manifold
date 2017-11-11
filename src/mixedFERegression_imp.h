#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeProjOnCovMatrix()
{
	//std::cout<<"Computing Projection Matrix"<<std::endl;
	UInt nlocations = this->inputData_.getNumberofObservations();

	//Set covariate matrix
	MatrixXr W(this->inputData_.getCovariates());
	if(this->inputData_.isLocationsByNodes())
	{
		MatrixXr W_reduced(this->inputData_.getNumberofObservations(), W.cols());
		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = this->inputData_.getObservationsIndices()[i];
			for (auto j=0; j<W.cols();++j)
			{
				W_reduced(i,j) = W(index_i,j);
			}
		}
		W = W_reduced;
	}

	//std::cout << W << std::endl;
	MatrixXr WTW(W.transpose()*W);
	P_.resize(W.rows(),W.rows());
	P_ = -W*WTW.ldlt().solve(W.transpose());
	for (int i=0; i<P_.rows();++i)
	{
		P_(i,i) += 1;
	}
	//std::cout << P_ << std::endl;
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues)
{
	UInt nnodes = this->mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = this->inputData_.getDirichletIndices();
	const std::vector<Real>& bc_values = this->inputData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	UInt id1,id3;
	for( auto i=0; i<nbc_indices; i++)
	 {
		id1=bcindex[i];
		id3=id1+nnodes;

		this->coeffmatrix_.coeffRef(id1,id1) = dirichlet_penalization;
		this->coeffmatrix_.coeffRef(id3,id3) = dirichlet_penalization;


		this->b_(id1) += bc_values[i]*dirichlet_penalization;
		this->b_(id3) = 0;
	 }
	this->coeffmatrix_.makeCompressed();
}

//construct NW block of the system matrix when basis evaluation is necessary
//!! Depends on computeBasisEvaluations and computeProjOnCovMatrix
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrix(SpMat& DMat)
{
		UInt nnodes = this->mesh_.num_nodes();
		DMat.resize(nnodes,nnodes);

		if (this->inputData_.getCovariates().rows() == 0)
			DMat = this->Psi_.transpose()*this->Psi_;
		else
			DMat = (SpMat(this->Psi_.transpose())*P_*this->Psi_).sparseView();
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations and computeProjOnCovMatrix
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = this->mesh_.num_nodes();
		UInt nlocations = this->inputData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		if (this->inputData_.getCovariates().rows() == 0)
		{
			DMat.reserve(1);
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index = this->inputData_.getObservationsIndices()[i];
				DMat.insert(index,index) = 1;
			}
		}
		else
		{
			//May be inefficient
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index_i = this->inputData_.getObservationsIndices()[i];
				for (auto j = 0; j<nlocations; ++j)
				{
					auto index_j = this->inputData_.getObservationsIndices()[j];
					DMat.insert(index_i,index_j) = P_(i,j);
				}
			}
		}
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = this->mesh_.num_nodes();
	UInt nlocations = this->inputData_.getNumberofObservations();
	//rightHandData.resize(nnodes);
	rightHandData = VectorXr::Zero(nnodes);

	if(this->inputData_.getCovariates().rows() == 0)
	{
		if(this->inputData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = this->inputData_.getObservationsIndices()[i];
				rightHandData(index_i) = this->inputData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData=this->Psi_.transpose()*this->inputData_.getObservations();
		}
	}
	else
	{
		if(this->inputData_.isLocationsByNodes())
		{
			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = this->inputData_.getObservationsIndices()[i];
				rightHandData(index_i) = (P_).row(i) * this->inputData_.getObservations();
			}
		}
		else
		{
			rightHandData=this->Psi_.transpose()*P_*this->inputData_.getObservations();
		}
	}
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = this->mesh_.num_nodes();
	UInt nlocations = this->inputData_.getNumberofObservations();

	Eigen::SparseLU<SpMat> solver;
	solver.compute(this->coeffmatrix_);
	SpMat I(this->coeffmatrix_.rows(),this->coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = solver.solve(I);


	Real degrees=0;

	if(this->inputData_.getCovariates().rows() == 0)
	{
		if(this->inputData_.isLocationsByNodes())
		{
			VectorXr d = coeff_inv.diagonal();

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = this->inputData_.getObservationsIndices()[i];
				degrees+=d(index_i);
			}
		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = this->Psi_*An*this->Psi_.transpose();
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}
	}
	else
	{	//Important
		degrees = this->inputData_.getCovariates().cols();

		if(this->inputData_.isLocationsByNodes())
		{
			MatrixXr S(coeff_inv);

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = this->inputData_.getObservationsIndices()[i];
				for (auto j=0; j<nlocations;++j)
				{
					auto index_j = this->inputData_.getObservationsIndices()[j];

					degrees+=S(index_i,index_j)*P_(j,i);
				}
			}

		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = this->Psi_*An*this->Psi_.transpose()*P_;
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}
	}


	//std::cout<<"TRACE "<<degrees<<std::endl;

	this->dof_[output_index] = degrees;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper)
{
	UInt nnodes=this->mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	if(!this->inputData_.isLocationsByNodes())
		this->computeBasisEvaluations();
	if(!(this->inputData_.getCovariates().rows() == 0))
		computeProjOnCovMatrix();
	if(!this->inputData_.isLocationsByNodes())
		computeDataMatrix(this->DMat_);
	else	
		computeDataMatrixByIndices(this->DMat_);

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, this->mesh_, fe, this->AMat_);
	Assembler::operKernel(mass, this->mesh_, fe, this->MMat_);
	VectorXr rightHandData;
	computeRightHandData(rightHandData);
	this->b_ = VectorXr::Zero(2*nnodes);
	this->b_.topRows(nnodes)=rightHandData;
	this->solution_.resize(this->inputData_.getLambda().size());
	this->dof_.resize(this->inputData_.getLambda().size());

	for(UInt i = 0; i<this->inputData_.getLambda().size(); ++i)
	{
		Real lambda = this->inputData_.getLambda()[i];
		SpMat AMat_lambda = (-lambda)*this->AMat_;
		SpMat MMat_lambda = (-lambda)*this->MMat_;
		this->buildCoeffMatrix(this->DMat_, AMat_lambda, MMat_lambda);
		//std::cout << coeffmatrix_ << std::endl;
		//Applying boundary conditions if necessary
		if(this->inputData_.getDirichletIndices().size() != 0)
			applyDirichletBC(this->inputData_.getDirichletIndices(), this->inputData_.getDirichletValues());

		this-> template solve<SpLU>(i);
		if(this->inputData_.computeDOF())
			computeDegreesOfFreedom(i);
		else
			this->dof_[i] = -1;
	}
}


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionData, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionData& regressionData):MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::apply(stiff);
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataElliptic, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataElliptic& regressionData):MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

	    const Real& c = this->inputData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->inputData_.getK();
	    const Eigen::Matrix<Real,2,1>& beta = this->inputData_.getBeta();

	    MixedFERegressionBase<RegressionDataElliptic, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
	}
};

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const RegressionDataEllipticSpaceVarying& regressionData):MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
	if(mydim!=2 || ndim !=2){

	#ifdef R_VERSION_
		Rprintf("ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2");
	#else
		std::cout << "ERROR: these dimensions are not yet implemented, for the moment smoothEllipticPDE is available for mydim=ndim=2\n";
	#endif

	}else{
		typedef EOExpr<Mass> ETMass;   Mass EMass;   ETMass mass(EMass);
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
		typedef EOExpr<Grad> ETGrad;   Grad EGrad;   ETGrad grad(EGrad);

		const Reaction& c = this->inputData_.getC();
		const Diffusivity& K = this->inputData_.getK();
		const Advection& beta = this->inputData_.getBeta();

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
	}
};

/////////////////////////////////////////////
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFERegression<FPCAData, Integrator, ORDER, mydim, ndim> : public MixedFERegressionBase<FPCAData, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFERegression(const MeshHandler<ORDER, mydim, ndim>& mesh, const FPCAData& regressionData):MixedFERegressionBase<FPCAData, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFERegressionBase<FPCAData, Integrator, ORDER, mydim, ndim>::apply(stiff);
	}
};


#endif
