#ifndef __MIXEDFEFPCA_IMP_HPP__
#define __MIXEDFEFPCA_IMP_HPP__

#include <iostream>
#include<iterator>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cmath>

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeBasisEvaluations(){

	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fpcaData_.getNumberofObservations();
	Real eps = 2.2204e-016,
		 tolerance = 100 * eps;

	Psi_.resize(nlocations, nnodes);
	//Psi_.reserve(Eigen::VectorXi::Constant(nlocations,ORDER*3));

	Triangle<ORDER*3, mydim, ndim> tri_activated;
	Eigen::Matrix<Real,ORDER * 3,1> coefficients;

	Real evaluator;
	for(UInt i=0; i<nlocations;i++)
	{
		tri_activated = mesh_.findLocationNaive(fpcaData_.getLocations()[i]);
		if(tri_activated.getId() == Identifier::NVAL)
		{
			#ifdef R_VERSION_
			Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			#else
			std::cout << "WARNING: Observation " << i+1 <<" is not in the domain\n";
			#endif
		}else
		{
			for(UInt node = 0; node < ORDER*3 ; ++node)
			{
				coefficients = Eigen::Matrix<Real,ORDER * 3,1>::Zero();
				coefficients(node) = 1; //Activates only current base
				evaluator = evaluate_point<ORDER, mydim, ndim>(tri_activated, fpcaData_.getLocations()[i], coefficients);
				Psi_.insert(i, tri_activated[node].getId()) = evaluator;
			}
		}
	}

	Psi_.prune(tolerance);
	Psi_.makeCompressed();
}


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
{
	UInt nnodes = mesh_.num_nodes();

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(DMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
	  }
	for (int k=0; k<MMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(MMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	coeffmatrix_.setZero();
	coeffmatrix_.resize(2*nnodes,2*nnodes);
	coeffmatrix_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	coeffmatrix_.makeCompressed();
}

//construct NW block of the system matrix when basis evaluation is necessary
//!! Depends on computeBasisEvaluations
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeDataMatrix(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		DMat.resize(nnodes,nnodes);
		DMat = Psi_.transpose()*Psi_;
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		UInt nlocations = fpcaData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		DMat.reserve(1);
		for (auto i = 0; i<nlocations; ++i)
		{
			auto index = fpcaData_.getObservationsIndices()[i];
			DMat.insert(index,index) = 1;
		}
}


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fpcaData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

		if(fpcaData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = fpcaData_.getObservationsIndices()[i];
				rightHandData(index_i) = FPCAinput.getObservationData()[i];
			}
		}
		else
		{
			rightHandData=Psi_.transpose()*FPCAinput.getObservationData();
		}
}





template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeVarianceExplained()
{	

	MatrixXr U_not_normalized(scores_mat_[0].size(),scores_mat_.size());
	for(UInt i=0;i<scores_mat_.size();i++)
		U_not_normalized.col(i)=scores_mat_[i];
	Eigen::HouseholderQR<MatrixXr> qr(U_not_normalized);
	MatrixXr R=qr.matrixQR().triangularView<Eigen::Upper>();
	variance_explained_.resize(fpcaData_.getNPC());
	for(UInt i=0;i<variance_explained_.size();i++)
	variance_explained_[i]=(R.diagonal()*R.diagonal().transpose()).diagonal()[i]/scores_mat_[0].size();
	}
	
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeCumulativePercentageExplained()
{	
	Eigen::BDCSVD<MatrixXr> svd(fpcaData_.getDatamatrix(),Eigen::ComputeThinU|Eigen::ComputeThinV);
	MatrixXr U_ALL(fpcaData_.getDatamatrix().rows(),fpcaData_.getDatamatrix().rows());
	for(UInt i=0;i<svd.singularValues().rows();i++)
		U_ALL.col(i)=svd.matrixU().col(i)*svd.singularValues().diagonal()[i]*std::sqrt(svd.matrixV().col(i).transpose()*MMat_*svd.matrixV().col(i));
	Real TotVar=(U_ALL.transpose()*U_ALL).trace()/fpcaData_.getDatamatrix().rows();
	
	cumsum_percentage_.resize(fpcaData_.getNPC());
	
	std::partial_sum(variance_explained_.begin(),variance_explained_.end(), cumsum_percentage_.begin());
	std::for_each(cumsum_percentage_.begin(), cumsum_percentage_.end(), [&TotVar](Real& i){i=i/TotVar;});
	}
	


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeIterations(MatrixXr & datamatrixResiduals_,FPCAObject & FPCAinput, UInt lambda_index, UInt nnodes)
{
		Real lambda = fpcaData_.getLambda()[lambda_index];
		SpMat AMat_lambda = (-lambda)*AMat_;
		SpMat MMat_lambda = (-lambda)*MMat_;
		buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);

		sparseSolver_.analyzePattern(coeffmatrix_);
		sparseSolver_.factorize(coeffmatrix_);
		solution_[lambda_index].resize(coeffmatrix_.rows());
		
		UInt niter=20;

		for(auto j=0;j<niter;j++)
		{	
			FPCAinput.setObservationData(datamatrixResiduals_);

			VectorXr rightHandData;
			computeRightHandData(rightHandData,FPCAinput);
			b_ = VectorXr::Zero(2*nnodes);
			b_.topRows(nnodes)=rightHandData;
			
			
			solution_[lambda_index]=sparseSolver_.solve(b_);
			if(sparseSolver_.info()!=Eigen::Success)
			{
			//std::cerr<<"solving failed!"<<std::endl;
			}
			
			if(fpcaData_.isLocationsByNodes())
				FPCAinput.setLoadings(nnodes, solution_[lambda_index]);
			else
				FPCAinput.setLoadingsPsi(nnodes, solution_[lambda_index],Psi_);
	
			
			FPCAinput.setScores(datamatrixResiduals_);

		}

}



template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::SetAndFixParameters()
{	
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	if(!fpcaData_.isLocationsByNodes())
	{
		computeBasisEvaluations();
		computeDataMatrix(DMat_);
	}else	
		computeDataMatrixByIndices(DMat_);

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	Assembler::operKernel(stiff, mesh_, fe, AMat_);
	Assembler::operKernel(mass, mesh_, fe, MMat_);
	
	scores_mat_.resize(fpcaData_.getNPC());
	loadings_mat_.resize(fpcaData_.getNPC());
	lambda_PC_.resize(fpcaData_.getNPC());
	
	
	datamatrixResiduals_ = fpcaData_.getDatamatrix();
	solution_.resize(fpcaData_.getLambda().size());
}	


///CLASS MIXEDFEFPCA
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCA<Integrator,ORDER, mydim, ndim>::apply()
{  

   MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::SetAndFixParameters();

   for(auto np=0;np<this->fpcaData_.getNPC();np++){
   
		UInt i=0;
		FPCAObject FPCAinput(this->datamatrixResiduals_);

		MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeIterations(this->datamatrixResiduals_,FPCAinput,i,this->mesh_.num_nodes());

	this->scores_mat_[np]=FPCAinput.getScores();
	this->loadings_mat_[np]=FPCAinput.getLoadings();
	this->lambda_PC_[np]=this->fpcaData_.getLambda()[i];
	
	//Devo settare la datamatrix togliendo i risultati ottenuti
	this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();
	
	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);
	
	this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;
	
	this->scores_mat_[np]=this->scores_mat_[np]*load_norm;	
	
	}
	
	MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeVarianceExplained();
	MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeCumulativePercentageExplained();
}

///CLASS MIXEDFEFPCAGCV
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<Integrator,ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = this->mesh_.num_nodes();
	UInt nlocations = this->fpcaData_.getNumberofObservations();

	SpMat I(this->coeffmatrix_.rows(),this->coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = this->sparseSolver_.solve(I);


	Real degrees=0;

		if(this->fpcaData_.isLocationsByNodes())
		{
			VectorXr d = coeff_inv.diagonal();

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = this->fpcaData_.getObservationsIndices()[i];
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

	//std::cout<<"TRACE "<<degrees<<std::endl;

	dof_[output_index] = degrees;
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<Integrator,ORDER, mydim, ndim>::computeGCV(FPCAObject& FPCAinput,UInt output_index)
{
	UInt s= this->fpcaData_.getNumberofObservations();
	VectorXr zhat=FPCAinput.getObservationData();
	Real norm_squared=(zhat-FPCAinput.getLoadings()).transpose()*(zhat-FPCAinput.getLoadings());
	if(s-dof_[output_index]<0){
		#ifdef R_VERSION_
			Rprintf("WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is: %d \n", this->fpcaData_.getLambda()[output_index]);
			#else
			std::cout << "WARNING: Some values of the trace of the matrix S('lambda') are inconstistent. This might be due to ill-conditioning of the linear system. Try increasing value of 'lambda'.Value of 'lambda' that produces an error is:" << this->fpcaData_.getLambda()[output_index] <<"\n";
			#endif
			}
	Real stderror=norm_squared/(s-dof_[output_index]);
	GCV_[output_index]=(s/(s-dof_[output_index]))*stderror;
	//GCV_[output_index]=norm_squared/(s-(1-1/s*dof_[output_index])*(1-1/s*dof_[output_index]));
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAGCV<Integrator,ORDER, mydim, ndim>::apply()
{
   MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::SetAndFixParameters();
   dof_.resize(this->fpcaData_.getLambda().size());
   GCV_.resize(this->fpcaData_.getLambda().size());
   loadings_lambda_.resize(this->fpcaData_.getLambda().size());
   scores_lambda_.resize(this->fpcaData_.getLambda().size());
   for(auto np=0;np<this->fpcaData_.getNPC();np++){

	for(auto i = 0; i<this->fpcaData_.getLambda().size(); ++i)
	{	
		FPCAObject FPCAinput(this->datamatrixResiduals_);
		MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeIterations(this->datamatrixResiduals_,FPCAinput,i,this->mesh_.num_nodes());
		loadings_lambda_[i]=FPCAinput.getLoadings();
		scores_lambda_[i]=FPCAinput.getScores();
		
			if(np==0) computeDegreesOfFreedom(i);
			computeGCV(FPCAinput,i);	
	}
	
	UInt index_best_GCV=std::distance(GCV_.begin(),std::min_element(GCV_.begin(),GCV_.end()));
	
	this->scores_mat_[np]=this->scores_lambda_[index_best_GCV];
	this->loadings_mat_[np]=this->loadings_lambda_[index_best_GCV];
	this->lambda_PC_[np]=this->fpcaData_.getLambda()[index_best_GCV];
	
	//Devo settare la datamatrix togliendo i risultati ottenuti
	this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();
	
	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);
	
	this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;
	
	this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
	
	}
	MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeVarianceExplained();
	MixedFEFPCABase<Integrator,ORDER, mydim,ndim>::computeCumulativePercentageExplained();
}


///CLASS MIXEDFEFPCAKFOLD
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAKFold<Integrator,ORDER, mydim, ndim>::computeKFolds(MatrixXr & datamatrixResiduals_, UInt lambda_index, UInt nnodes, UInt nFolds)
{
		Real lambda = this->fpcaData_.getLambda()[lambda_index];
		SpMat AMat_lambda = (-lambda)*this->AMat_;
		SpMat MMat_lambda = (-lambda)*this->MMat_;
		this->buildCoeffMatrix(this->DMat_, AMat_lambda, MMat_lambda);

		this->sparseSolver_.analyzePattern(this->coeffmatrix_);
		this->sparseSolver_.factorize(this->coeffmatrix_);
		this->solution_[lambda_index].resize(this->coeffmatrix_.rows());
		
		UInt niter=20;
		std::vector<UInt> indices_valid;
		for(auto k=0;k<nFolds;k++)
		{
			UInt length_chunk=floor(datamatrixResiduals_.rows()/nFolds);
			indices_valid.resize(datamatrixResiduals_.rows());
			std::chrono::high_resolution_clock::time_point t15= std::chrono::high_resolution_clock::now();
			std::iota(indices_valid.begin(),indices_valid.begin()+length_chunk,k*length_chunk);
			if(k==0) std::iota(indices_valid.begin()+length_chunk,indices_valid.end(),(k+1)*length_chunk);
			else {
				std::iota(indices_valid.begin()+length_chunk,indices_valid.begin()+(k+1)*length_chunk,0);
				std::iota(indices_valid.begin()+(k+1)*length_chunk,indices_valid.end(),(k+1)*length_chunk);
			}
			
			std::chrono::high_resolution_clock::time_point t16= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration8 =t16-t15;
	
	std::cout<<"Time elapsed for computing iota : "<<duration8.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t17= std::chrono::high_resolution_clock::now();
		
			VectorXi indices_v=Eigen::Map<VectorXi,Eigen::Unaligned> (indices_valid.data(),indices_valid.size());
			Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(indices_v);
		
			MatrixXr X_train=(perm*datamatrixResiduals_).bottomRows(datamatrixResiduals_.rows()-length_chunk);
			MatrixXr X_clean_train=(perm*datamatrixResiduals_).bottomRows(datamatrixResiduals_.rows()-length_chunk);
			MatrixXr X_valid=(perm*datamatrixResiduals_).topRows(length_chunk);
		
			std::chrono::high_resolution_clock::time_point t18= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration9 =t18-t17;
	
	std::cout<<"Time elapsed for computing permutation : "<<duration9.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t21= std::chrono::high_resolution_clock::now();
		
			VectorXr X_train_mean=X_train.colwise().mean();
			VectorXr X_clean_train_mean=X_clean_train.colwise().mean();
			VectorXr X_valid_mean=X_valid.colwise().mean();
		
			VectorXr ones=VectorXr::Constant(X_train.rows(),1,1);
			X_train=X_train-ones*X_train_mean.transpose();
			X_clean_train=X_clean_train-ones*X_clean_train_mean.transpose();
			X_valid=X_valid-ones*X_valid_mean.transpose();
			
			std::chrono::high_resolution_clock::time_point t22= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration11 =t22-t21;
	
	std::cout<<"Time elapsed for computing permutation : "<<duration11.count()<<std::endl;
		
		
			FPCAObject FPCAinputKF(X_clean_train);
			for(auto j=0;j<niter;j++)
			{	
				FPCAinputKF.setObservationData(X_train);

				VectorXr rightHandData;
				this->computeRightHandData(rightHandData,FPCAinputKF);
				this->b_ = VectorXr::Zero(2*nnodes);
				this->b_.topRows(nnodes)=rightHandData;
			
			
				this->solution_[lambda_index]=this->sparseSolver_.solve(this->b_);
				if(this->sparseSolver_.info()!=Eigen::Success)
				{
				//std::cerr<<"solving failed!"<<std::endl;
				}
			
				if(this->fpcaData_.isLocationsByNodes())
					FPCAinputKF.setLoadings(nnodes, this->solution_[lambda_index]);
				else
					FPCAinputKF.setLoadingsPsi(nnodes, this->solution_[lambda_index],this->Psi_);
	
			
				FPCAinputKF.setScores(X_clean_train);

			}
			
			std::chrono::high_resolution_clock::time_point t19= std::chrono::high_resolution_clock::now();
			Real U_hat_const=FPCAinputKF.getLoadings().squaredNorm() + lambda* (this->solution_[lambda_index].bottomRows(nnodes)).transpose()*this->MMat_*this->solution_[lambda_index].bottomRows(nnodes);
			VectorXr U_hat_valid=(X_valid*FPCAinputKF.getLoadings())/U_hat_const;
		Real diffCV=(X_valid-U_hat_valid*FPCAinputKF.getLoadings().transpose()).squaredNorm()/(X_valid.rows()*X_valid.cols());
		Real sumCV=(X_valid+U_hat_valid*FPCAinputKF.getLoadings().transpose()).squaredNorm()/(X_valid.rows()*X_valid.cols());
		KFold_[lambda_index]+=std::min(diffCV,sumCV);
		
		std::chrono::high_resolution_clock::time_point t20= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration10 =t20-t19;
	
	std::cout<<"Time elapsed for computing CV for fold K : "<<duration10.count()<<std::endl;
		}

}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCAKFold<Integrator,ORDER, mydim, ndim>::apply()
{
   MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::SetAndFixParameters();
   nFolds=this->fpcaData_.getNFolds();
   KFold_.resize(this->fpcaData_.getLambda().size());
   for(auto np=0;np<this->fpcaData_.getNPC();np++){
	std::fill(KFold_.begin(),KFold_.end(),0);
	for(auto i = 0; i<this->fpcaData_.getLambda().size(); ++i)
	{	
		FPCAObject FPCAinput(this->datamatrixResiduals_);
		computeKFolds(this->datamatrixResiduals_, i, this->mesh_.num_nodes(), nFolds);	
	}
	UInt index_best_KF=std::distance(KFold_.begin(),std::min_element(KFold_.begin(),KFold_.end()));
	FPCAObject FPCAinput(this->datamatrixResiduals_);
	MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeIterations(this->datamatrixResiduals_,FPCAinput,index_best_KF,this->mesh_.num_nodes());
	
	this->scores_mat_[np]=FPCAinput.getScores();
	this->loadings_mat_[np]=FPCAinput.getLoadings();
	this->lambda_PC_[np]=this->fpcaData_.getLambda()[index_best_KF];
	
	//Devo settare la datamatrix togliendo i risultati ottenuti
	this->datamatrixResiduals_=this->datamatrixResiduals_-this->scores_mat_[np]*this->loadings_mat_[np].transpose();
	
	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(this->loadings_mat_[np].transpose()*this->MMat_*this->loadings_mat_[np]);
	
	this->loadings_mat_[np]=this->loadings_mat_[np]/load_norm;
	
	this->scores_mat_[np]=this->scores_mat_[np]*load_norm;
	}
	
	MixedFEFPCABase<Integrator,ORDER, mydim, ndim>::computeVarianceExplained();
	MixedFEFPCABase<Integrator,ORDER, mydim,ndim>::computeCumulativePercentageExplained();
}


#endif
