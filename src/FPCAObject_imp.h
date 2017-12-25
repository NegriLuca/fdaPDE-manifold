#ifndef __FPCAOBJECT_IMP_HPP__
#define __FPCAOBJECT_IMP_HPP__

#include<chrono>

FPCAObject::FPCAObject(const MatrixXr& datamatrix_)
{
	//Initialize loadings vector
	RedSVD::RedSVD<MatrixXr> svd(datamatrix_,1); 
		
	loadings_=svd.matrixV().col(0);
	scores_=svd.matrixU().col(0);
}


void FPCAObject::printScores(std::ostream & out) const
{

	for(auto i=0;i<scores_.size(); i++)
	{
		out<<scores_(i)<<"\t";
	}
	out<<std::endl;
}

void FPCAObject::printLoadings(std::ostream & out) const
{

	for(auto i=0;i<loadings_.size(); i++)
	{
		out<<loadings_(i)<<"\t";
	}
	out<<std::endl;
}

void FPCAObject::printObservationData(std::ostream & out) const
{

	for(auto i=0;i<ObservationData_.size(); i++)
	{
		out<<ObservationData_(i)<<"\t";
	}
	out<<std::endl;
}

void FPCAObject::setScores(const MatrixXr& datamatrix_)
{
	scores_=datamatrix_*loadings_;
	scores_=scores_/scores_.norm();
}

/*void FPCAObject::setObservationData(const MatrixXr& datamatrix_, const SpMat& psi_)
{
	ObservationData_=psi_.transpose()*datamatrix_.transpose()*scores_;
}*/

void FPCAObject::setObservationData(const MatrixXr& datamatrix_)
{
	ObservationData_=datamatrix_.transpose()*scores_;
}

void FPCAObject::setLoadingsPsi(UInt nnodes, const VectorXr& f_sol,const SpMat& psi_)
{	
	//VectorXr load_=psi_*f_sol.topRows(nnodes);
	VectorXr load_=psi_.transpose()*f_sol.topRows(nnodes);
	loadings_=load_;
	//std::cout<<"Load dim:"<<loadings_.size()<<std::endl;
}

void FPCAObject::setLoadings(UInt nnodes, const VectorXr& f_sol,const std::vector<UInt>& obs_indices)
{
	VectorXr loadings_full_=f_sol.topRows(nnodes);
	for(auto i=0;i<obs_indices.size();i++)  loadings_(i)=loadings_full_(obs_indices[i]);
}

void FPCAObject::finalizeLoadings(const std::vector<UInt>& obs_indices,UInt nlocations)
{
	VectorXr finalize_=VectorXr::Zero(nlocations);
	for(auto i=0;i<obs_indices.size();i++)
		finalize_(obs_indices[i])=loadings_(i);
	loadings_=finalize_;
}

#endif

