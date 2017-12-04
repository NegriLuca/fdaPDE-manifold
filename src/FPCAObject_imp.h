#ifndef __FPCAOBJECT_IMP_HPP__
#define __FPCAOBJECT_IMP_HPP__

#include<chrono>

FPCAObject::FPCAObject(const MatrixXr& datamatrix_)
{
	//Initialize loadings vector
	std::chrono::high_resolution_clock::time_point t13= std::chrono::high_resolution_clock::now();

	RedSVD::RedSVD<MatrixXr> svd(datamatrix_,1); 
	std::chrono::high_resolution_clock::time_point t14= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration7 =t14-t13;
	
	std::cout<<"Time elapsed for computing svd : "<<duration7.count()<<std::endl;	
	loadings_=svd.matrixV().col(0);
	scores_=svd.matrixU().col(0);
}

/*
void FPCAObject::newDatamatrix(const VectorXr& scores,const VectorXr& loadings)
{    	
	//std::cout<<"RIGHE DATA: "<<datamatrix_.rows()<<std::endl;
	//std::cout<<"COLONNE DATA: "<<datamatrix_.cols()<<std::endl;
	
	//MatrixXr aa=loadings*scores.transpose();
	
	//std::cout<<"RIGHE: "<<aa.rows()<<std::endl;
	//std::cout<<"COLONNE: "<<aa.cols()<<std::endl;
	
	
	datamatrix_=datamatrix_- scores*loadings.transpose();
	
	//Initialize loadings vector
	Eigen::JacobiSVD<MatrixXr> svd(datamatrix_,Eigen::ComputeThinU|Eigen::ComputeThinV);
	loadings_=svd.matrixV().col(0);
	scores_=svd.matrixU().col(0);
}
*/


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

/*
void FPCAObject::printDatamatrix(std::ostream & out) const
{

	for(auto i=0;i<datamatrix_.rows(); i++)
	{
		for(auto j=0;j<datamatrix_.cols(); i++)
		{
			out<<datamatrix_(i,j)<<"\t";
		}
	}
	out<<std::endl;
}
*/


void FPCAObject::setScores(const MatrixXr& datamatrix_)
{
	scores_=datamatrix_*loadings_;
	scores_=scores_/scores_.norm();
}

void FPCAObject::setObservationData(const MatrixXr& datamatrix_)
{
	ObservationData_=datamatrix_.transpose()*scores_;
}

void FPCAObject::setLoadingsPsi(UInt nnodes, const VectorXr& f_sol,const SpMat& psi_)
{
	loadings_=psi_.transpose()*f_sol.topRows(nnodes);
}

void FPCAObject::setLoadings(UInt nnodes, const VectorXr& f_sol)
{
	loadings_=f_sol.topRows(nnodes);
}



#endif

