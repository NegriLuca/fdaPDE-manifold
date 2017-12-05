#ifndef __FPCAOBJECT_IMP_HPP__
#define __FPCAOBJECT_IMP_HPP__

#include<chrono>
#include<armadillo>


FPCAObject::FPCAObject(MatrixXr& datamatrix_)
{
	//Initialize loadings vector
	std::chrono::high_resolution_clock::time_point t13= std::chrono::high_resolution_clock::now();

	RedSVD::RedSVD<MatrixXr> svd0(datamatrix_,1); 
	std::chrono::high_resolution_clock::time_point t14= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration7 =t14-t13;
	
	std::cout<<"Time elapsed for computing svd Red: "<<duration7.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t15= std::chrono::high_resolution_clock::now();

	Eigen::JacobiSVD<MatrixXr> svd1(datamatrix_,Eigen::ComputeThinU|Eigen::ComputeThinV); 
	std::chrono::high_resolution_clock::time_point t16= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration8 =t16-t15;
	
	std::cout<<"Time elapsed for computing svd Eigen: "<<duration8.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t17= std::chrono::high_resolution_clock::now();

	Eigen::BDCSVD<MatrixXr> svd2(datamatrix_,Eigen::ComputeThinU|Eigen::ComputeThinV); 
	std::chrono::high_resolution_clock::time_point t18= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration9 =t18-t17;
	
	std::cout<<"Time elapsed for computing svd Eigen BDCSVD: "<<duration9.count()<<std::endl;
	
	
	arma::mat matrix=arma::mat(datamatrix_.data(),datamatrix_.rows(),datamatrix_.cols(),false,false);
	
	arma::mat U;
	arma::vec S;
	arma::mat V;
	
	std::chrono::high_resolution_clock::time_point t21= std::chrono::high_resolution_clock::now();
	
	arma::svd(U,S,V,matrix,"dc");

	std::chrono::high_resolution_clock::time_point t22= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration11 =t22-t21;
	
	std::cout<<"Time elapsed for computing svd Armadillo dc: "<<duration11.count()<<std::endl;
	
	
	arma::mat matrix1=arma::mat(datamatrix_.data(),datamatrix_.rows(),datamatrix_.cols(),false,false);
	
	arma::mat U1;
	arma::vec S1;
	arma::mat V1;
	
	std::chrono::high_resolution_clock::time_point t23= std::chrono::high_resolution_clock::now();
	
	arma::svd(U1,S1,V1,matrix1,"std");

	std::chrono::high_resolution_clock::time_point t24= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration12 =t24-t23;
	
	std::cout<<"Time elapsed for computing svd Armadillo std: "<<duration12.count()<<std::endl;
	
	arma::mat matrix2=arma::mat(datamatrix_.data(),datamatrix_.rows(),datamatrix_.cols(),false,false);
	
	arma::mat U2;
	arma::vec S2;
	arma::mat V2;
	
	std::chrono::high_resolution_clock::time_point t25= std::chrono::high_resolution_clock::now();
	
	arma::svd_econ(U2,S2,V2,matrix2,"both","std");

	std::chrono::high_resolution_clock::time_point t26= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration13 =t26-t25;
	
	std::cout<<"Time elapsed for computing svd econ Armadillo std: "<<duration13.count()<<std::endl;
	
	arma::mat matrix3=arma::mat(datamatrix_.data(),datamatrix_.rows(),datamatrix_.cols(),false,false);
	
	arma::mat U3;
	arma::vec S3;
	arma::mat V3;
	
	std::chrono::high_resolution_clock::time_point t27= std::chrono::high_resolution_clock::now();
	
	arma::svd_econ(U3,S3,V3,matrix3,"both","dc");

	std::chrono::high_resolution_clock::time_point t28= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration14 =t28-t27;
	
	std::cout<<"Time elapsed for computing svd econ Armadillo dc: "<<duration14.count()<<std::endl;
	
	arma::mat matrix4=arma::mat(datamatrix_.data(),datamatrix_.rows(),datamatrix_.cols(),false,false);
	
	arma::mat U4;
	arma::vec S4;
	arma::mat V4;
	
	std::chrono::high_resolution_clock::time_point t29= std::chrono::high_resolution_clock::now();
	
	svd_econ(U4,S4,V4,matrix4,1);

	std::chrono::high_resolution_clock::time_point t30= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration15 =t30-t29;
	
	std::cout<<"Time elapsed for computing svds: "<<duration15.count()<<std::endl;
		
	loadings_=svd1.matrixV().col(0);
	scores_=svd1.matrixU().col(0);
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

