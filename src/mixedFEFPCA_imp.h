#ifndef __MIXEDFEFPCA_IMP_HPP__
#define __MIXEDFEFPCA_IMP_HPP__

#include <iostream>
#include<iterator>
#include <numeric>


//construct NW block of the system matrix when basis evaluation is necessary
//!! Depends on computeBasisEvaluations
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrix(SpMat& DMat)
{
		UInt nnodes = this->mesh_.num_nodes();
		DMat.resize(nnodes,nnodes);
		DMat = this->Psi_.transpose()*this->Psi_;
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = this->mesh_.num_nodes();
		UInt nlocations = this->inputData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		DMat.reserve(1);
		for (auto i = 0; i<nlocations; ++i)
		{
			auto index = this->inputData_.getObservationsIndices()[i];
			DMat.insert(index,index) = 1;
		}
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput)
{
	UInt nnodes = this->mesh_.num_nodes();
	UInt nlocations = this->inputData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

		if(this->inputData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = this->inputData_.getObservationsIndices()[i];
				rightHandData(index_i) = FPCAinput.getObservationData()[i];
			}
		}
		else
		{
			rightHandData=this->Psi_.transpose()*FPCAinput.getObservationData();
		}
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = this->mesh_.num_nodes();
	UInt nlocations = this->inputData_.getNumberofObservations();

	Eigen::SparseLU<SpMat> solver;
	solver.compute(this->coeffmatrix_);
	SpMat I(this->coeffmatrix_.rows(),this->coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = solver.solve(I);


	Real degrees=0;

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

	//std::cout<<"TRACE "<<degrees<<std::endl;

	this->dof_[output_index] = degrees;
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeVarianceExplained()
{	

	MatrixXr U_not_normalized(scores_mat_[0].size(),scores_mat_.size());
	for(UInt i=0;i<scores_mat_.size();i++)
		U_not_normalized.col(i)=scores_mat_[i];
	Eigen::HouseholderQR<MatrixXr> qr(U_not_normalized);
	MatrixXr R=qr.matrixQR().triangularView<Eigen::Upper>();
	variance_explained_.resize(this->inputData_.getNPC());
	for(UInt i=0;i<variance_explained_.size();i++)
	variance_explained_[i]=(R.diagonal()*R.diagonal().transpose()).diagonal()[i]/scores_mat_[0].size();
	}
	
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeCumulativePercentageExplained()
{	
	Eigen::JacobiSVD<MatrixXr> svd(this->inputData_.getDatamatrix(),Eigen::ComputeThinU|Eigen::ComputeThinV);
	MatrixXr U_ALL(this->inputData_.getDatamatrix().rows(),this->inputData_.getDatamatrix().rows());
	for(UInt i=0;i<svd.singularValues().rows();i++)
		U_ALL.col(i)=svd.matrixU().col(i)*svd.singularValues().diagonal()[i]*std::sqrt(svd.matrixV().col(i).transpose()*this->MMat_*svd.matrixV().col(i));
	Real TotVar=(U_ALL.transpose()*U_ALL).trace()/this->inputData_.getDatamatrix().rows();
	/*std::cout<<(U_ALL.transpose()*U_ALL)<<std::endl;
	std::cout<<(U_ALL.transpose()*U_ALL).trace()<<std::endl;
	std::cout<<this->inputData_.getDatamatrix().rows()<<std::endl;
	std::cout<<TotVar<<std::endl;
	*/
	
	cumsum_percentage_.resize(this->inputData_.getNPC());
	
	std::partial_sum(variance_explained_.begin(),variance_explained_.end(), cumsum_percentage_.begin());
	std::for_each(cumsum_percentage_.begin(), cumsum_percentage_.end(), [&TotVar](Real& i){i=i/TotVar;});
	}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper)
{
	UInt nnodes=this->mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	if(!this->inputData_.isLocationsByNodes())
	{
		this->computeBasisEvaluations();
		computeDataMatrix(this->DMat_);
	}else	
		computeDataMatrixByIndices(this->DMat_);

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, this->mesh_, fe, this->AMat_);
	Assembler::operKernel(mass, this->mesh_, fe, this->MMat_);
	
	///Fino a qui devo farlo una volta sola
	
	UInt niter=20;
	
	scores_mat_.resize(this->inputData_.getNPC());
	loadings_mat_.resize(this->inputData_.getNPC());
	lambda_PC_.resize(this->inputData_.getNPC());
	
	FPCAObject FPCAinput(this->inputData_.getDatamatrix());
	/*
	const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,Eigen::DontAlignCols,",","\n");
	
	std::string name1("MASSFPCA.csv");
	std::string name2("STIFFFPCA.csv");
	
	std::ofstream file1(name1.c_str());
	file1<<MatrixXr(this->MMat_).format(CSVFormat);
	
	std::ofstream file2(name2.c_str());
	file2<<MatrixXr(this->AMat_).format(CSVFormat);
	*/
	
	
for(UInt np=0;np<this->inputData_.getNPC();np++){
	//std::cout<<"Datamatrix"<<std::endl;
	
	//FPCAinput.printDatamatrix(std::cout);
	/*std::cout<<"Scores";
	FPCAinput.printScores(std::cout);
	std::cout<<"Loadings";
	FPCAinput.printLoadings(std::cout);
	std::cout<<"ObservationData";
	FPCAinput.printObservationData(std::cout);
	*/
	this->solution_.resize(this->inputData_.getLambda().size());
	this->dof_.resize(this->inputData_.getLambda().size());
	for(UInt i = 0; i<this->inputData_.getLambda().size(); ++i)
	{
		Real lambda = this->inputData_.getLambda()[i];
		//std::cout<<"Lambda: "<<lambda<<std::endl;
		SpMat AMat_lambda = (-lambda)*this->AMat_;
		SpMat MMat_lambda = (-lambda)*this->MMat_;
		this->buildCoeffMatrix(this->DMat_, AMat_lambda, MMat_lambda);
		
		
		/*
		std::string name("Coefmatrix.csv");
		std::ofstream file(name.c_str());
		file<<MatrixXr(this->coeffmatrix_).format(CSVFormat);
		*/
		for(UInt j=0;j<niter;j++)
		{
			FPCAinput.setObservationData();
			VectorXr rightHandData;
			computeRightHandData(rightHandData,FPCAinput);
			this->b_ = VectorXr::Zero(2*nnodes);
			this->b_.topRows(nnodes)=rightHandData;
			/*if(j==0){
			std::cout<<"Size: "<<std::endl;
			std::cout<<rightHandData.size()<<std::endl;
			std::cout<<"rightHandData: "<<std::endl;
			std::cout<<rightHandData<<std::endl;
			
			std::cout<<"B: "<<std::endl;
			std::cout<<this->b_<<std::endl;
			}*/
			/*if(j==0)
			{std::string name3("RIGHTHANDDATA.csv");
	
			std::ofstream file3(name3.c_str());
			file1<<this->b_.format(CSVFormat);
			}*/
			
			this-> template solve<SpLU>(i);
			if(this->inputData_.isLocationsByNodes())
				FPCAinput.setLoadings(nnodes, this->solution_[i]);
			else
				FPCAinput.setLoadingsPsi(nnodes, this->solution_[i],this->Psi_);
	
			FPCAinput.setScores();
			/*if(j==0){
			std::cout<<"Scores";
			FPCAinput.printScores(std::cout);
			std::cout<<"Loadings";
			FPCAinput.printLoadings(std::cout);
			std::cout<<"ObservationData";
			FPCAinput.printObservationData(std::cout);
			}*/
		}
		if(this->inputData_.computeDOF())
			computeDegreesOfFreedom(i);
		else
			this->dof_[i] = -1;
	}
	scores_mat_[np]=FPCAinput.getScores();
	loadings_mat_[np]=FPCAinput.getLoadings();
	lambda_PC_[np]=this->inputData_.getLambda()[0];
	
	//Devo settare la datamatrix togliendo i risultati ottenuti
	FPCAinput.newDatamatrix(scores_mat_[np],loadings_mat_[np]);
	
	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(loadings_mat_[np].transpose()*this->MMat_*loadings_mat_[np]);
	
	loadings_mat_[np]=loadings_mat_[np]/load_norm;
	
	scores_mat_[np]=scores_mat_[np]*load_norm;
	
	//std::cout<<"Scores*load transp"<<std::endl;
	//std::cout<<scores_mat_[np]*loadings_mat_[np].transpose()<<std::endl;
	//std::cout<<"NewDataMAT"<<std::endl;
	//std::cout<<this->inputData_.getDatamatrix()-scores_mat_[np]*loadings_mat_[np].transpose()<<std::endl;
	}
	computeVarianceExplained();
	computeCumulativePercentageExplained();
}

template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCA<FPCAData, Integrator, ORDER, mydim, ndim> : public MixedFEFPCABase<FPCAData, Integrator, ORDER, mydim, ndim>
{
public:
	MixedFEFPCA(const MeshHandler<ORDER, mydim, ndim>& mesh, const FPCAData& regressionData):MixedFEFPCABase<FPCAData, Integrator, ORDER, mydim, ndim>(mesh, regressionData){};

	void apply()
	{
		typedef EOExpr<Stiff> ETStiff; Stiff EStiff; ETStiff stiff(EStiff);
	    MixedFEFPCABase<FPCAData, Integrator, ORDER, mydim, ndim>::apply(stiff);
	}
};

#endif
