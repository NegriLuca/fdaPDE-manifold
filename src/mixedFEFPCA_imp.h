#ifndef __MIXEDFEFPCA_IMP_HPP__
#define __MIXEDFEFPCA_IMP_HPP__

#include <iostream>
#include<iterator>
#include <numeric>

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeBasisEvaluations(){

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


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
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
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrix(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		DMat.resize(nnodes,nnodes);
		DMat = Psi_.transpose()*Psi_;
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrixByIndices(SpMat& DMat)
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


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeRightHandData(VectorXr& rightHandData,FPCAObject& FPCAinput)
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


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = fpcaData_.getNumberofObservations();

	Eigen::SparseLU<SpMat> solver;
	solver.compute(coeffmatrix_);
	SpMat I(coeffmatrix_.rows(),coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = solver.solve(I);


	Real degrees=0;

		if(fpcaData_.isLocationsByNodes())
		{
			VectorXr d = coeff_inv.diagonal();

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = fpcaData_.getObservationsIndices()[i];
				degrees+=d(index_i);
			}
		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = Psi_*An*Psi_.transpose();
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}

	//std::cout<<"TRACE "<<degrees<<std::endl;

	dof_[output_index] = degrees;
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeVarianceExplained()
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
	
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::computeCumulativePercentageExplained()
{	
	Eigen::JacobiSVD<MatrixXr> svd(fpcaData_.getDatamatrix(),Eigen::ComputeThinU|Eigen::ComputeThinV);
	MatrixXr U_ALL(fpcaData_.getDatamatrix().rows(),fpcaData_.getDatamatrix().rows());
	for(UInt i=0;i<svd.singularValues().rows();i++)
		U_ALL.col(i)=svd.matrixU().col(i)*svd.singularValues().diagonal()[i]*std::sqrt(svd.matrixV().col(i).transpose()*MMat_*svd.matrixV().col(i));
	Real TotVar=(U_ALL.transpose()*U_ALL).trace()/fpcaData_.getDatamatrix().rows();
	/*std::cout<<(U_ALL.transpose()*U_ALL)<<std::endl;
	std::cout<<(U_ALL.transpose()*U_ALL).trace()<<std::endl;
	std::cout<<this->inputData_.getDatamatrix().rows()<<std::endl;
	std::cout<<TotVar<<std::endl;
	*/
	
	cumsum_percentage_.resize(fpcaData_.getNPC());
	
	std::partial_sum(variance_explained_.begin(),variance_explained_.end(), cumsum_percentage_.begin());
	std::for_each(cumsum_percentage_.begin(), cumsum_percentage_.end(), [&TotVar](Real& i){i=i/TotVar;});
	}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFEFPCABase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper)
{	
	UInt nnodes=mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;

	if(!fpcaData_.isLocationsByNodes())
	{
		computeBasisEvaluations();
		computeDataMatrix(DMat_);
	}else	
		computeDataMatrixByIndices(DMat_);

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, AMat_);
	Assembler::operKernel(mass, mesh_, fe, MMat_);
	
	///Fino a qui devo farlo una volta sola
	UInt niter=20;
	
	scores_mat_.resize(fpcaData_.getNPC());
	loadings_mat_.resize(fpcaData_.getNPC());
	lambda_PC_.resize(fpcaData_.getNPC());
	
	FPCAObject FPCAinput(fpcaData_.getDatamatrix());
	/*
	const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,Eigen::DontAlignCols,",","\n");
	
	std::string name1("MASSFPCA.csv");
	std::string name2("STIFFFPCA.csv");
	
	std::ofstream file1(name1.c_str());
	file1<<MatrixXr(this->MMat_).format(CSVFormat);
	
	std::ofstream file2(name2.c_str());
	file2<<MatrixXr(this->AMat_).format(CSVFormat);
	*/
	
	//std::chrono::high_resolution_clock::time_point t11= std::chrono::high_resolution_clock::now();
	
for(auto np=0;np<fpcaData_.getNPC();np++){
	//std::cout<<"Datamatrix"<<std::endl;
	
	//FPCAinput.printDatamatrix(std::cout);
	/*std::cout<<"Scores";
	FPCAinput.printScores(std::cout);
	std::cout<<"Loadings";
	FPCAinput.printLoadings(std::cout);
	std::cout<<"ObservationData";
	FPCAinput.printObservationData(std::cout);
	*/
	
	//std::chrono::high_resolution_clock::time_point t9= std::chrono::high_resolution_clock::now();
	
	solution_.resize(fpcaData_.getLambda().size());
	dof_.resize(fpcaData_.getLambda().size());
	for(auto i = 0; i<fpcaData_.getLambda().size(); ++i)
	{
		Real lambda = fpcaData_.getLambda()[i];
		//std::cout<<"Lambda: "<<lambda<<std::endl;
		SpMat AMat_lambda = (-lambda)*AMat_;
		SpMat MMat_lambda = (-lambda)*MMat_;
		buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);
		
		/*Eigen::BiCGSTAB<SpMat,Eigen::IncompleteLUT<Real>> solver;
		Eigen::BiCGSTAB<SpMat> solver2;
		solver.analyzePattern(this->coeffmatrix_);
		solver.factorize(this->coeffmatrix_);
		solver2.analyzePattern(this->coeffmatrix_);
		solver2.factorize(this->coeffmatrix_);
		*/
		Eigen::SparseLU<SpMat> solver;
		solver.analyzePattern(coeffmatrix_);
		solver.factorize(coeffmatrix_);
		/*
		std::string name("Coefmatrix.csv");
		std::ofstream file(name.c_str());
		file<<MatrixXr(this->coeffmatrix_).format(CSVFormat);
		*/
		for(auto j=0;j<niter;j++)
		{	
			FPCAinput.setObservationData();

			
			
			VectorXr rightHandData;
			computeRightHandData(rightHandData,FPCAinput);
			b_ = VectorXr::Zero(2*nnodes);
			b_.topRows(nnodes)=rightHandData;
			
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
/*			
			std::chrono::high_resolution_clock::time_point t19= std::chrono::high_resolution_clock::now();
			
	this->solution_[i].resize(this->coeffmatrix_.rows());
		this->solution_[i]=solver.solve(this->b_);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
			
			std::chrono::high_resolution_clock::time_point t20= std::chrono::high_resolution_clock::now();
			
			std::chrono::duration<double> duration10 =t20-t19;
	
	std::cout<<"Time elapsed for solving sistem BiCGSTAB ILUT: "<<duration10.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t31= std::chrono::high_resolution_clock::now();
			
	this->solution_[i].resize(this->coeffmatrix_.rows());
		this->solution_[i]=solver2.solve(this->b_);
		if(solver2.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
			
			std::chrono::high_resolution_clock::time_point t32= std::chrono::high_resolution_clock::now();
			
			std::chrono::duration<double> duration16 =t32-t31;
	
	std::cout<<"Time elapsed for solving sistem BiCGSTAB: "<<duration16.count()<<std::endl;
	*/
	/*std::chrono::high_resolution_clock::time_point t33= std::chrono::high_resolution_clock::now();*/
			
		solution_[i].resize(coeffmatrix_.rows());
		solution_[i]=solver.solve(b_);
		if(solver.info()!=Eigen::Success)
		{
		//std::cerr<<"solving failed!"<<std::endl;
		}
			
		/*	std::chrono::high_resolution_clock::time_point t34= std::chrono::high_resolution_clock::now();
			
			std::chrono::duration<double> duration17 =t34-t33;
	
	std::cout<<"Time elapsed for solving sistem SpLU: "<<duration17.count()<<std::endl;*/
	/*
	std::chrono::high_resolution_clock::time_point t25= std::chrono::high_resolution_clock::now();
	
	
			this-> template solve<BiCGSTAB>(i);
			
			std::chrono::high_resolution_clock::time_point t26= std::chrono::high_resolution_clock::now();
			
			std::chrono::duration<double> duration13 =t26-t25;
	
	std::cout<<"Time elapsed for solving sistem BiCGSTAB: "<<duration13.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t15= std::chrono::high_resolution_clock::now();
	
	
			this-> template solve<SpLU>(i);
			
			std::chrono::high_resolution_clock::time_point t16= std::chrono::high_resolution_clock::now();
			
			std::chrono::duration<double> duration8 =t16-t15;
	
	std::cout<<"Time elapsed for solving sistem SpLU : "<<duration8.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t17= std::chrono::high_resolution_clock::now();
	
	
			this-> template solve<SpConjGrad>(i);
			
			std::chrono::high_resolution_clock::time_point t18= std::chrono::high_resolution_clock::now();
			
			std::chrono::duration<double> duration9 =t18-t17;
	
	std::cout<<"Time elapsed for solving sistem ConjGrad: "<<duration9.count()<<std::endl;
	*/
			
			if(fpcaData_.isLocationsByNodes())
				FPCAinput.setLoadings(nnodes, solution_[i]);
			else
				FPCAinput.setLoadingsPsi(nnodes, solution_[i],Psi_);
	
			
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
		/*if(this->inputData_.computeDOF())
			computeDegreesOfFreedom(i);
		else
			this->dof_[i] = -1;
			*/
	}
	
/*	std::chrono::high_resolution_clock::time_point t10= std::chrono::high_resolution_clock::now();
	
	std::chrono::duration<double> duration5 =t10-t9;
	
	std::cout<<"Time elapsed for computing a single NPC: "<<duration5.count()<<std::endl;*/
	
	scores_mat_[np]=FPCAinput.getScores();
	loadings_mat_[np]=FPCAinput.getLoadings();
	lambda_PC_[np]=fpcaData_.getLambda()[0];
	
	//Devo settare la datamatrix togliendo i risultati ottenuti
	FPCAinput.newDatamatrix(scores_mat_[np],loadings_mat_[np]);
	
	//Normalize the loadings and unnormalize the scores
	Real load_norm=std::sqrt(loadings_mat_[np].transpose()*MMat_*loadings_mat_[np]);
	
	loadings_mat_[np]=loadings_mat_[np]/load_norm;
	
	scores_mat_[np]=scores_mat_[np]*load_norm;
	
	//std::cout<<"Scores*load transp"<<std::endl;
	//std::cout<<scores_mat_[np]*loadings_mat_[np].transpose()<<std::endl;
	//std::cout<<"NewDataMAT"<<std::endl;
	//std::cout<<this->inputData_.getDatamatrix()-scores_mat_[np]*loadings_mat_[np].transpose()<<std::endl;
	}
	
	/*std::chrono::high_resolution_clock::time_point t12= std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration6 =t12-t11;
	
	std::cout<<"Time elapsed for computing all NPC: "<<duration6.count()<<std::endl;
	
	std::chrono::high_resolution_clock::time_point t13= std::chrono::high_resolution_clock::now();*/
	computeVarianceExplained();
	computeCumulativePercentageExplained();
/*	std::chrono::high_resolution_clock::time_point t14= std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration7 =t14-t13;
	
	std::cout<<"Time elapsed for computing variances: "<<duration7.count()<<std::endl;
*/
}

template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
template <typename P>
void MixedFEFPCABase<InputHandler,Integrator,ORDER,mydim,ndim>::solve(UInt output_index)
{
	//std::cout<<this->coeffmatrix_;
	this->solution_[output_index].resize(this->coeffmatrix_.rows());
	P::solve(this->coeffmatrix_,this->b_,this->solution_[output_index]);
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
