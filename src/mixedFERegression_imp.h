#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "R_ext/Print.h"

//#include <libseq/mpi.h>
#include "../inst/include/dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::addDirichletBC()
{
	UInt id1,id3;

	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bc_indices[i];
			id3=id1+nnodes;

			A_.coeffRef(id1,id1)=pen;
			A_.coeffRef(id3,id3)=pen;


			_b(id1)+=bc_values[i]*pen;
			_b(id3)=0;
	 }

	A_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setPsi(){

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	psi_.resize(nlocations, nnodes);
	if (regressionData_.isLocationsByNodes()){

		std::vector<coeff> tripletAll;
		auto k = regressionData_.getObservationsIndices();
		tripletAll.reserve(k.size());
		for (int i = 0; i< k.size(); ++i){
			tripletAll.push_back(coeff(i,k[i],1.0));
		}
		psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
		psi_.makeCompressed();
	}
	else {
		Triangle<3*ORDER+mydim%2, mydim, ndim> tri_activated;
		Eigen::Matrix<Real,3*ORDER+mydim%2,1> coefficients;

		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else
			{
				for(UInt node = 0; node < 3*ORDER+mydim%2 ; ++node)
				{
					coefficients = Eigen::Matrix<Real,3*ORDER+mydim%2,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<ORDER>(tri_activated, regressionData_.getLocations()[i], coefficients);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}

		psi_.makeCompressed();
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::LeftMultiplybyQ(const MatrixXr& u)
{	
	if (regressionData_.getCovariates().rows() == 0){
		return u;
	}
	else{
		MatrixXr W(this->regressionData_.getCovariates());
		if (isWTWfactorized_ == false ){
			WTWinv_.compute(W.transpose()*W);
			isWTWfactorized_=true;
		}
		MatrixXr Pu= W*WTWinv_.solve(W.transpose()*u);
		return u-Pu;
	}

}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildA(const SpMat& Psi,  const SpMat& R1,  const SpMat& R0) {

	UInt nnodes = mesh_.num_nodes();
	
	SpMat DMat = Psi.transpose()*Psi;

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*R1.nonZeros() + R0.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
		for (SpMat::InnerIterator it(DMat,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
		}
	for (int k=0; k<R0.outerSize(); ++k)
		for (SpMat::InnerIterator it(R0,k); it; ++it)
		{
			tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
		}
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<R1.outerSize(); ++k)
	  for (SpMat::InnerIterator it(R1,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	A_.setZero();
	A_.resize(2*nnodes,2*nnodes);
	A_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	A_.makeCompressed();
	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::system_factorize() {

	UInt nnodes = mesh_.num_nodes();

	// First phase: Factorization of matrix A
	Adec_.compute(A_);
	
	if (regressionData_.getCovariates().rows() != 0) {
		// Second phase: factorization of matrix  G =  C + [V * A^-1 * U]

		// Definition of matrix U = [ psi * W | 0 ]^T
		MatrixXr W(this->regressionData_.getCovariates());
		U_ = MatrixXr::Zero(2*nnodes, W.cols());
		U_.topRows(nnodes) = psi_.transpose()*W;

		// D = U^T * A^-1 * U
		MatrixXr D = U_.transpose()*Adec_.solve(U_);
		// G = C + D
		MatrixXr G = -W.transpose()*W + D;
		Gdec_.compute(G);
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename Derived>
MatrixXr MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::system_solve(const Eigen::MatrixBase<Derived> &b) {

	// Resolution of the system A * x1 = b
	MatrixXr x1 = Adec_.solve(b);
	
	if (regressionData_.getCovariates().rows() != 0) {
		// Resolution of G * x2 = U^T * x1
		MatrixXr x2 = Gdec_.solve(U_.transpose()*x1);
		// Resolution of the system A * x3 = U * x2
		x1 -= Adec_.solve(U_*x2);
	}

	return x1;
}


 template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::setQ()
 {
 	//std::cout<<"Computing Orthogonal Space Projection Matrix"<<std::endl;
 	Q_.resize(H_.rows(),H_.cols());
 	Q_ = -H_;
 	for (int i=0; i<H_.rows();++i)
 	{
 		Q_(i,i) += 1;
 	}
 }
 
 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim, ndim>::setH()
 {
 	//std::cout<<"Computing Projection Matrix"<<std::endl;
 	//UInt nnodes = mesh_.num_nodes();
 	UInt nlocations = regressionData_.getNumberofObservations();
 
 	//regressionData_.printCovariates(std::cout);
 	MatrixXr W(this->regressionData_.getCovariates());
 	//std::cout<<"W "<< W <<std::endl;
 	//total number of mesh nodes
 	//UInt nnodes = mesh_.num_nodes();
 	if(regressionData_.isLocationsByNodes())
 	{
 		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
 		for (auto i=0; i<nlocations;++i)
 		{
 			auto index_i = regressionData_.getObservationsIndices()[i];
 			for (auto j=0; j<W.cols();++j)
 			{
 				W_reduced(i,j) = W(index_i,j);
 			}
 		}
 		W = W_reduced;
 	}
 
 
 	MatrixXr WTW(W.transpose()*W);
 
 	H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix
 }


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
 {
 	//I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
 	//_coeffmatrix.setFromTriplets(tripletA.begin(),tripletA.end());
 
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
 
 	_coeffmatrix.setZero();
 	_coeffmatrix.resize(2*nnodes,2*nnodes);
 	_coeffmatrix.setFromTriplets(tripletAll.begin(),tripletAll.end());
 	_coeffmatrix.makeCompressed();
 	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
 }

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::getDataMatrix(SpMat& DMat)
 {
 		UInt nnodes = mesh_.num_nodes();
 		//UInt nlocations = regressionData_.getNumberofObservations();
 
 		DMat.resize(nnodes,nnodes);
 
 		if (regressionData_.getCovariates().rows() == 0)
 			DMat = psi_.transpose()*psi_;
 		else
 		{
 			DMat = (SpMat(psi_.transpose())*Q_*psi_).sparseView();
 		}
 }
 
 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::getDataMatrixByIndices(SpMat& DMat)
 {
 		UInt nnodes = mesh_.num_nodes();
 		UInt nlocations = regressionData_.getNumberofObservations();
 
 		DMat.resize(nnodes,nnodes);
 
 		if (regressionData_.getCovariates().rows() == 0)
 		{
 			DMat.reserve(1);
 			for (auto i = 0; i<nlocations; ++i)
 			{
 				auto index = regressionData_.getObservationsIndices()[i];
 				DMat.insert(index,index) = 1;
 			}
 		}
 		else
 		{
 			//May be inefficient
 			for (auto i = 0; i<nlocations; ++i)
 			{
 				auto index_i = regressionData_.getObservationsIndices()[i];
 				for (auto j = 0; j<nlocations; ++j)
 				{
 					auto index_j = regressionData_.getObservationsIndices()[j];
 					DMat.insert(index_i,index_j) = Q_(i,j);
 				}
 			}
 		}
 }

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	rightHandData = VectorXr::Zero(nnodes);

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData=psi_.transpose()*regressionData_.getObservations();
		}
	}
	else
	{
		rightHandData=psi_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedom(UInt output_index, Real lambda)
{
	int GCVmethod = regressionData_.getGCVmethod();
	switch (GCVmethod) {
		case 1:
			computeDegreesOfFreedomExact(output_index, lambda);
			break;
		case 2:
			computeDegreesOfFreedomStochastic(output_index, lambda);
			break;
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedomExact(UInt output_index, Real lambda)
{

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	Real degrees=0;

	// Case 1: MUMPS
	if (regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() == 0 )
	{
		auto k = regressionData_.getObservationsIndices();
		DMUMPS_STRUC_C id;
		int myid, ierr;
        int argc=0;
        char ** argv= NULL;
        //MPI_Init(&argc,&argv);
		//ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

		id.sym=0;
		id.par=1;
		id.job=JOB_INIT;
		id.comm_fortran=USE_COMM_WORLD;
		dmumps_c(&id);

		std::vector<int> irn;
		std::vector<int> jcn;
		std::vector<double> a;
		std::vector<int> irhs_ptr;
		std::vector<int> irhs_sparse;
		double* rhs_sparse= (double*)malloc(nlocations*sizeof(double));
		
		//if( myid==0){
			id.n=2*nnodes;
			for (int j=0; j<A_.outerSize(); ++j){
				for (SpMat::InnerIterator it(A_,j); it; ++it){
					irn.push_back(it.row()+1);
					jcn.push_back(it.col()+1);
					a.push_back(it.value());
				}
			}
		//}
		id.nz=irn.size();
		id.irn=irn.data();
		id.jcn=jcn.data();
		id.a=a.data();
		id.nz_rhs=nlocations;
		id.nrhs=2*nnodes;
		int j = 1;
		irhs_ptr.push_back(j);
		for (int l=0; l<k[0]-1; ++l) {
			irhs_ptr.push_back(j);
		}
		for (int i=0; i<k.size()-1; ++i) {
			++j;
			for (int l=0; l<k[i+1]-k[i]; ++l) {
				irhs_ptr.push_back(j);
			}
			
		}
		++j;
		for (int i=k[k.size()-1]; i < id.nrhs; ++i) {
			irhs_ptr.push_back(j);
		}
		for (int i=0; i<nlocations; ++i){
			irhs_sparse.push_back(k[i]+1);
		}
		id.irhs_sparse=irhs_sparse.data();
		id.irhs_ptr=irhs_ptr.data();
		id.rhs_sparse=rhs_sparse;

		#define ICNTL(I) icntl[(I)-1]
		//Output messages suppressed
		id.ICNTL(1)=-1;
		id.ICNTL(2)=-1;
		id.ICNTL(3)=-1;
		id.ICNTL(4)=0;
		id.ICNTL(20)=1;
		id.ICNTL(30)=1;
		id.ICNTL(14)=200;

		id.job=6;
		dmumps_c(&id);
		id.job=JOB_END;
		dmumps_c(&id);

		//if (myid==0){
			for (int i=0; i< nlocations; ++i){
				//std::cout << "rhs_sparse" << rhs_sparse[i] << std::endl;
				degrees+=rhs_sparse[i];
			}
		//}
		free(rhs_sparse);

		//MPI_Finalize();
	}
	// Case 2: Eigen
	else{
		MatrixXr X1 = psi_.transpose() * LeftMultiplybyQ(psi_);

		if (isRcomputed_ == false ){
			isRcomputed_ = true;
			Eigen::SparseLU<SpMat> solver;
			solver.compute(R0_);
			auto X2 = solver.solve(R1_);
			R_ = R1_.transpose() * X2;
		}

		MatrixXr X3 = X1 + lambda * R_;
		Eigen::LDLT<MatrixXr> Dsolver(X3);

		auto k = regressionData_.getObservationsIndices();

		if(regressionData_.isLocationsByNodes() && regressionData_.getCovariates().rows() != 0) {
			degrees += regressionData_.getCovariates().cols();

			// Setup rhs B
			MatrixXr B;
			B = MatrixXr::Zero(nnodes,nlocations);
			// B = I(:,k) * Q
			for (auto i=0; i<nlocations;++i) {
				VectorXr ei = VectorXr::Zero(nlocations);
				ei(i) = 1;
				VectorXr Qi = LeftMultiplybyQ(ei);
				for (int j=0; j<nlocations; ++j) {
					B(k[i], j) = Qi(j);
				}
			}
			// Solve the system TX = B
			MatrixXr X;
			X = Dsolver.solve(B);
			// Compute trace(X(k,:))
			for (int i = 0; i < k.size(); ++i) {
				degrees += X(k[i], i);
			}
		}

		if (!regressionData_.isLocationsByNodes()){
			MatrixXr X;
			X = Dsolver.solve(MatrixXr(X1));

			if (regressionData_.getCovariates().rows() != 0) {
				degrees += regressionData_.getCovariates().cols();
			}
			for (int i = 0; i<nnodes; ++i) {
				degrees += X(i,i);
			}
		}
	}
	_dof[output_index] = degrees;
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::computeDegreesOfFreedomStochastic(UInt output_index, Real lambda)
{	
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	std::default_random_engine generator;
	// Creation of the random matrix
	std::bernoulli_distribution distribution(0.5);
	UInt nrealizations = regressionData_.getNrealizations();
	MatrixXr u(nlocations, nrealizations);
	for (int j=0; j<nrealizations; ++j) {
		for (int i=0; i<nlocations; ++i) {
			if (distribution(generator)) {
				u(i,j) = 1.0;
			}
			else {
				u(i,j) = -1.0;
			}
		}
	}

	// Define the first right hand side : | I  0 |^T * psi^T * Q * u
	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	b.topRows(nnodes) = psi_.transpose()* LeftMultiplybyQ(u);

	// Resolution of the system
	//MatrixXr x = system_solve(b);
	Eigen::SparseLU<SpMat> solver;
	solver.compute(_coeffmatrix);
	auto x = solver.solve(b);

	MatrixXr uTpsi = u.transpose()*psi_;
	VectorXr edf_vect(nrealizations);
	Real q = 0;

	// Degrees of freedom = q + E[ u^T * psi * | I  0 |* x ]
	if (regressionData_.getCovariates().rows() != 0) {
		q = regressionData_.getCovariates().cols();
	}
	// For any realization we calculate the degrees of freedom
	for (int i=0; i<nrealizations; ++i) {
		edf_vect(i) = uTpsi.row(i).dot(x.col(i).head(nnodes)) + q;
	}

	// Estimates: sample mean, sample variance
	Real mean = edf_vect.sum()/nrealizations;
	_dof[output_index] = mean;
}


/*template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeBasisEvaluations(){

	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	Real eps = 2.2204e-016,
		 tolerance = 100 * eps;

	Psi_.resize(nlocations, nnodes);
	//Psi_.reserve(Eigen::VectorXi::Constant(nlocations,ORDER*3));

	Triangle<3*ORDER+mydim%2,mydim, ndim> tri_activated;
	Eigen::Matrix<Real,3*ORDER+mydim%2,1> coefficients;

	Real evaluator;
	for(UInt i=0; i<nlocations;i++)
	{
		tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
		if(tri_activated.getId() == Identifier::NVAL)
		{
			#ifdef R_VERSION_
			Rprintf("WARNING: Observation %d is not in the domain, remove point and re-perform smoothing\n", i+1);
			#else
			std::cout << "WARNING: Observation " << i+1 <<" is not in the domain\n";
			#endif
		}else
		{
			for(UInt node = 0; node < 3*ORDER+mydim%2 ; ++node)
			{
				coefficients = Eigen::Matrix<Real,3*ORDER+mydim%2,1>::Zero();
				coefficients(node) = 1; //Activates only current base
				evaluator = evaluate_point<ORDER, mydim, ndim>(tri_activated, regressionData_.getLocations()[i], coefficients);
				Psi_.insert(i, tri_activated[node].getId()) = evaluator;
			}
		}
	}

	Psi_.prune(tolerance);
	Psi_.makeCompressed();
}

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeProjOnCovMatrix()
{
	//std::cout<<"Computing Projection Matrix"<<std::endl;
	UInt nlocations = regressionData_.getNumberofObservations();

	//Set covariate matrix
	MatrixXr W(regressionData_.getCovariates());
	if(regressionData_.isLocationsByNodes())
	{
		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = regressionData_.getObservationsIndices()[i];
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
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
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


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::applyDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues)
{
	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	UInt id1,id3;
	for( auto i=0; i<nbc_indices; i++)
	 {
		id1=bcindex[i];
		id3=id1+nnodes;

		coeffmatrix_.coeffRef(id1,id1) = dirichlet_penalization;
		coeffmatrix_.coeffRef(id3,id3) = dirichlet_penalization;


		b_(id1) += bc_values[i]*dirichlet_penalization;
		b_(id3) = 0;
	 }
	coeffmatrix_.makeCompressed();
}

//construct NW block of the system matrix when basis evaluation is necessary
//!! Depends on computeBasisEvaluations and computeProjOnCovMatrix
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrix(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		DMat.resize(nnodes,nnodes);

		if (regressionData_.getCovariates().rows() == 0)
			DMat = Psi_.transpose()*Psi_;
		else
			DMat = (SpMat(Psi_.transpose())*P_*Psi_).sparseView();
}

//construct NW block of the system matrix in Ramsay when locations of observations are
//a subset of the meshe's nodes
//!! Depends on computeBasisEvaluations and computeProjOnCovMatrix
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		UInt nlocations = regressionData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		if (regressionData_.getCovariates().rows() == 0)
		{
			DMat.reserve(1);
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index = regressionData_.getObservationsIndices()[i];
				DMat.insert(index,index) = 1;
			}
		}
		else
		{
			//May be inefficient
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				for (auto j = 0; j<nlocations; ++j)
				{
					auto index_j = regressionData_.getObservationsIndices()[j];
					DMat.insert(index_i,index_j) = P_(i,j);
				}
			}
		}
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	//rightHandData.resize(nnodes);
	rightHandData = VectorXr::Zero(nnodes);

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData=Psi_.transpose()*regressionData_.getObservations();
		}
	}
	else
	{
		if(regressionData_.isLocationsByNodes())
		{
			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = (P_).row(i) * regressionData_.getObservations();
			}
		}
		else
		{
			rightHandData=Psi_.transpose()*P_*regressionData_.getObservations();
		}
	}
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::computeDegreesOfFreedom(UInt output_index)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	Eigen::SparseLU<SpMat> solver;
	solver.compute(coeffmatrix_);
	SpMat I(coeffmatrix_.rows(),coeffmatrix_.cols());
	I.setIdentity();
	SpMat coeff_inv = solver.solve(I);


	Real degrees=0;

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{
			VectorXr d = coeff_inv.diagonal();

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
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
	}
	else
	{	//Important
		degrees = regressionData_.getCovariates().cols();

		if(regressionData_.isLocationsByNodes())
		{
			MatrixXr S(coeff_inv);

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				for (auto j=0; j<nlocations;++j)
				{
					auto index_j = regressionData_.getObservationsIndices()[j];

					degrees+=S(index_i,index_j)*P_(j,i);
				}
			}

		}
		else
		{
			MatrixXr An(coeff_inv.topLeftCorner(nnodes, nnodes));
			MatrixXr S = Psi_*An*Psi_.transpose()*P_;
			for (auto i=0; i<nlocations;++i)
			{
				degrees+=S(i,i);
			}
		}
	}


	//std::cout<<"TRACE "<<degrees<<std::endl;

	dof_[output_index] = degrees;
}
*/

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
template<typename A>
void MixedFERegressionBase<InputHandler,Integrator,ORDER, mydim, ndim>::apply(EOExpr<A> oper)
{
	UInt nnodes=mesh_.num_nodes();
	FiniteElement<Integrator, ORDER, mydim, ndim> fe;
	
	setPsi();

	if(!regressionData_.getCovariates().rows() == 0)
 	{
 		setH();
 		setQ();
 	}
 
	if(!regressionData_.isLocationsByNodes())
	{
		getDataMatrix(DMat_);
	}
	else
	{
		getDataMatrixByIndices(DMat_);
	}

	typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
	Assembler::operKernel(oper, mesh_, fe, R1_);
	Assembler::operKernel(mass, mesh_, fe, R0_);
	VectorXr rightHandData;
	getRightHandData(rightHandData);
	this->_b = VectorXr::Zero(2*nnodes);
	this->_b.topRows(nnodes)=rightHandData;
	this->_solution.resize(regressionData_.getLambda().size());
	this->_dof.resize(regressionData_.getLambda().size());

	for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
		Real lambda = regressionData_.getLambda()[i];
		SpMat R1_lambda = (-lambda)*R1_;
		SpMat R0_lambda = (-lambda)*R0_;
		this->buildA(psi_, R1_lambda, R0_lambda);
		this->buildCoeffMatrix(DMat_, R1_lambda, R0_lambda);
		//std::cout << coeffmatrix_ << std::endl;
		//Applying boundary conditions if necessary
		if(regressionData_.getDirichletIndices().size() != 0)
			addDirichletBC();

		system_factorize();
		_solution[i] = this->template system_solve(this->_b);
		if(regressionData_.computeDOF())
			computeDegreesOfFreedom(i,lambda);
		else
			_dof[i] = -1;
	}
}

/*template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
template <typename P>
void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim,ndim>::solve(UInt output_index)
{
	//std::cout<<this->coeffmatrix_;
	this->solution_[output_index].resize(this->coeffmatrix_.rows());
	P::solve(this->coeffmatrix_,this->b_,this->solution_[output_index]);
}
*/


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

	    const Real& c = this->regressionData_.getC();
	    const Eigen::Matrix<Real,2,2>& K = this->regressionData_.getK();
	    const Eigen::Matrix<Real,2,1>& beta = this->regressionData_.getBeta();

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

		const Reaction& c = this->regressionData_.getC();
		const Diffusivity& K = this->regressionData_.getK();
		const Advection& beta = this->regressionData_.getBeta();

		MixedFERegressionBase<RegressionDataEllipticSpaceVarying, Integrator, ORDER, mydim, ndim>::apply(c*mass+stiff[K]+dot(beta,grad));
	}
	}
};

#endif
