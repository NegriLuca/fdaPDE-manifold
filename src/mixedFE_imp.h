#ifndef __MIXEDFE_IMP_HPP__
#define __MIXEDFE_IMP_HPP__

#include <iostream>

template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFE<InputHandler,Integrator,ORDER, mydim, ndim>::computeBasisEvaluations(){

	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = inputData_.getNumberofObservations();
	Real eps = 2.2204e-016,
		 tolerance = 100 * eps;

	Psi_.resize(nlocations, nnodes);
	//Psi_.reserve(Eigen::VectorXi::Constant(nlocations,ORDER*3));

	Triangle<ORDER*3, mydim, ndim> tri_activated;
	Eigen::Matrix<Real,ORDER * 3,1> coefficients;

	Real evaluator;
	for(UInt i=0; i<nlocations;i++)
	{
		tri_activated = mesh_.findLocationNaive(inputData_.getLocations()[i]);
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
				evaluator = evaluate_point<ORDER, mydim, ndim>(tri_activated, inputData_.getLocations()[i], coefficients);
				Psi_.insert(i, tri_activated[node].getId()) = evaluator;
			}
		}
	}

	Psi_.prune(tolerance);
	Psi_.makeCompressed();
}


template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void MixedFE<InputHandler,Integrator,ORDER, mydim, ndim>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
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

//solve sparse system with P method

template<typename InputHandler, typename Integrator, UInt ORDER,UInt mydim, UInt ndim>
template <typename P>
void MixedFE<InputHandler,Integrator,ORDER,mydim,ndim>::solve(UInt output_index)
{
	//std::cout<<this->coeffmatrix_;
	this->solution_[output_index].resize(this->coeffmatrix_.rows());
	P::solve(this->coeffmatrix_,this->b_,this->solution_[output_index]);
}

#endif
