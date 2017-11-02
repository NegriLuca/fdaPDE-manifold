#ifndef __FPCADATA_IMP_HPP__
#define __FPCADATA_IMP_HPP__

FPCAData::FPCAData(std::vector<Point>& locations, MatrixXr& datamatrix, UInt order, std::vector<Real> lambda ,UInt nPC, bool DOF):InputData(locations, order, lambda,DOF),
 datamatrix_(datamatrix), nPC_(nPC)
{
	if(locations.size()==0)
	{
		for(int i = 0; i<datamatrix_.cols();++i) observations_indices_.push_back(i);
	}
}

#ifdef R_VERSION_
FPCAData::FPCAData(SEXP Rlocations, SEXP Rdatamatrix, SEXP Rorder, SEXP Rlambda, SEXP RnPC, SEXP DOF):InputData(Rlocations, Rorder, Rlambda, DOF)
{
	setDatamatrix(Rdatamatrix);

	nPC_ = INTEGER(RnPC)[0];
}


void FPCAData::setDatamatrix(SEXP Rdatamatrix)
{
	n_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[1];
	datamatrix_.resize(n_,p_);
	observations_indices_.reserve(p_);
	
	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			datamatrix_(i,j)=REAL(Rdatamatrix)[i+ n_*j];
		}
	}

	if(this->getLocations().size() == 0)
	{
		for(auto i=0;i<p_;++i) observations_indices_.push_back(i);
	}

}

#endif


void FPCAData::printDatamatrix(std::ostream & out) const
{

	for(auto i=0;i<datamatrix_.rows(); i++)
	{
		for(auto j=0;j<datamatrix_.cols();j++)
		{
		out<<datamatrix_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

/*void newDatamatrix(const VectorXr& scores_,const VectorXr& loadings_)
{    	datamatrix_=scores_*loadings_.transpose();
}
*/	


#endif

