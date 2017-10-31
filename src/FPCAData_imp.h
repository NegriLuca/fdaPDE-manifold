#ifndef __FPCADATA_IMP_HPP__
#define __FPCADATA_IMP_HPP__

FPCAData::FPCAData(std::vector<Point>& locations, MatrixXr& datamatrix, UInt order, std::vector<Real> lambda ,UInt nPC, bool DOF):
					locations_(locations), datamatrix_(datamatrix), order_(order), lambda_(lambda), nPC_(nPC), DOF_(DOF)
{
	if(locations_.size()==0)
	{
		locations_by_nodes_= true;
		for(int i = 0; i<datamatrix_.rows();++i) datamatrix_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_= false;
	}
}

#ifdef R_VERSION_
FPCAData::FPCAData(SEXP Rlocations, SEXP Rdatamatrix, SEXP Rorder, SEXP Rlambda, SEXP RnPC, SEXP DOF)
{
	setLocations(Rlocations);
	setDatamatrix(Rdatamatrix);

	order_ =  INTEGER(Rorder)[0];
	DOF_ = INTEGER(DOF)[0];
	nPC_ = INTEGER(RnPC)[0];
	
    	UInt length_lambda = Rf_length(Rlambda);
    	for (UInt i = 0; i<length_lambda; ++i)  lambda_.push_back(REAL(Rlambda)[i]);

}


void FPCAData::setDatamatrix(SEXP Rdatamatrix)
{
	n_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[1];
	datamatrix_.resize(n_,p_);
	datamatrix_indices_.reserve(n_);
	
	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			datamatrix_(i,j)=REAL(Rdatamatrix)[i+ n_*j];
		}
	}

	if(locations_.size() == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_;++i) datamatrix_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}

}

void FPCAData::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0){
		int ndim = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[1];

	  if (ndim == 2){
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
			}
		}else{
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1],REAL(Rlocations)[i+ n_*2]);
			}
		}
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


void FPCAData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

#endif

