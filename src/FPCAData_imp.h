#ifndef __FPCADATA_IMP_HPP__
#define __FPCADATA_IMP_HPP__

FPCAData::FPCAData(std::vector<Point>& locations, MatrixXr& datamatrix, UInt order, std::vector<Real> lambda ,UInt nPC, UInt nFolds):locations_(locations), order_(order),lambda_(lambda),
 datamatrix_(datamatrix), nPC_(nPC),nFolds_(nFolds)
{
	if(locations.size()==0)
	{
		locations_by_nodes_= true;
		for(int i = 0; i<datamatrix_.cols();++i) observations_indices_.push_back(i);
	} else
		locations_by_nodes_= false;
}

#ifdef R_VERSION_
FPCAData::FPCAData(SEXP Rlocations, SEXP Rdatamatrix, SEXP Rorder, SEXP Rlambda, SEXP RnPC, SEXP RnFolds,SEXP RGCVmethod, SEXP Rnrealizations, SEXP RRNGstate, SEXP Rsolver, SEXP Rnprocessors, SEXP Rhosts)
{
	
	setLocations(Rlocations);
	setDatamatrix(Rdatamatrix);

	setNrealizations(Rnrealizations);
	
	setSolver(Rsolver);
	setnprocessors(Rnprocessors);
	setRNGstate(RRNGstate);
	setHosts(Rhosts);
	
	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	
    	UInt length_lambda = Rf_length(Rlambda);
    	for (UInt i = 0; i<length_lambda; ++i)  lambda_.push_back(REAL(Rlambda)[i]);

	nPC_ = INTEGER(RnPC)[0];
	
	nFolds_=INTEGER(RnFolds)[0];
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


void FPCAData::setDatamatrix(SEXP Rdatamatrix)
{
	n_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rdatamatrix, R_DimSymbol))[1];
	datamatrix_.resize(n_,p_);
	observations_indices_.reserve(p_);
	VectorXr auxiliary_row_;
	auxiliary_row_.resize(p_);
	if(locations_.size() == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0; i<n_; ++i)
		{	
			UInt count=0;
			for(auto j=0; j<p_ ; ++j)
			{
				if(!ISNA(REAL(Rdatamatrix)[i+n_*j]))
				{
					auxiliary_row_[count]=REAL(Rdatamatrix)[i+n_*j];
					count++;
					if(i==0) observations_indices_.push_back(j);
				}
			}
			datamatrix_.row(i)=auxiliary_row_;
		}
		datamatrix_.conservativeResize(Eigen::NoChange,observations_indices_.size());
	} else {
		locations_by_nodes_ = false;
		for(auto i=0; i<n_; ++i)
		{
			for(auto j=0; j<p_ ; ++j)
			{
				datamatrix_(i,j)=REAL(Rdatamatrix)[i+n_*j];
			}
		}
	}

}



void FPCAData::setRNGstate(SEXP RRNGstate) {
	RNGstate_.assign(CHAR(STRING_ELT(RRNGstate, 0)));
}

void FPCAData::setNrealizations(SEXP Rnrealizations) {
	nrealizations_ = INTEGER(Rnrealizations)[0];
}

void FPCAData::setnprocessors(SEXP Rnprocessors) {
	nprocessors_ = INTEGER(Rnprocessors)[0];
}
void FPCAData::setSolver(SEXP Rsolver) {
	solver_.assign(CHAR(STRING_ELT(Rsolver, 0)));
}

void FPCAData::setHosts(SEXP Rhosts) {
	hosts_.assign(CHAR(STRING_ELT(Rhosts, 0)));
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
void FPCAData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}



#endif

