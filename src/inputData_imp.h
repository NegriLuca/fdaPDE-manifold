#ifndef __INPUTDATA_IMP_HPP__
#define __INPUTDATA_IMP_HPP__

InputData::InputData(std::vector<Point>& locations, UInt order, std::vector<Real> lambda, bool DOF):
					locations_(locations), order_(order),lambda_(lambda),DOF_(DOF)
{
	if(locations_.size()==0)
		locations_by_nodes_= true;
	else
		locations_by_nodes_= false;
}

#ifdef R_VERSION_
InputData::InputData(SEXP Rlocations, SEXP Rorder, SEXP Rlambda, SEXP DOF)
{
	setLocations(Rlocations);

	order_ =  INTEGER(Rorder)[0];
	DOF_ = INTEGER(DOF)[0];

    UInt length_lambda = Rf_length(Rlambda);
    for (UInt i = 0; i<length_lambda; ++i)  lambda_.push_back(REAL(Rlambda)[i]);

}

void InputData::setLocations(SEXP Rlocations)
{
	UInt n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0){
		locations_by_nodes_=false;
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
	}else
	locations_by_nodes_=true;
}

#endif


void InputData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

#endif
