#ifndef __MIXEDFEFPCAFACTORY_HPP__
#define __MIXEDFEFPCAFACTORY_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "solver.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "mixedFEFPCA.h"

#include <memory>

//! A LinearSystem class: A class for the linear system construction and resolution.
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCAfactory
{	
	public:
	static std::unique_ptr<MixedFEFPCABase<Integrator, ORDER,  mydim,  ndim>> createFPCAsolver(const std::string &validation, const MeshHandler<ORDER,mydim,ndim>& mesh, const FPCAData& fpcaData){
	if(validation=="GCV") return make_unique<MixedFEFPCAGCV<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	if(validation=="KFold") return make_unique<MixedFEFPCAKFold<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	if(validation=="NoValidation") return make_unique<MixedFEFPCA<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	}

};

#endif
