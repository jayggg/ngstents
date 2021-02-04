#ifndef TENTSOLVER_IMPL
#define TENTSOLVER_IMPL
#include "tentsolver.hpp"

template<typename T>
FlatVector<> AsFV (T & mat)
{
  return FlatVector<>( mat.Height()*mat.Width(), &mat(0,0) );
};

template <typename T>
void SAT<T>::PropagateTent(const int tentnr, BaseVector & hu,
			   const BaseVector & hu0, LocalHeap & lh)
{
  static Timer tproptent ("Propagate Tent SAT", 2);
  ThreadRegionTimer reg(tproptent, TaskManager::GetThreadId());

  const Tent & tent = tcl->tps->GetTent(tentnr);
  tent.fedata = new (lh) TentDataFE(tent, *(tcl->fes), *(tcl->ma), lh);

  int ndof = tent.fedata->nd;
  FlatMatrixFixWidth<COMP> local_uhat(ndof,lh);
  FlatMatrixFixWidth<COMP> local_u0(ndof,lh);
  FlatMatrixFixWidth<COMP> local_u0temp(ndof,lh);
  hu.GetIndirect(tent.fedata->dofs, AsFV(local_uhat));
  hu0.GetIndirect(tent.fedata->dofs, AsFV(local_u0));
  local_u0temp = local_u0;
  
  FlatMatrixFixWidth<COMP> local_uhat1(ndof,lh);
  FlatMatrixFixWidth<COMP> local_u(ndof,lh);
  FlatMatrixFixWidth<COMP> local_help(ndof,lh);
  
  double taustar = 1.0/substeps;
  for (int j = 0; j < substeps; j++)
    {
      local_uhat1 = local_uhat;
      double fac = 1.0;
      local_u0 = local_u0temp;
      for(int k : Range(1,stages+1))
  	{
  	  tcl->Cyl2Tent(tentnr, j*taustar, local_uhat1, local_u, lh);
  	  tcl->CalcFluxTent(tentnr, local_u, local_u0, local_uhat1, j*taustar, lh);
  	  local_uhat1 *= 1.0/k;
  	  fac *= taustar;
  	  local_uhat += fac*local_uhat1;           
  
  	  if(k < stages)
  	    {
  	      tcl->ApplyM1(tentnr, j*taustar, local_u, local_help, lh);
  	      local_uhat1 += local_help;
  	    }
  	  local_u0 = 0.0;
  	}
    }
  hu.SetIndirect(tent.fedata->dofs, AsFV(local_uhat));
  tent.fedata = nullptr; 
};

#endif //TENTSOLVER_IMPL
