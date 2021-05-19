#ifndef TENTSOLVER_IMPL
#define TENTSOLVER_IMPL
#include "tentsolver.hpp"

template<typename T>
FlatVector<> AsFV (T & mat)
{
  return FlatVector<>( mat.Height()*mat.Width(), &mat(0,0) );
};

////// structure-aware Taylor time stepping //////
template <typename TCONSLAW>
void SAT<TCONSLAW>::PropagateTent(const Tent & tent, BaseVector & hu,
				  const BaseVector & hu0, LocalHeap & lh)
{
  static Timer tproptent ("SAT::Propagate Tent", 2);
  ThreadRegionTimer reg(tproptent, TaskManager::GetThreadId());

  tent.fedata = new (lh) TentDataFE(tent, *(tcl->fes), lh);
  tent.InitTent(tcl->gftau);

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
      for(int k : Range(stages))
  	{
  	  tcl->Cyl2Tent(tent, j*taustar, local_uhat1, local_u, lh);
  	  tcl->CalcFluxTent(tent, local_u, local_u0, local_uhat1, j*taustar, k, lh);
  	  local_uhat1 *= 1.0/(k+1);
  	  fac *= taustar;
  	  local_uhat += fac*local_uhat1;           
  
  	  if(k < stages-1)
  	    {
  	      tcl->ApplyM1(tent, j*taustar, local_u, local_help, lh);
  	      local_uhat1 += local_help;
  	    }
  	  local_u0 = 0.0;
  	}
    }
  hu.SetIndirect(tent.fedata->dofs, AsFV(local_uhat));
  tent.fedata = nullptr;
  tent.SetFinalTime();
};

////// structure-aware Runge-Kutta time stepping //////
template <typename TCONSLAW>
void SARK<TCONSLAW>::PropagateTent(const Tent & tent, BaseVector & hu,
				   const BaseVector & hu0, LocalHeap & lh)
{
  static Timer tproptent ("SARK::Propagate Tent", 2);
  ThreadRegionTimer reg(tproptent, TaskManager::GetThreadId());

  static Timer tres ("calc residual", 2);
  static Timer tnu ("calc nu", 2);
  static Timer tvisc ("apply viscosity", 2);

  tent.fedata = new (lh) TentDataFE(tent, *(tcl->fes), lh);
  tent.InitTent(tcl->gftau);

  const int ndof = tent.fedata->nd;
  FlatMatrixFixWidth<COMP> local_u0(ndof,lh);
  FlatMatrixFixWidth<COMP> local_Gu0(ndof,lh);
  FlatMatrixFixWidth<COMP> local_init(ndof,lh);

  hu.GetIndirect(tent.fedata->dofs, AsFV(local_Gu0));
  hu0.GetIndirect(tent.fedata->dofs, AsFV(local_init));

  FlatMatrixFixWidth<COMP> local_u1(ndof,lh);
  FlatMatrixFixWidth<COMP> local_u1_half(ndof,lh);
  FlatMatrixFixWidth<COMP> local_u(ndof,lh);
  FlatMatrixFixWidth<COMP> local_help(ndof,lh);
  FlatMatrixFixWidth<COMP> local_flux(ndof,lh), local_flux1(ndof,lh);

  FlatMatrixFixWidth<COMP> U0(ndof,lh);
  FlatMatrixFixWidth<COMP> U1(ndof,lh);
  FlatMatrixFixWidth<COMP> U2(ndof,lh);
  FlatMatrixFixWidth<COMP> U3(ndof,lh);
  FlatMatrixFixWidth<COMP> u0(ndof,lh);
  FlatMatrixFixWidth<COMP> u1(ndof,lh);
  FlatMatrixFixWidth<COMP> u2(ndof,lh);
  FlatMatrixFixWidth<COMP> u3(ndof,lh);
  FlatMatrixFixWidth<COMP> M1u0(ndof,lh);
  FlatMatrixFixWidth<COMP> M1u1(ndof,lh);
  FlatMatrixFixWidth<COMP> M1u2(ndof,lh);
  FlatMatrixFixWidth<COMP> M1u3(ndof,lh);
  FlatMatrixFixWidth<COMP> fu0(ndof,lh);
  FlatMatrixFixWidth<COMP> fu1(ndof,lh);
  FlatMatrixFixWidth<COMP> fu2(ndof,lh);
  FlatMatrixFixWidth<COMP> fu3(ndof,lh);
  FlatMatrixFixWidth<COMP> Uhat(ndof,lh);

  shared_ptr<BaseVector> hres = (ECOMP > 0) ? tcl->gfres->GetVectorPtr() : nullptr;
  FlatMatrixFixWidth<COMP> dUhatdt(ndof,lh);
  FlatMatrixFixWidth<ECOMP> res(ndof,lh);
  FlatVector<> local_nu(tent.els.Size(),lh);

  const double tau_tent = tent.ttop - tent.tbot;
  double h_tent = 0;
  for (int j : Range(tent.els))
    h_tent = max2(h_tent, tent.fedata->mesh_size[j]);

  const int order = max(1,tcl->fes->GetOrder());

  const double tau_visc1 = sqr (h_tent / sqr(order));

  // // calc |u|_M0 on advancing front
  // tcl->Cyl2Tent (tent, 0, local_Gu0, local_u0, lh);
  // tcl->Tent2Cyl(tent, 0, local_u0, local_help, false, lh);
  // double norm_bot = InnerProduct(AsFV(local_u0),AsFV(local_help));

  const double taustar = 1.0/substeps;
  //TODO: implement SARK for different number to stages
  for (int j = 0; j < substeps; j++)
    {
      // third order
      U0 = local_Gu0;
      tcl->Cyl2Tent (tent, j*taustar, U0, u0, lh);
      tcl->ApplyM1(tent, j*taustar, u0, M1u0, lh);
      tcl->CalcFluxTent(tent, u0, local_init, fu0, j*taustar, 0, lh);
      U1 = U0 + 0.5*taustar*(M1u0+fu0);

      tcl->Cyl2Tent (tent, j*taustar, U1, u1, lh);
      tcl->ApplyM1(tent, j*taustar, u1, M1u1, lh);
      tcl->CalcFluxTent(tent, u1, local_init, fu1, (j+0.5)*taustar, 0, lh);
      U2 = U0 + taustar*(4*M1u1-3*M1u0+2*fu1-fu0);

      tcl->Cyl2Tent (tent, j*taustar, U2, u2, lh);
      tcl->CalcFluxTent(tent, u2, local_init, fu2, (j+1)*taustar, 0, lh);

      Uhat = U0 + taustar * (1.0/6.0 * fu0 + 2.0/3.0 * fu1 + 1.0/6.0 * fu2);
      local_Gu0 = Uhat;
      // for viscosity
      dUhatdt = fu0;
      if (ECOMP > 0)
	{
	  /////// use dUhatdt as approximation at the final time
	  // U0 = local_Gu0;
	  // tcl->Cyl2Tent (tent, (j+1)*taustar, local_Gu0, local_help, lh);
	  // tcl->CalcFluxTent(tent, local_help, local_init, dUhatdt,
	  //              (j+1)*taustar, 0, lh);
	  // tcl->CalcEntropyResidualTent(tent, U0, dUhatdt, res, local_init,
	  //                         (j+1)*taustar, lh);
	  // hres->SetIndirect(tent.dofs,AsFV(res));
	  // double nu_tent = tcl->CalcViscosityCoefficientTent(
	  //                        tent, U0, res, (j+1)*taustar, lh);
	  /////// use dUhatdt as approximation at the initial time
	  tres.Start();
	  tcl->CalcEntropyResidualTent(tent, U0, dUhatdt, res, local_init, j*taustar, lh);
	  tres.Stop();
	  hres->SetIndirect(tent.fedata->dofs,AsFV(res));
	  tnu.Start();
	  double nu_tent = tcl->CalcViscosityCoefficientTent(tent, U0, res,j*taustar, lh);
	  tnu.Stop();
	  local_nu = nu_tent;
	  double steps_visc = (40*tau_tent*nu_tent/tau_visc1)/substeps;
	  tvisc.Start();
	  if (steps_visc > 0.2)
	    {
	      steps_visc = max(1.0,ceil(steps_visc));
	      double tau_visc = taustar/steps_visc;
	      // store boundary conditions in local_help
	      tcl->Cyl2Tent (tent, (j+1)*taustar, local_Gu0, local_u, lh);
	      local_help = local_u;
	      for (int k = 0; k < steps_visc; k++)
		{
		  tcl->CalcViscosityTent (tent, local_u, local_help, local_nu, local_flux, lh);
		  local_u -= tau_visc * local_flux;
		}
	      tcl->Tent2Cyl(tent, (j+1)*taustar, local_u, local_Gu0, true, lh);
	    }
	  tvisc.Stop();
	}
    }

  // // calc |u|_M1 norm on advancing front
  // tcl->Cyl2Tent (tent, 1, local_Gu0, local_u0, lh);
  // tcl->Tent2Cyl(tent, 1, local_u0, local_help, false, lh);
  // double norm_top = InnerProduct(AsFV(local_u0),AsFV(local_help));
  // if(norm_bot>0)
  //   *testout << "top/bot = " << norm_top/norm_bot << endl;
  // else
  //   *testout << "bot, top = " << norm_bot << ", " << norm_top << endl;
  // if(norm_top/norm_bot > 1.0)
  //   *testout << "bot, top : " << norm_bot << ", " << norm_top << endl;
  // if(norm_top/norm_bot < 0.9)
  //   *testout << "bot, top : " << norm_bot << ", " << norm_top << endl;

  hu.SetIndirect(tent.fedata->dofs, AsFV(local_Gu0));
  tent.fedata = nullptr;
  tent.SetFinalTime();
};


#endif //TENTSOLVER_IMPL
