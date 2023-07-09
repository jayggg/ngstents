#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

typedef CoefficientFunction CF;

template <int D, int COMP, int ECOMP>
class SymbolicConsLaw :
  public T_ConservationLaw<SymbolicConsLaw<D,COMP,ECOMP>, D, COMP, ECOMP, true>
{
  typedef T_ConservationLaw<SymbolicConsLaw<D, COMP, ECOMP>, D, COMP, ECOMP, true>
  BASE;

  shared_ptr<CF> cf_flux = nullptr;
  shared_ptr<CF> cf_numflux = nullptr;
  shared_ptr<CF> cf_invmap = nullptr;

  // cf's for entropy residual
  shared_ptr<CF> cf_entropy = nullptr;
  shared_ptr<CF> cf_entropyflux = nullptr;
  shared_ptr<CF> cf_numentropyflux = nullptr;
  shared_ptr<CF> cf_visccoeff = nullptr;

  // compiled differentials
  shared_ptr<CF> ddu_invmap = nullptr;
  shared_ptr<CF> ddphi_invmap = nullptr;
  shared_ptr<CF> ddu_entropy = nullptr;

  using BASE::proxy_u;
  using BASE::proxy_uother;
  using BASE::proxy_graddelta;
  using BASE::proxy_res;
  using BASE::tps;
public:
  SymbolicConsLaw (const shared_ptr<GridFunction> & agfu,
		   const shared_ptr<TentPitchedSlab> & atps,
		   const shared_ptr<ProxyFunction> & aproxy_u,
		   const shared_ptr<ProxyFunction> & aproxy_uother,
		   const shared_ptr<CF> & acf_flux,
		   const shared_ptr<CF> & acf_numflux,
		   const shared_ptr<CF> & acf_invmap,
		   const shared_ptr<CF> & acf_entropy,
		   const shared_ptr<CF> & acf_entropyflux,
		   const shared_ptr<CF> & acf_numentropyflux,
		   const bool compile)
    : BASE (agfu, atps, "symbolic"),
      cf_flux{acf_flux}, cf_numflux{acf_numflux}, cf_invmap{acf_invmap},
      cf_entropy{acf_entropy}, cf_entropyflux{acf_entropyflux},
      cf_numentropyflux{acf_numentropyflux}
  {
    // set proxies
    proxy_u = aproxy_u;
    proxy_uother = aproxy_uother;

    if(cf_entropy)
      {
	bool wait = false;
	// precompute derivatives for entropy residual
	ddu_invmap = cf_invmap->Diff(proxy_u.get(), proxy_uother);
	ddu_invmap = Compile(ddu_invmap, compile, 0, wait);

	ddphi_invmap = cf_invmap->Diff(BASE::tps->cfgradphi.get(), proxy_graddelta);
	ddphi_invmap = Compile(ddphi_invmap, compile, 0, wait);

	auto temp = cf_entropy - cf_entropyflux*tps->cfgradphi;
	ddu_entropy = temp->Diff(proxy_u.get(), proxy_uother);
	ddu_entropy = Compile(ddu_entropy, compile, 0, wait);
      }
  }

  using BASE::Flux;
  using BASE::NumFlux;
  using BASE::InverseMap;


  /// Given values of grad(φ) (at some fixed pseudotime τ) in "gradphi", 
  /// return the values of "u" at points of a mapped integration rule "mir"
  /// using the saved inverse map y –to–> u where y = g(u) - grad(φ) f(u).  
  
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<SIMD<double>> gradphi,
		  FlatMatrix<SIMD<double>> u) const   {

    // Load the values of u & grad(φ) into the ProxyUserData object of "mir"
    ProxyUserData & ud =
      *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = u;                
    ud.GetAMemory(BASE::tps->cfgradphi.get()) = gradphi;

    // Compute y –to–> u using the inverse map data member (evaluation of its
    // expression tree will use the data we loaded into  "ud" above.
    cf_invmap->Evaluate(mir, u);
  }

  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<SIMD<double>> gradphi,
		  FlatMatrix<SIMD<double>> graddelta,
		  FlatMatrix<SIMD<double>> u,
		  FlatMatrix<SIMD<double>> ut) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = u;         // set values for u
    ud.GetAMemory(proxy_uother.get()) = ut;   // abuse other proxy for derivatives
    ud.GetAMemory(BASE::tps->cfgradphi.get()) = gradphi; // set values for grad(phi)
    ud.GetAMemory(proxy_graddelta.get()) = graddelta;    // set values for graddelta

    STACK_ARRAY(SIMD<double>, mem, COMP*mir.Size());
    FlatMatrix<SIMD<double>> temp(COMP, mir.Size(), mem);

    // map derivative
    ddu_invmap->Evaluate(mir, ut);
    ddphi_invmap->Evaluate(mir, temp);
    ut += temp;
    // map function value
    cf_invmap->Evaluate(mir, u);
  }

  // flux f(u)
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = u; // set values for u
    cf_flux->Evaluate(mir, flux);
  }

  // numerical flux
  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = ul; // set values for ul
    ud.GetAMemory(proxy_uother.get()) = ur; // set values for ur
    cf_numflux->Evaluate(mir,fna);
  }

  // calc \d\hat{t}(\hat{E}) for \hat{E} = E(u) - F(u)*grad(phi) and F(u)
  void CalcEntropy (const SIMD_BaseMappedIntegrationRule & mir,
		    FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> ut,
		    FlatMatrix<SIMD<double>> gradphi, FlatMatrix<SIMD<double>> graddelta,
		    FlatMatrix<SIMD<double>> dEdt, FlatMatrix<SIMD<double>> F) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = u;                  // set values for u
    ud.GetAMemory(proxy_uother.get()) = ut;            // abuse other proxy for derivatives
    ud.GetAMemory(BASE::tps->cfgradphi.get()) = gradphi;   // set values for grad(phi)
    ud.GetAMemory(proxy_graddelta.get()) = graddelta;      // set values for graddelta

    ddu_entropy->Evaluate(mir, dEdt);
    cf_entropyflux->Evaluate(mir, F);
    // add linear part to derivative
    for( size_t i : Range(dEdt.Width()))
      dEdt(0,i) -= InnerProduct(F.Col(i), graddelta.Col(i));
  }

  // numerical entropy flux
  void NumEntropyFlux(const SIMD_BaseMappedIntegrationRule & mir,
		      FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
		      FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = ul;      // set values for ul
    ud.GetAMemory(proxy_uother.get()) = ur; // set values for ur
    cf_numentropyflux->Evaluate(mir,fna);
  }

  void SetViscosityCoefficient(shared_ptr<CoefficientFunction> cf_visc)
  {
    cf_visccoeff = cf_visc;
  }

  void SetNumEntropyFlux(shared_ptr<CoefficientFunction> cf_numentropyflux)
  {
    BASE::cf_numentropyflux = cf_numentropyflux;
  }

  // compute the viscosity coefficient
  void CalcViscCoeffEl(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<SIMD<double>> u,
                       FlatMatrix<SIMD<double>> res,
                       const double hi, double & coeff) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(proxy_u.get()) = u; // set values for u
    ud.GetAMemory(proxy_res.get()) = res;
    cf_visccoeff->Evaluate(mir,res);
    coeff = 0.0;
    for(size_t i : Range(res.Width()))
      for(size_t j : Range(SIMD<double>::Size()))
	if(res(0,i)[j] > coeff)
	  coeff = res(0,i)[j];
  }
};

/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateSymbolicConsLaw (const shared_ptr<GridFunction> & gfu,
						   const shared_ptr<TentPitchedSlab> & tps,
						   const shared_ptr<ProxyFunction> & proxy_u,
						   const shared_ptr<ProxyFunction> & proxy_uother,
						   const shared_ptr<CF> & flux,
						   const shared_ptr<CF> & numflux,
						   const shared_ptr<CF> & invmap,
						   const shared_ptr<CF> & entropy,
						   const shared_ptr<CF> & entropyflux,
						   const shared_ptr<CF> & numentropyflux,
						   const bool compile)
{
  const int dim = tps->ma->GetDimension();
  constexpr int MAXCOMP = 6;
  const int comp_space = gfu->GetFESpace()->GetDimension();
  const auto ecomp = (entropy && entropyflux && numentropyflux) ? 1 : 0;

  shared_ptr<ConservationLaw> cl = nullptr;
  Switch<3>(dim-1, [&](auto DIM) {
      Switch<7>(comp_space, [&](auto COMP) {
	  Switch<2>(ecomp, [&](auto ECOMP) {
	      cl = make_shared<SymbolicConsLaw<DIM.value+1, COMP.value, ECOMP>>(gfu, tps, proxy_u, proxy_uother,
									      flux, numflux, invmap,
									      entropy, entropyflux, numentropyflux,
									      compile);
	    });
	});
    });
  if(cl)
    return cl;
  else
    throw Exception ("Illegal dimension for SymbolicConsLaw");
}
