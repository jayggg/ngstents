#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

template <int D, int COMP>
class SymbolicConsLaw : public T_ConservationLaw<SymbolicConsLaw<D,COMP>, D, COMP, 0, true>
{
  shared_ptr<CoefficientFunction> bfield = nullptr;

  typedef T_ConservationLaw<SymbolicConsLaw<D, COMP>, D, COMP, 0, true> BASE;

  shared_ptr<CoefficientFunction> cf_flux = nullptr;
  shared_ptr<CoefficientFunction> cf_numflux = nullptr;
  shared_ptr<CoefficientFunction> cf_invmap = nullptr;
public:
  ProxyFunction * flux_proxy;
  Array<ProxyFunction*> numflux_proxies;
  ProxyFunction * invmap_proxy;
public:
  SymbolicConsLaw (const shared_ptr<GridFunction> & agfu,
		   const shared_ptr<TentPitchedSlab> & atps,
		   const shared_ptr<CoefficientFunction> & acf_flux,
		   const shared_ptr<CoefficientFunction> & acf_numflux,
		   const shared_ptr<CoefficientFunction> & acf_invmap)
    : BASE (agfu, atps, "symbolic"),
      cf_flux{acf_flux}, cf_numflux{acf_numflux}, cf_invmap{acf_invmap}
  {
    if(cf_flux)
      {
	cf_flux->TraverseTree
	  ( [&] (CoefficientFunction & nodecf)
	    {
	      auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
	      if (proxy)
		{
		  if (!proxy->IsTestFunction())
		    flux_proxy = proxy;
		}
	    });
	cout << "flux_proxy: " << flux_proxy << endl;
      }
    if(cf_numflux)
      {
	cf_numflux->TraverseTree
	  ( [&] (CoefficientFunction & nodecf)
	    {
	      auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
	      if (proxy)
		{
		  if (!proxy->IsTestFunction())
		    if (!numflux_proxies.Contains(proxy))
		      numflux_proxies.Append (proxy);
		}
	    });
	cout << "numflux_proxies: " << numflux_proxies << endl;
      }
    if(cf_invmap)
      {
	cf_invmap->TraverseTree
	  ( [&] (CoefficientFunction & nodecf)
	    {
	      auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
	      if (proxy)
		{
		  if (!proxy->IsTestFunction())
		    invmap_proxy = proxy;
		}
	    });
	cout << "invmap_proxy: " << invmap_proxy << endl;
      }
  };

  using BASE::Flux;
  using BASE::NumFlux;
  using BASE::InverseMap;

  void SetVectorField(shared_ptr<CoefficientFunction> cf) { bfield = cf; }
  
  // solve for û: Û = ĝ(x̂, t̂, û) - ∇̂ φ(x̂, t̂) ⋅ f̂(x̂, t̂, û)
  // at all points in an integration rule
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<SIMD<double>> grad, FlatMatrix<SIMD<double>> u) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(invmap_proxy) = u;                  // set values for u
    ud.GetAMemory(BASE::tps->cfgradphi.get()) = grad; // set values for grad(phi)
    cf_invmap->Evaluate(mir, u);
  }

  // flux f(u)
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(flux_proxy) = u; // set values for u
    cf_flux->Evaluate(mir, flux);
  }

  // numerical flux
  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    ProxyUserData & ud = *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    ud.GetAMemory(numflux_proxies[0]) = ul; // set values for ul
    ud.GetAMemory(numflux_proxies[1]) = ur; // set values for ur
    cf_numflux->Evaluate(mir,fna);
  }
};

/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateSymbolicConsLaw (const shared_ptr<GridFunction> & gfu,
						   const shared_ptr<TentPitchedSlab> & tps,
						   const shared_ptr<CoefficientFunction> & flux,
						   const shared_ptr<CoefficientFunction> & numflux,
						   const shared_ptr<CoefficientFunction> & invmap)
{
  const int dim = tps->ma->GetDimension();
  switch(dim){
  case 1:
    Switch<MAXCOMP>(comp_space, [&](auto COMP) {
	cl = make_shared<SymbolicConsLaw<1, COMP>>(gfu, tps, flux, numflux, invmap);
      });
    break;
  case 2:
    Switch<MAXCOMP>(comp_space, [&](auto COMP) {
	cl = make_shared<SymbolicConsLaw<2, COMP>>(gfu, tps, flux, numflux, invmap);
      });
    break;
  case 3:
    Switch<MAXCOMP>(comp_space, [&](auto COMP) {
	cl = make_shared<SymbolicConsLaw<3, COMP>>(gfu, tps, flux, numflux, invmap);
      });
    break;
  }
  if(cl)
    return cl;
  else
    throw Exception ("Illegal dimension for SymbolicConsLaw");
}
