#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

template <int D>
class SymbolicConsLaw : public T_ConservationLaw<SymbolicConsLaw<D>, D, 1, 0, true>
{
  shared_ptr<CoefficientFunction> bfield = nullptr;

  typedef T_ConservationLaw<SymbolicConsLaw<D>, D, 1, 0, true> BASE;

  shared_ptr<CoefficientFunction> cf_flux = nullptr;
  shared_ptr<CoefficientFunction> cf_numflux = nullptr;
  shared_ptr<CoefficientFunction> cf_invmap = nullptr;
public:
  Array<ProxyFunction*> flux_proxies;
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
		    if (!flux_proxies.Contains(proxy))
		      flux_proxies.Append (proxy);
		}
	    });
	cout << "flux_proxies: " << flux_proxies << endl;
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
  //
  // in this case, û = Û / [ 1 - b̂ ⋅∇̂ φ(x̂, t̂) ]
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<SIMD<double>> grad, FlatMatrix<SIMD<double>> u) const
  {
    // STACK_ARRAY(SIMD<double>, mem, D*mir.Size());
    // FlatMatrix<SIMD<double>> bmat(D, mir.Size(), mem);
    // 
    // bfield->Evaluate (mir, bmat);
    // for (size_t i : Range(mir))
    //   {
    //     SIMD<double> ip = 1.0;
    //     for (size_t j : Range(D))
    //       ip -= grad(j,i) * bmat(j,i);
    //     u(0,i) *= (1/ip);
    //   }
    // cout << "grad_phi = " << endl << grad << endl;
    cf_invmap->Evaluate(mir, u);
  }

  // Flux f(u) = b*u on element 
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    // bfield->Evaluate(mir, flux);
    // for (size_t i : Range(mir))
    //   flux.Col(i) *= u(0,i);
    cf_flux->Evaluate(mir, flux);
  }

  // numerical flux
  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    // STACK_ARRAY(SIMD<double>, mem, D*mir.Size());
    // FlatMatrix<SIMD<double>> bmat(D, mir.Size(), mem);
    // bfield->Evaluate(mir, bmat);
    // 
    // for(size_t i : Range(mir))
    //   {
    //     SIMD<double> bn = 0.0;
    //     for(size_t j : Range(D))
    //       bn += bmat(j,i)*normals(j,i);
    // 
    //     fna(0,i) = IfPos(bn, bn*ul(0,i), bn*ur(0,i));
    //   }
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
    return make_shared<SymbolicConsLaw<1>>(gfu, tps, flux, numflux, invmap);
  case 2:
    return make_shared<SymbolicConsLaw<2>>(gfu, tps, flux, numflux, invmap);
  case 3:
    return make_shared<SymbolicConsLaw<3>>(gfu, tps, flux, numflux, invmap);
  }
  throw Exception ("Illegal dimension for SymbolicConsLaw");
}
