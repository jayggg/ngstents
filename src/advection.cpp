#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

template <int D>
class Advection : public T_ConservationLaw<Advection<D>, D, 1, 0>
{
  shared_ptr<CoefficientFunction> bfield = nullptr;

  typedef T_ConservationLaw<Advection<D>, D, 1, 0> BASE;
  
public:
  Advection (const shared_ptr<GridFunction> & agfu,
	     const shared_ptr<TentPitchedSlab> & atps)
    : BASE (agfu, atps, "advection")
  { };

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
    STACK_ARRAY(SIMD<double>, mem, D*mir.Size());
    FlatMatrix<SIMD<double>> bmat(D, mir.Size(), mem);

    bfield->Evaluate (mir, bmat);
    for (size_t i : Range(mir))
      {
        SIMD<double> ip = 1.0;
        for (size_t j : Range(D))
          ip -= grad(j,i) * bmat(j,i);
        u(0,i) *= (1/ip);
      }
  }

  // Flux f(u) = b*u on element 
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    bfield->Evaluate(mir, flux);
    for (size_t i : Range(mir))
      flux.Col(i) *= u(0,i);
  }

  // numerical flux
  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    STACK_ARRAY(SIMD<double>, mem, D*mir.Size());
    FlatMatrix<SIMD<double>> bmat(D, mir.Size(), mem);
    bfield->Evaluate(mir, bmat);

    for(size_t i : Range(mir))
      {
        SIMD<double> bn = 0.0;
        for(size_t j : Range(D))
          bn += bmat(j,i)*normals(j,i);

        fna(0,i) = IfPos(bn, bn*ul(0,i), bn*ur(0,i));
      }
  }
};

/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateAdvection (const shared_ptr<GridFunction> & gfu,
					     const shared_ptr<TentPitchedSlab> & tps)
{
  const int dim = tps->ma->GetDimension();
  switch(dim){
  case 1:
    return make_shared<Advection<1>>(gfu, tps);
  case 2:
    return make_shared<Advection<2>>(gfu, tps);
  case 3:
    return make_shared<Advection<3>>(gfu, tps);
  }
  throw Exception ("Illegal dimension for Advection");
}
