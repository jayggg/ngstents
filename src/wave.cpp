#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

template <int D>
class Wave : public T_ConservationLaw<Wave<D>, D, D+1, 0>
{
  //whether to use te constitutive parameters mu/epsilon
  //as, e.g., permeability/permittivity in electromagnetics
  bool use_mu_eps = false;
  shared_ptr<CoefficientFunction> cf_mu = nullptr;
  shared_ptr<CoefficientFunction> cf_eps = nullptr;
  typedef T_ConservationLaw<Wave<D>, D, D+1, 0> BASE;
  
public:
  Wave (const shared_ptr<GridFunction> & agfu,
	const shared_ptr<TentPitchedSlab> & atps)
    : BASE (agfu, atps, "wave")
  { };

  using BASE::Flux;
  using BASE::NumFlux;
  using BASE::u_reflect;
  using BASE::u_transparent;
  
  void SetMaterialParameters(shared_ptr<CoefficientFunction> mu,
			     shared_ptr<CoefficientFunction> eps)
  {
    use_mu_eps = true;
    cf_mu = mu;
    cf_eps = eps;
  }

  // solve for û: Û = ĝ(x̂, t̂, û) - ∇̂ φ(x̂, t̂) ⋅ f̂(x̂, t̂, û)
  // at all points in an integration rule
  template <typename T>
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<T> grad, FlatMatrix<T> u) const
  {
    Matrix<T> mu(1, mir.Size());
    Matrix<T> eps(1, mir.Size());
    if(use_mu_eps)
      {
        cf_mu->Evaluate(mir,mu);
        cf_eps->Evaluate(mir,eps);
      }
    for (int i : Range(grad.Width()))
      {
        T mueps = T(1.0);
        if(use_mu_eps)
          {
            mueps = mu(i)*eps(i);
            for(int j : Range(D))
              u(j,i) *= 1.0/mu(i);
          }
	T prod = T(0.0);
	T norm = T(0.0);
	for(int j : Range(D))
	  {
	    prod += grad(j,i) * u(j,i);
	    norm += grad(j,i) * grad(j,i);
	  }

	auto fac = 1.0/(mueps-norm);
	auto p = fac * (u(D,i) + prod);
        for(int j : Range(D))
          u(j,i) += p * grad(j,i);
	u(D,i) = (use_mu_eps) ? mu(i)*p : p;
      }
  }

  template <typename T>
  Mat<D+1,D,typename T::TELEM> Flux (const T & u) const
  {

    Mat<D+1,D,typename T::TELEM> flux;
    flux.Rows(0,D) = u(D)*Id<D>();//0.0*Id<D>();
    flux.Row(D) = u.Range(0,D);

    return flux;
  }
  
  // Flux on element
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    for(int i : Range(mir))
      {
	Mat<D+1,D,SIMD<double>> fluxmat = Flux(u.Col(i));
	flux.Col(i) = fluxmat.AsVector();
      }
  }

  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    for (size_t i : Range(ul.Width()))
      {
        auto nvi = normals.Col(i);
        Vec<D+1,SIMD<double>> flux = 0.5 * (Flux(ul.Col(i))+Flux(ur.Col(i))) * nvi;
	SIMD<double> norm = 0.0;
	SIMD<double> prod = 0.0;
	for (size_t j : Range(D))
	  {
	    norm += nvi(j) * nvi(j);
	    prod += (ul(j,i) - ur(j,i)) * nvi(j);
	  }
	prod /= (2.0*sqrt(norm));
	for (size_t j : Range(D))
          fna(j,i) = flux(j) + prod * nvi(j);
	fna(D,i) = flux(D) + 0.5*sqrt(norm)*(ul(D,i) - ur(D,i));
      }
  }
  
  void u_reflect(const SIMD_BaseMappedIntegrationRule & mir,
		 FlatMatrix<SIMD<double>> u,
                 FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const
  {
    /////////////////// standing wave ///////////////////////
    // u_refl.Rows(0,D) = u.Rows(0,D);
    // u_refl.Row(D) = -u.Row(D);
    /////////////////// reflect wave ///////////////////////
    u_refl.Rows(0,D) = u.Rows(0,D);
    for (int i : Range(u.Width()))
      {
        SIMD<double> prod = 0.0;
        SIMD<double> norm2 = 0.0;
        for(int j : Range(D))
          {
            norm2 += normals(j,i) * normals(j,i); 
            prod += normals(j,i) * u(j,i);
          }
	u_refl.Rows(0,D).Col(i) -= 2/norm2 * prod * normals.Col(i);
      }
    u_refl.Row(D) = u.Row(D);
  }

  void u_transparent(const SIMD_BaseMappedIntegrationRule & mir,
                     FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                     FlatMatrix<SIMD<double>> u_transp) const
  {
    Matrix<SIMD<double>> mu(1, mir.Size());
    Matrix<SIMD<double>> eps(1, mir.Size());
    if(use_mu_eps)
      {
        cf_mu->Evaluate(mir,mu);
        cf_eps->Evaluate(mir,eps);
      }
    u_transp.Rows(0,D) = u.Rows(0,D);
    for (int i : Range(u.Width()))
      {
        SIMD<double> prod = 0.0;
        for(int j : Range(D))
          prod += normals(j,i) * u(j,i);
        u_transp(D,i) = (use_mu_eps) ? sqrt(mu(i)/eps(i)) * prod : prod;
      }
  }
};

/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateWave(const shared_ptr<GridFunction> & gfu,
				       const shared_ptr<TentPitchedSlab> & tps)
{
  const int dim = tps->ma->GetDimension();
  switch(dim){
  case 1:
    return make_shared<Wave<1>>(gfu, tps);
  case 2:
    return make_shared<Wave<2>>(gfu, tps);
  case 3:
    return make_shared<Wave<3>>(gfu, tps);
  }
  throw Exception ("Illegal dimension for Wave");
}
