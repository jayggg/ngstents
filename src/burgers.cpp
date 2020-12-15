#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"
#include <python_ngstd.hpp>

template <int D>
class Burgers : public T_ConservationLaw<Burgers<D>,D,1,1,false>
{
  typedef T_ConservationLaw<Burgers<D>,D,1,1,false> BASE;

public:
  Burgers (const shared_ptr<TentPitchedSlab> & tps, const int & order)
    : BASE (tps, "burgers", order) { ; }

  // these two were private
  using BASE::gfnu;
  using BASE::gfres;

  using BASE::CalcViscCoeffEl;
  using BASE::Flux;
  using BASE::u_reflect;
  using BASE::TransformBackIR;
  using BASE::CalcEntropy;


  // solve for û: Û = ĝ(x̂, t̂, û) - ∇̂ φ(x̂, t̂) ⋅ f̂(x̂, t̂, û)
  //
  // in this case, û = 2 Û / [ 1 + √(1 - 2 Û ⋅∇̂ φ(x̂, t̂)) ]
  //
  template <typename T>
  void TransformBackIR(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<T> grad, FlatMatrix<T> u) const
  {
    for (int i : Range(grad.Width()))
      {
	auto sum = T(0.0);
	for(size_t j : Range(D))
	  sum += grad(j,i);
	u(0,i) = 2 * u(0,i) / (1.0 + sqrt(1.0 - 2.0*sum*u(0,i)) );
      }
  }

  // Flux on element
  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    for(size_t i : Range(mir))
      flux.Col(i) = 0.5 * u(0,i)*u(0,i);
  }

  // Numerical Flux on a facet.  ul and ur are the values at integration points
  // of the two elements adjacent to an internal facet
  // of the spatial mesh of a tent.
  void Flux(FlatMatrix<SIMD<double>> ula, FlatMatrix<SIMD<double>> ura,
            FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    for(size_t i : Range(ula.Width()))
      {
        SIMD<double> sumn = 0.0;
        for(size_t j : Range(D))
          sumn += normals(j,i);

        SIMD<double> um = 0.5*(ula(0,i) + ura(0,i));
        SIMD<double> fprimen = um * sumn; // f'(u) * n gives characteristic speed
        fna(0,i) = IfPos(fprimen,
                         0.5*(ula(0,i)*ula(0,i))*sumn,
                         0.5*(ura(0,i)*ura(0,i))*sumn);
      }
  }

  void u_reflect(FlatMatrix<SIMD<double>> u,
                 FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const
  {
    cout << "no reflecting B.C. for Burgers Equation" << endl;
  }

  void CalcEntropy(FlatMatrix<AutoDiff<1,SIMD<double>>> adu,
                   FlatMatrix<AutoDiff<1,SIMD<double>>> grad,
		   FlatMatrix<SIMD<double>> dEdt, FlatMatrix<SIMD<double>> F) const
  {
    // E = u^2/2,
    // F = u^3/3 [ 1 ]*D
    for(size_t i : Range(adu.Width()))
      {
        auto ui = adu(0,i);
        auto sumgrad = grad(0,i);
        auto Fi = ui*ui*ui/3.0;
        for(size_t j : Range(1,D))
          sumgrad += grad(j,i);

        AutoDiff<1,SIMD<double>> adE = ui*ui/2.0 - Fi*sumgrad; // E - (F,grad(phi))
        dEdt(0,i) = adE.DValue(0);
        for(size_t j : Range(D))
          F(j,i) = Fi.Value();
      }
  }


  void EntropyFlux (FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
                    FlatMatrix<SIMD<double>> n,
                    FlatMatrix<SIMD<double>> flux) const
  {
    for(size_t i : Range(ul.Width()))
      {
        SIMD<double> sumn = 0.0;
        for(size_t j : Range(D))
          sumn += n(j,i);
	//////// upwind flux ////////
	auto um = 0.5 * (ul(0,i) + ur(0,i));
	auto Fprimen = um*um * sumn;
	flux(0,i) = IfPos(Fprimen,
			  1./3*ul(0,i)*ul(0,i)*ul(0,i)*sumn,
			  1./3*ur(0,i)*ur(0,i)*ur(0,i)*sumn);
      }
  }

  // Given the entropy, entropy residual and 'hi' in the cylinder
  // at points in a MIR, where
  // hi = ([Jacobian determinant for element]/DIM)^{1/DIM}
  //
  // compute the viscosity coefficient
  // ν = max_{T ∈  Tᵢ} min(ν_*ᵀ, νₑᵀ), where
  // ν_*ᵀ = κ₂ diam(T)||Dᵤ f̂||_{L^∞(T)},
  // vₑᵀ = c_X² ||Rₕ||_{L^∞ (T)} / |E̅|, where E̅ is the mean value of entropy
  // c_X² is an effective local grid size of X, Rₕ is entropy residual
  //
  // betai corresponds to v_*ᵀ and visci * hi/Emean corresponds to vₑᵀ.
  void CalcViscCoeffEl(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<SIMD<double>> elu_ipts,
                       FlatMatrix<SIMD<double>> res_ipts,
                       const double hi, double & coeff) const
  {
    int nipt = mir.IR().GetNIP();

    double betai = 0.0;
    double visci = 0.0;
    double Emean = 0.0;

    for(size_t i : Range(elu_ipts.Width()))
      {
        auto res = res_ipts(0,i);
        auto ui = elu_ipts(0,i);
        for(size_t j : Range(SIMD<double>::Size()))
          {
            if(fabs(res[j]) > visci)
              visci = fabs(res[j]);
            Emean += ui[j]*ui[j]/2.0;
            betai = max(betai, fabs(ui[j]));
          }
      }
    Emean /= nipt;
    coeff = min (visci * hi/Emean, betai);
    coeff *= hi;
  }
};

/////////////////////////////////////////////////////////////////////////
//                 EXPORT TO PYTHON
/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateBurgers(const shared_ptr<TentPitchedSlab> & tps, const int & order)
{
  int dim = tps->ma->GetDimension();
  switch(dim){
  case 1:
    return make_shared<Burgers<1>>(tps,order);
  case 2:
    return make_shared<Burgers<2>>(tps,order);
  }
  throw Exception ("Burgers only available for 1D and 2D");
}
