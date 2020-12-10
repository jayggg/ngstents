#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"
#include <python_ngstd.hpp>

template <int D>
class Wave : public T_ConservationLaw<Wave<D>,D,D+1,0,false>
{
  shared_ptr<GridFunction> gfU;
  //whether to use te constitutive parameters mu/epsilon
  //as, e.g., permeability/permittivity in electromagnetics
  bool use_mu_eps = false;
  shared_ptr<CoefficientFunction> cf_mu = nullptr;
  shared_ptr<CoefficientFunction> cf_eps = nullptr;
  typedef T_ConservationLaw<Wave<D>,D,D+1,0,false> BASE;
  
public:
  Wave (const shared_ptr<TentPitchedSlab> &atps, const int & order)
    : BASE (atps, "wave", order)
  {
    shared_ptr<FESpace> fesvel =
      CreateFESpace("l2ho",this->tps->ma, Flags().SetFlag("order",order).SetFlag("dim",D).SetFlag("all_dofs_together"));
    fesvel->Update();
    fesvel->FinalizeUpdate();
    gfU = CreateGridFunction(fesvel,"U",Flags());
    gfU->Update();
  }

  // these two were private
  using BASE::gfnu;
  using BASE::gfres;

  using BASE::CalcViscCoeffEl;
  using BASE::Flux;
  using BASE::u_reflect;
  using BASE::u_transparent;
  using BASE::TransformBackIR;
  using BASE::CalcEntropy;

  virtual void SetMaterialParameters(shared_ptr<CoefficientFunction> mu,
                                     shared_ptr<CoefficientFunction> eps)
  {
    use_mu_eps = true;
    cf_mu = mu;
    cf_eps = eps;
  }

  // solve for û: Û = ĝ(x̂, t̂, û) - ∇̂ φ(x̂, t̂) ⋅ f̂(x̂, t̂, û)
  //
  // in this case, û = 2 Û / [ 1 + √(1 - 2 Û ⋅∇̂ φ(x̂, t̂)) ]
  //
  template <typename MIP=BaseMappedIntegrationPoint, typename SCAL=double>
  void TransformBack(const MIP & mip,
                     const Vec<D,SCAL > & grad,
                     const FlatVec<1,SCAL > u) const
  {
    const double p = (u(D)+InnerProduct(grad, u.Range(0,D)))/(1.0-L2Norm2(grad));
    u.Range(0,D) += p*grad;
    u(D) = p;
  }

  // solve for û at all points in an integration rule
  template <typename T>
  void TransformBackIR(const SIMD_BaseMappedIntegrationRule & mir,
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

  template <typename T>//, typename SCAL=double>
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
	Mat<D+1,D,SIMD<double>> fluxmat;
        fluxmat = Flux(u.Col(i));
	size_t l = 0;
	for(size_t j : Range(D))
	  for(size_t k : Range(D+1))
	    flux(l++,i) = fluxmat(k,j);
      }
  }

   Vec<D+1> Flux (const Vec<D+1> & ul, const Vec<D+1> & ur, const Vec<D> & nv) const
  {
    /*
    Mat<D+1,D> flux;
    flux.Rows(0,D) = (ul(D)+ur(D))*Id<D>();// 0.0*Id<D>();
    flux.Row(D) = ul.Range(0,D)+ur.Range(0,D);
    flux.Row(D) += 0.05*(ul(D)-ur(D)) * nv;
    return 0.5*flux*nv;
    */
    Vec<D+1> flux;
    flux.Range(0,D) = 0.5 * (ul(D)+ur(D)) * nv;
    flux(D) = 0.5 * InnerProduct (ul.Range(0,D)+ur.Range(0,D), nv);

    double N2 = L2Norm2(nv);
    flux.Range(0,D) += 1/N2 * InnerProduct(ul.Range(0,D)-ur.Range(0,D), nv) * nv;
    flux(D) += 1*(ul(D)-ur(D));
    return flux;
  }

  void Flux(FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
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
  
  void u_reflect(FlatMatrix<SIMD<double>> u,
                 FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const override
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
        u_refl.Col(i) -= 2/norm2 * prod * normals.Col(i);
      }
    u_refl.Row(D) = u.Row(D);
  }

  void u_transparent(SIMD_BaseMappedIntegrationRule & mir,
                     FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                     FlatMatrix<SIMD<double>> u_transp) const override
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

  virtual void UpdateVisualisation(const BaseVector & hu, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    int n = gfU->GetVector().FVDouble().Size()/D;
    FlatMatrixFixWidth<D+1> u(n, &hu.FV<double>()(0));
    FlatMatrixFixWidth<D> U(n,&gfU->GetVector().FVDouble()(0));
    U = u.Cols(0,D);
  }
};

/////////////////////////////////////////////////////////////////////////
//                 EXPORT TO PYTHON
/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateWave(const shared_ptr<TentPitchedSlab> & tps, const int & order)
{
  const int dim = tps->ma->GetDimension();
  switch(dim){
  case 1:
    return make_shared<Wave<1>>(tps,order);
  case 2:
    return make_shared<Wave<2>>(tps,order);
  case 3:
    return make_shared<Wave<3>>(tps,order);
  }
  throw Exception ("Illegal dimension for Wave");
}
