#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

template <int D>
class Maxwell : public T_ConservationLaw<Maxwell<D>, D, 2*D, 0> 
{
  // shared_ptr<GridFunction> gfE;
  // shared_ptr<GridFunction> gfH;
  typedef T_ConservationLaw<Maxwell<D>, D, 2*D, 0> BASE;
  
public:
  Maxwell (const shared_ptr<TentPitchedSlab> & atps, const int & order)
    : BASE (atps, "maxwell", order)
  { 
    // shared_ptr<FESpace> fesD =
    //   CreateFESpace("l2ho", ma, Flags().SetFlag("order",order).SetFlag("dim",D).SetFlag("all_dofs_together"));
    // fesD->Update();
    // fesD->FinalizeUpdate();
    // gfE = CreateGridFunction(fesD,"E",Flags());
    // gfE->Update();
    // gfH = CreateGridFunction(fesD,"H",Flags());
    // gfH->Update();
  }
  
  using BASE::Flux;
  using BASE::NumFlux;
  using BASE::InverseMap;
  
  template<typename T>
  INLINE auto skew (T vec) const
  {
    if(D!=3) throw Exception("not implemented for D!=3");
    Mat<D,D,typename T::TSCAL> mat;
    mat(0,0) = 0.0;
    mat(0,1) = -vec(2);
    mat(0,2) = vec(1);
    mat(1,0) = vec(2);
    mat(1,1) = 0.0;
    mat(1,2) = -vec(0);
    mat(2,0) = -vec(1);
    mat(2,1) = vec(0);
    mat(2,2) = 0.0;
    return mat;
  }

  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    for(size_t i : Range(u.Width()))
      {
	Mat<2*D,D,SIMD<double>> fluxmat;
        fluxmat.Rows(0,D) = skew(u.Col(i).Range(D,2*D));
	fluxmat.Rows(D,2*D) = -skew(u.Col(i).Range(0,D));

        size_t l = 0;
	for(size_t j : Range(D))
	  for(size_t k : Range(2*D))
	    flux(l++,i) = fluxmat(k,j);
      }
  }

  Vec<2*D> Flux (const FlatVec<2*D> & ul, const FlatVec<2*D> & ur, const Vec<D> & nv) const
  {
    Vec<2*D> flux = 0.5*(Flux(ul)+Flux(ur))*nv;
    // double N2 = L2Norm2(nv);
    // flux.Range(0,D) += 1/N2 * InnerProduct(ul.Range(0,D)-ur.Range(0,D), nv) * nv;
    // flux.Range(D,2*D) += 1/N2 * InnerProduct(ul.Range(0,D)-ur.Range(0,D), nv) * nv;

    //// alternating flux
    // Vec<2*D> flux;
    // flux.Rows(0,D) = Flux(ur).Rows(0,D)*nv;
    // flux.Rows(D,2*D) = Flux(ul).Rows(D,2*D)*nv;
    return flux;
  }

  auto NumFlux (const Vec<6,SIMD<double>> & ul,
                const Vec<6,SIMD<double>> & ur,
                const Vec<3,SIMD<double>> nvi) const -> Vec<6,SIMD<double>>
  {
    // central flux
    Vec<3,SIMD<double>> el(ul(0), ul(1), ul(2));
    Vec<3,SIMD<double>> hl(ul(3), ul(4), ul(5));
    Vec<3,SIMD<double>> er(ur(0), ur(1), ur(2));
    Vec<3,SIMD<double>> hr(ur(3), ur(4), ur(5));
    Vec<3,SIMD<double>> emean = 0.5 * (el+er);
    Vec<3,SIMD<double>> hmean = 0.5 * (hl+hr);
    auto exn = Cross (emean, nvi);
    auto hxn = Cross (hmean, nvi);

    // jump terms
    Vec<3,SIMD<double>> jumpe = el-er;
    Vec<3,SIMD<double>> jumph = hl-hr;
    SIMD<double> invN2 = SIMD<double>(1.0) / sqrt(L2Norm2(nvi));
    Vec<3,SIMD<double>> jumpe_tau = 0.5 * invN2 * Cross(Cross (jumpe, nvi), nvi);
    Vec<3,SIMD<double>> jumph_tau = 0.5 * invN2 * Cross(Cross (jumph, nvi), nvi);
    
    Vec<6,SIMD<double>> flux;
    for (size_t j = 0; j < 3; j++)
      {
        flux(j) = hxn(j) - jumpe_tau(j);
        flux(j+3) = -exn(j) - jumph_tau(j);
      }
    return flux;
  }

  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    for (size_t i = 0; i < ul.Width(); i++)
      fna.Col(i) = NumFlux(ul.Col(i), ur.Col(i), normals.Col(i));
  }
  
  void TransformBack(const BaseMappedIntegrationPoint & mip,
		     const Vec<D> & gradphi,
		     const FlatVec<2*D> u) const
  {
    Mat<2*D> mat;
    mat.Rows(0,D).Cols(0,D) = Id<D>();
    mat.Rows(0,D).Cols(D,2*D) = skew(gradphi);
    mat.Rows(D,2*D).Cols(0,D) = -skew(gradphi);
    mat.Rows(D,2*D).Cols(D,2*D) = Id<D>();
    CalcInverse(mat);
    u = mat * u;
  }

  
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<SIMD<double>> grad, FlatMatrix<SIMD<double>> u) const
  {
    /* MapBack: solves uhat = u - f(u)*grad for u
       For this linear flux function f, we can write the above equation as uhat = M * u,
       where M = [I,skew(grad);-skew(grad),I] and I the 3x3 identity matrix.
       MapBack applies the M^-1 on the fly.
     */
    auto MapBack = [&](int i, auto uhat)
      {
        auto a = grad(0,i);
        auto b = grad(1,i);
        auto c = grad(2,i);
	auto tmp1 = a * a + b * b + c * c;
	for(int j : Range(tmp1.Size()))
	  if((1.0 - tmp1[j]) < 1e-8)
	    cout << "norm gradphi = " << tmp1[j] << ", choose higher wave speed" << endl;
	
        auto tmp = 1.0/(a * a + b * b + c * c - 1);

	Vec<2*D,decltype(a)> res;
        Vec<2*D,decltype(a)> uin = tmp * uhat;
        
        res(0) = (a*a-1) * uin(0) + b*a * uin(1) + c*a * uin(2) - c * uin(4) + b * uin(5);
        res(1) = b*a * uin(0) + (b*b-1) * uin(1) + c*b * uin(2) + c * uin(3) - a * uin(5);
        res(2) = c*a * uin(0) + c*b * uin(1) + (c*c-1) * uin(2) - b * uin(3) + a * uin(4);
        
        res(3) = c * uin(1) - b * uin(2) + (a*a-1) * uin(3) + b*a * uin(4) + c*a * uin(5);
        res(4) =-c * uin(0) + a * uin(2) + b*a * uin(3) + (b*b-1) * uin(4) + c*b * uin(5);
        res(5) = b * uin(0) - a * uin(1) + c*a * uin(3) + c*b * uin(4) + (c*c-1) * uin(5);
	return res;
      };

    for (int i : Range(mir))
      u.Col(i) = MapBack (i, u.Col(i));
  }

  void u_reflect(const SIMD_BaseMappedIntegrationRule & mir,
		 FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const
  {
    // dirichlet bcs
    u_refl.Rows(0,D) = -u.Rows(0,D);
    u_refl.Rows(D,2*D) = u.Rows(D,2*D);
  }

//   virtual void UpdateVisualisation(const BaseVector & hu, LocalHeap & lh) const
//   {
//     HeapReset hr(lh);
//     int n = gfE->GetVector().FVDouble().Size()/D;
//     FlatMatrixFixWidth<2*D> u(n, &hu.FV<double>()(0));
//     FlatMatrixFixWidth<D> E(n,&gfE->GetVector().FVDouble()(0));
//     FlatMatrixFixWidth<D> H(n,&gfH->GetVector().FVDouble()(0));
//     E = u.Cols(0,D);
//     H = u.Cols(D,2*D);
//   }
}; 

/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateMaxwell(const shared_ptr<TentPitchedSlab> & tps, const int & order)
{
  const int dim = tps->ma->GetDimension();
  if(dim == 3)
    return make_shared<Maxwell<3>>(tps,order);
  else
    throw Exception("Maxwell equations not implemented for D != 3");
}
