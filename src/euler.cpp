#include <solve.hpp>
using namespace ngsolve;

#include "tconservationlaw_tp_impl.hpp"

double dim_ = 5; // degrees of freedom of gas molecules
double gamma_ = 1.4; // (dim_+2)/dim_

double erf(double x)
{
  /* erf(z) = 2/sqrt(pi) * Integral(0..x) exp( -t^2) dt
     erf(0.01) = 0.0112834772 erf(3.7) = 0.9999998325
     Abramowitz/Stegun: p299, |erf(z)-erf| <= 1.5*10^(-7)
  */
  if (x < 0) return -erf(-x);
  
  double y = 1.0 / ( 1.0 + 0.3275911 * x);
  return 1 - (((((
		  + 1.061405429 * y
		  - 1.453152027) * y
		 + 1.421413741) * y
		- 0.284496736) * y
	       + 0.254829592) * y)
    * exp (-x * x);
}

inline void Int_x_infty (double x, double & int0, double & int1, double & int2, double & int3)
{
  // int_x^\infty  exp(-v^2) v^i dv
  
  int0 = sqrt(M_PI)/2 * (1-erf(x));
  int1 = 0.5 * exp(-x*x);
  int2 = 0.5 * int0 + x * int1;
  int3 = int1 * (1+x*x);
}

template <int D>
class Euler : public T_ConservationLaw<Euler<D>, D, D+2, 1> 
{
  typedef T_ConservationLaw<Euler<D>, D, D+2, 1> BASE;
  
  shared_ptr<GridFunction> gfU;
  shared_ptr<GridFunction> gfrho;
  shared_ptr<GridFunction> gfp;
  shared_ptr<GridFunction> gfT;
  shared_ptr<GridFunction> gfmach;
  using BASE::gfnu;
  using BASE::gfres;
  using BASE::fes;
  using BASE::ma;
  using BASE::pylh;

public: 
  Euler (const shared_ptr<GridFunction> & agfu,
	 const shared_ptr<TentPitchedSlab> & atps)
    : BASE (agfu, atps, "euler")
    {
      int order = fes->GetOrder();
      shared_ptr<FESpace> fesvel =
	CreateFESpace("l2ho", ma, Flags().SetFlag("order",order).SetFlag("dim",D).SetFlag("all_dofs_together"));
      fesvel->Update();
      fesvel->FinalizeUpdate();
      gfU = CreateGridFunction(fesvel,"U",Flags());
      gfU->Update();
      
      shared_ptr<FESpace> fes_scal = gfres->GetFESpace();
      gfrho = CreateGridFunction(fes_scal,"rho",Flags());
      gfrho->Update();
      gfp = CreateGridFunction(fes_scal,"p",Flags());
      gfp->Update();
      gfT = CreateGridFunction(fes_scal,"T",Flags());
      gfT->Update();
      gfmach = CreateGridFunction(fes_scal,"mach",Flags());
      gfmach->Update();
    }
  
  using BASE::u_reflect;
  using BASE::InverseMap;
  using BASE::CalcEntropy;

  template <typename T>//SCAL=double>
  Mat<D+2,D,typename T::TELEM> Flux (const T & U) const
  {
    Mat<D+2,D,typename T::TELEM> flux;
    auto rho = U(0);
    Vec<D,typename T::TELEM> u = 1/rho * U.Range(1,D+1);
    auto E = U(D+1);
    auto e = E/rho - 0.5 * InnerProduct(u,u); // L2Norm2(u) wouldn't work for AutoDiff
    auto temp = 4.0 * e / dim_;

    flux.Row(0) = rho * u;
    flux.Rows(1,D+1) = 0.5 *rho*temp*Id<D>() + rho *u*Trans(u);
    flux.Row(D+1) = ( E + 0.5*rho*temp ) * u;
    return flux;
  }

  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    for(int i : Range(u.Width()))
      {
        Mat<D+2,D,SIMD<double>> fluxmat = Flux(u.Col(i));
	flux.Col(i) = fluxmat.AsVector();
      }
  }

  Vec<D+2> NumFlux (Vec<D+2> ul, Vec<D+2> ur, Vec<D> n) const
    {
      /* // cout << "ul, ur = " << ul << ", " << ur << endl;
      // left value
      double rhol = ul(0);
      Vec<D> Ul = ul.Range(1,D+1); // rho * u
      double Tl = 2*(2*rhol*ul(D+1) - L2Norm2(Ul)) / (dim_*rhol*rhol);
      double betal = L2Norm(Ul) + sqrt(gamma_*Tl);
      // cout << Tl << ", " << betal << endl;
      // right value
      double rhor = ur(0);
      Vec<D> Ur = ur.Range(1,D+1); // rho * u
      double Tr = 2*(2*rhor*ur(D+1) - L2Norm2(Ur)) / (dim_*rhor*rhor);
      double betar = L2Norm(Ur) + sqrt(gamma_*Tr);
      // cout << Tr << ", " << betar << endl;
      // cout << "betamax = " << max(betar,betal) << endl;
      return 0.5 * (Flux(ul)*n + Flux(ur)*n) + L2Norm(n) * 0.5*max(betar,betal) * (ul-ur);
      // return 0.5 * (Flux(ul)*n + Flux(ur)*n) + L2Norm(n) * (ul-ur);
      */
      Vec<D+2> h = 0.0;
      
      double dim = dim_;
      
      for (int side = 0; side < 2; side++)
	{
	  Vec<D+2> ut = (side == 0) ? ul : ur;
	  Vec<D> ntn = (side == 0) ? n : Vec<D>(-n);
	  double fac = (side == 0) ? 1 : -1;
	  
	  double len = L2Norm (n);
	  ntn /= len;
	  
	  double rho = ut(0);
	  Vec<D> U = 1/rho * ut.Range(1,D+1);
	  
	  double un = InnerProduct (U, ntn);
	  double normu2 = L2Norm2 (U);
	  
	  double e = ut(D+1)/rho - 0.5 * normu2;
	  
	  double T = 4.0/dim * e;
	  
	  if (T <= 0)
	    {
	      // cout << "T correction" << endl;
	      T = 1e-10;
	    }
	  
	  if (rho <= 0)
	    {
	      // throw Exception ("rho negative");
	      // cout << "rho negative" << endl;
	    }
	  
	  double sqrtT = sqrt(T);
	  
	  double i0, i1, i2, i3;
	  Int_x_infty (-un/sqrtT, i0, i1, i2, i3);
	  
	  fac *= len * rho / sqrt(M_PI);
	  
	  h(0) += fac * (un*i0 + sqrtT*i1);
	  
	  for (int j = 0; j < D; j++)
	    h(1+j) += fac * (U(j)*un*i0 + sqrtT * (un*ntn(j) + U(j)) * i1 + T * ntn(j) * i2);
	  
	  h(D+1) += 0.5 * fac * (un* (normu2+0.5*(dim-1)*T) * i0 +
				 (2*un*un+normu2 + 0.5*(dim-1)*T) * sqrtT* i1 +
				 3*un * T * i2 + T * sqrtT * i3);
	}
      
      return h;
    }

  void NumFlux(FlatMatrix<SIMD<double>> ula, FlatMatrix<SIMD<double>> ura,
            FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    // TODO: proper use of SIMD!
    for (size_t i : Range(ula.Width()))
      for (size_t j : Range(SIMD<double>::Size()))
        {
          Vec<D+2> ul, ur;
          Vec<D> n;
          for (size_t k : Range(D+2))
            {
              ul(k) = ula(k,i)[j];
              ur(k) = ura(k,i)[j];
            }
          for (size_t k : Range(D))
            n(k) = normals(k,i)[j];
          Vec<D+2> fn = NumFlux(ul,ur,n);
          for (size_t k : Range(D+2))
            reinterpret_cast<double*>(&fna(k,i))[j] = fn(k);
        }
  }

  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
            FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
            FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    NumFlux (ul, ur, normals, fna);
  }

  void u_reflect(Vec<D+2> U, Vec<D> n, Vec<D+2> & U_refl) const
  {
    double rho = U(0);
    Vec<D> u = 1/rho * U.Range(1,D+1);
    double E = U(D+1) / rho;
    double e = E - 0.5 * L2Norm2(u);
    
    Vec<D> nn = n;
    nn /= L2Norm(n);
    
    Vec<D> u_refl = u - 2 * InnerProduct(nn, u) * nn;

    U_refl(0) = rho; 
    U_refl.Range(1,D+1) = rho * u_refl;
    U_refl(D+1) = rho * ( e + 0.5 * L2Norm2(u_refl));
  }

  void u_reflect(const SIMD_BaseMappedIntegrationRule & mir,
		 FlatMatrix<SIMD<double>> u,
		 FlatMatrix<SIMD<double>> normals,
		 FlatMatrix<SIMD<double>> u_refl) const
  {
    for(auto i : Range(u.Width()))
      {
        for (size_t j : Range(SIMD<double>::Size()))
          {
            Vec<D+2> uvec, uvec_refl;
            for(auto k : Range(D+2))
              uvec(k) = u(k,i)[j];
            Vec<D> n;
            for(auto k : Range(D))
              n(k) = normals(k,i)[j];
            u_reflect(uvec,n,uvec_refl);
            for(auto k : Range(D+2))
              reinterpret_cast<double*>(&u_refl(k,i))[j] = uvec_refl(k);
          }
      }
  }

  void CalcEntropy(FlatMatrix<AutoDiff<1,SIMD<double>>> adu,
                   FlatMatrix<AutoDiff<1,SIMD<double>>> grad, 
		   FlatMatrix<SIMD<double>> dEdt, FlatMatrix<SIMD<double>> F) const
  {
    auto InnerProduct = [](auto a, auto b)
      {
        typename decltype(a)::TELEM res(0.0);
        for (size_t i : Range(a))
          res += a(i)*b(i);
        return res;
      };
    
    for(size_t i : Range(adu.Width()))
      {
        auto ucol = adu.Col(i);
        AutoDiff<1,SIMD<double>> adnormu2 = InnerProduct(ucol.Range(1,D+1),ucol.Range(1,D+1));
        AutoDiff<1,SIMD<double>> adT =  2*(2*adu(0,i)*adu(D+1,i) - adnormu2) / (dim_*adu(0,i)*adu(0,i));
        adT.Value() = IfPos((-1)*adT,SIMD<double>(1e-10),adT.Value());
        
        AutoDiff<1,SIMD<double>> adS = log (adu(0,i)) - dim_/2.0 * log (adT) - dim_/2.0 * log (M_PI) - dim_/2.0;
        AutoDiff<1,SIMD<double>> adrS = adu(0,i) * adS;
        Vec<D,AutoDiff<1,SIMD<double>> > adF = adS*ucol.Range(1,D+1);
        adrS -= InnerProduct(grad.Col(i),adF);
        dEdt(0,i) = adrS.DValue(0);
        for(size_t j : Range(D))
          F(j,i) = adF(j).Value();
      }
  }
  
  void NumEntropyFlux (FlatMatrix<SIMD<double>> ula, FlatMatrix<SIMD<double>> ura,
		       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> flux) const
  {
    for(size_t i : Range(ula.Width()))
      {        
        auto rhol = ula(0,i);
        Vec<D,SIMD<double>> ul = ula.Col(i).Range(1,D+1);
        auto Tl = 2.0*(2.0*rhol*ula(D+1,i) - L2Norm2(ul)) / (dim_*rhol*rhol);
        Tl = IfPos(Tl,Tl,SIMD<double>(1e-10));
        auto Sl = (log(rhol) - dim_/2.0 * log(Tl) - dim_/2.0 * log (M_PI) - dim_/2.0);
        
        auto rhor = ura(0,i);
        Vec<D,SIMD<double>> ur = ura.Col(i).Range(1,D+1);
        auto Tr = 2.0*(2.0*rhor*ura(D+1,i) - L2Norm2(ur)) / (dim_*rhor*rhor);
        Tr = IfPos(Tr,Tr,SIMD<double>(1e-10));
        auto Sr = (log(rhor) - dim_/2.0 * log(Tr) - dim_/2.0 * log (M_PI) - dim_/2.0);

        Vec<D,SIMD<double>> n = normals.Col(i);
        auto uln = InnerProduct (ul, n);
        auto urn = InnerProduct (ur, n);
    
        auto upwind = IfPos(uln, uln * Sl, urn * Sr);
        flux.Col(i) = upwind;
      }
  }

  void CalcViscCoeffEl(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<SIMD<double>> elu_ipts,
                       FlatMatrix<SIMD<double>> res_ipts,
                       const double hi, double & coeff) const
  {
    double betaeli = 0.0;
    double rhoeli = 0.0;
    double visci = 0.0;

    Vec<D+2,SIMD<double>> eluk;

    for(size_t k : Range(elu_ipts.Width()))
      {
	eluk = elu_ipts.Col(k);
        for(size_t l : Range(SIMD<double>::Size()))
          {
            if(eluk(0)[l] > rhoeli) rhoeli = eluk(0)[l];
            if(fabs(res_ipts(0,k)[l]) > visci) visci = fabs(res_ipts(0,k)[l]);
          }

        SIMD<double> normu = 0.0;
        for(size_t l : Range(D))
          normu += eluk(l+1)*eluk(l+1);
        normu /= sqr(eluk(0));

        auto e = eluk(D+1) / eluk(0) - 0.5 * normu;
	auto elT = 4.0 * e / dim_;
        
	auto betaipt = sqrt(normu) + sqrt(gamma_*elT);
	for(size_t l : Range(SIMD<double>::Size()))
          if(betaipt[l] > betaeli) betaeli = betaipt[l];
      }
    // coeff = min (1.0/8*betaeli*rhoeli, 1.0/8*visci * hi);
    coeff = min (0.05*betaeli*rhoeli, 5*visci * hi);
    coeff *= hi;
  }
  
  void u2pUT(FlatMatrix<double> u, FlatMatrix<double> pUT, LocalHeap & lh) const
  {
    int ne = ma->GetNE();
    
    for (int i = 0; i < ne; i++)
      {
	HeapReset hr(lh);
	ElementId ei(VOL,i);

	const ScalarFiniteElement<D> & fel = dynamic_cast<const ScalarFiniteElement<D>&> (fes->GetFE (ei, lh));
	const IntegrationRuleTP<D> ir(ma->GetTrafo(ei, lh), 2*fel.Order()); // , &lh);
	MappedIntegrationRule<D,D> mir(ir,ma->GetTrafo(ei,lh),lh);
	IntRange dn = fes->GetElementDofs (i);
	  
	int ndof = fel.GetNDof();
	int nipt = ir.Size();
	  
	FlatVector<> elp(nipt,lh), elT(nipt, lh), elmach(nipt, lh), elpUT_col(nipt, lh);
	FlatVector<> elpUT_coef(ndof,lh);
	  
	  
	FlatMatrix<> elu_coef(D+2,ndof,lh);
	FlatMatrix<> elu_ipts(D+2,nipt,lh), elU(D,nipt,lh), elpUT_ipts(D+3,nipt,lh);
	  
	elu_coef = Trans(u.Rows(dn));
	  
	for(int k = 0; k < D+2; k++)
	  fel.Evaluate (ir, elu_coef.Row(k), elu_ipts.Row(k));
	  
	FlatVector<> eluk(D+2,lh);
	  
	for (int k = 0; k < nipt; k++) 
	  {
	    double fac = mir[k].GetWeight();
	    eluk = elu_ipts.Col(k);
	    double rho = eluk(0);
	      
	    for(int l = 0; l < D; l++)
	      elU(l,k) = eluk(l+1) / rho;
	      
	    double normU2 = L2Norm2(elU.Col(k));
	      
	    elU.Col(k) *= fac;
	      
	    double e = eluk(D+1) / rho - 0.5 * normU2;
	      
	    elT(k) = fac * 4.0 * e / dim_;
	    elp(k) = 0.5 * rho *elT(k) ;  // fac in elT(k) included
	      
	    elmach(k) = fac * sqrt(2 * normU2 / (gamma_ * elT(k) / fac) );
	      
	  }
	  
	elpUT_ipts.Row(0) = elp;
	elpUT_ipts.Rows(1,D+1) = elU;
	elpUT_ipts.Row(D+1) = elT;
	elpUT_ipts.Row(D+2) = elmach;
	  
	  
	for(int k = 0; k < D+3; k++)
	  {
	    fel.EvaluateTrans(ir,elpUT_ipts.Row(k),elpUT_coef);
	    pUT.Col(k).Range(dn) = elpUT_coef;
	  }
	  
	// SolveM(i,pUT.Rows(dn),lh);
	BASE::SolveM(i,pUT.Rows(dn),lh);
      }
  }
  
  
  void UpdateVisualisation(const BaseVector & hu, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    int ndof_total = gfrho->GetVector().FVDouble().Size();
    
    FlatMatrixFixWidth<D+2> u (ndof_total, &hu.FV<double>()(0));
    FlatMatrixFixWidth<D+3> pUT (ndof_total,lh);
    
    u2pUT(u, pUT, lh);
    
    FlatVector<> vecU(2*ndof_total, lh);
    for(int l = 0; l < ndof_total; l++)
      for(int m = 0; m < D; m++)
	vecU(D*l+m) = pUT(l,m+1);
    
    gfrho->GetVector().FVDouble() = u.Col(0);
    gfp->GetVector().FVDouble() = pUT.Col(0);
    gfU->GetVector().FVDouble() = vecU;
    gfT->GetVector().FVDouble() = pUT.Col(D+1);
    gfmach->GetVector().FVDouble() = pUT.Col(D+2);
  }

  template <typename MIP=BaseMappedIntegrationPoint, typename TA, typename TB>
  void InverseMap(const MIP & mip, const TA & grad,
		  const TB & u) const
  {
    auto InnerProduct = [](auto a, auto b)
      {
        typename decltype(a)::TELEM res(0.0);
        for (size_t i : Range(a))
          res += a(i)*b(i);
        return res;
      };

    // version diploma thesis
    auto temp = u(0) - InnerProduct(u.Range(1,D+1),grad);
    auto temp2 = 2*u(D+1)*u(0) - InnerProduct(u.Range(1,D+1),u.Range(1,D+1));
    auto rhoe = temp2 /(temp + sqrt( sqr(temp) - 4*(dim_+1)/sqr(dim_)*InnerProduct(grad,grad)*temp2));
    
    auto rho = sqr(u(0))/(temp - 2*InnerProduct(grad,grad)*rhoe/dim_);
    Vec<D,typename TA::TELEM> m = rho/u(0) * (u.Range(1,D+1) + 2*rhoe/dim_ * grad);
    auto E = (rho*u(D+1)+2*rhoe*InnerProduct(m,grad)/dim_)/u(0);

    // version in diss
    // auto a1 = dim_/2.0*(u(0)-InnerProduct(u.Range(1,D+1),grad));
    // auto a2 = 2*u(0)*u(D+1)-InnerProduct(u.Range(1,D+1),u.Range(1,D+1));
    // auto normgrad = InnerProduct(grad,grad);
    // auto p = a2/(a1 + sqrt(sqr(a1) - (dim_+1)*normgrad*a2));
    // auto rho = sqr(u(0))/(u(0) - (InnerProduct(u.Range(1,D+1),grad) + p*normgrad));
    // Vec<D,typename TA::TELEM> m = rho/u(0) * (u.Range(1,D+1) + p* grad);
    // auto E = (rho*u(D+1) + p*InnerProduct(m,grad))/u(0);

    u(0) = rho;
    u.Range(1,D+1) = m;
    u(D+1) = E;
  }

  template <typename T>
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<T> grad, FlatMatrix<T> u) const
  {
    for (auto i : Range(mir))
      InverseMap(mir[i], grad.Col(i), u.Col(i));
  }

};

/////////////////////////////////////////////////////////////////////////

shared_ptr<ConservationLaw> CreateEuler(const shared_ptr<GridFunction> & gfu,
					const shared_ptr<TentPitchedSlab> & tps)
{
  int dim = tps->ma->GetDimension();
  switch(dim){
  case 1:
    return make_shared<Euler<1>>(gfu, tps);
  case 2:
    return make_shared<Euler<2>>(gfu, tps);
  }
  throw Exception ("euler only available for 1D and 2D");
}
