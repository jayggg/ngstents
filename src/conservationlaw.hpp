#ifndef CONSERVATIONLAW_HPP
#define CONSERVATIONLAW_HPP

#ifdef USE_NETGEN_GUI
#include "myvisual.hpp"
#endif
#include "tents.hpp"
#include "tentsolver.hpp"
#include "vis3d.hpp"
#include <atomic>

class ConservationLaw
{
public:
  shared_ptr<MeshAccess> ma = nullptr;
  shared_ptr<TentPitchedSlab> tps = nullptr;

  const int order = {};
  const string equation = {};
  shared_ptr<FESpace> fes = nullptr;
  shared_ptr<GridFunction> gfu = nullptr;
  shared_ptr<GridFunction> gfres = nullptr;
  shared_ptr<GridFunction> gfnu = nullptr;

  shared_ptr<LocalHeap> pylh = nullptr;

  shared_ptr<BaseVector> u = nullptr;     // u(n)
  shared_ptr<BaseVector> uinit = nullptr; // initial data, also used for bc

  shared_ptr<TentSolver> tentsolver;

  shared_ptr<GridFunction> gftau = nullptr;  // advancing front (used for time-dependent bc)
  shared_ptr<CoefficientFunction> cftau = nullptr;  // CF representing gftau

  shared_ptr<ProxyFunction> proxy_u = nullptr;
  shared_ptr<ProxyFunction> proxy_uother = nullptr;

  // instance for 3D visualization of tent slab solutions based on 2D meshes.
  shared_ptr<Visualization3D> vis3d = nullptr;

  shared_ptr<ProxyFunction> proxy_graddelta = nullptr;
  shared_ptr<ProxyFunction> proxy_res = nullptr;
public:
  ConservationLaw (const shared_ptr<GridFunction> & agfu,
		   const shared_ptr<TentPitchedSlab> & atps,
		   const string & eqn)
    : gfu{agfu}, fes{agfu->GetFESpace()}, tps {atps}, ma {atps->ma},
      equation {eqn}, order{agfu->GetFESpace()->GetOrder()}
  { };
  
  virtual ~ConservationLaw() { ; }
  
  virtual void SetBC(int bcnr, const BitArray & region) = 0;

  virtual void CheckBC() = 0;

  virtual void SetBoundaryCF(int bcnr, shared_ptr<CoefficientFunction> cf) = 0;
  
  virtual void SetVectorField(shared_ptr<CoefficientFunction> cf) = 0;

  virtual void SetMaterialParameters(shared_ptr<CoefficientFunction> cf_mu,
                                     shared_ptr<CoefficientFunction> cf_eps) = 0;

  virtual void SetViscosityCoefficient(shared_ptr<CoefficientFunction> cf_visc) = 0;

  virtual void SetNumEntropyFlux(shared_ptr<CoefficientFunction> cf_numentropyflux) = 0;

  virtual void SetTentSolver(string method, int stages, int substeps) = 0;

  // virtual void Propagate(LocalHeap & lh) = 0;

  virtual void Propagate(LocalHeap & lh, shared_ptr<GridFunction> hdgf) = 0;

};


////////////////////////////////////////////////////////////////////////////
///
/// TEMPLATE FOR CONSERVATION LAW SPECIFICATION IN MTP SCHEMES.
///
/// Any conservation law to be used in the MTP scheme should derive
/// from this class. The template parameters are as follows: 
///   EQUATION: the class representing the equation
///   DIM: spatial dimension in which the conservation law is prescribed upon
///   COMP: number of state variables/equations
///   ECOMP: number of state variables/equations for entropy residual (non-linear eqs)
///
/// The equation is
///
///      d_t u + div_x F(u) = 0,     on Omega, 
///
/// where
///
///   -  Omega = An DIM-space dimensional domain,
///   -      u = Solution, an COMP x 1 vector function,
///   -      F = Flux, an COMP x DIM matrix function,
///   -    d_t = time derivative,
///   -   div_x = row-wise spatial divergence.
///

template <typename EQUATION, int DIM, int COMP, int ECOMP, bool SYMBOLIC=false>
class T_ConservationLaw :
  public ConservationLaw,
  public enable_shared_from_this<T_ConservationLaw<EQUATION,DIM,COMP,ECOMP,SYMBOLIC>>
{
protected:

  bool def_bcnr = false; ///< check if the array below is properly set
  Array<int> bcnr;       ///< array of boundary condition numbers
  Array<shared_ptr<CoefficientFunction>> cf_bnd; ///< CF used for boundary values
  bool cf_bnd_deriv = false;
  BitArray scale_deriv; ///< scale time-dependent boundary CF by tent height (for SAT)
  FlatVector<> nu;  ///< viscosity coefficient (for nonlinear case)
  shared_ptr<CoefficientFunction> cf_numentropyflux = nullptr;
  /// collection of tents in timeslab
  Table<int> & tent_dependency = tps->tent_dependency;

  const EQUATION & Cast() const {return static_cast<const EQUATION&> (*this);}

public:
  enum { NCOMP = COMP };
  enum { NECOMP = ECOMP };

  T_ConservationLaw (const shared_ptr<GridFunction> & gfu,
		     const shared_ptr<TentPitchedSlab> & tps,
		     const string & eqn)
    : ConservationLaw(gfu, tps, eqn)
  {
    size_t heapsize = 10*1000000;
    pylh = make_shared<LocalHeap>(heapsize,"ConsLaw - py main heap",true);

    // store boundary condition numbers
    bcnr = FlatArray<int>(ma->GetNFacets(),*pylh);
    bcnr = -1;
    scale_deriv = BitArray(ma->GetNRegions(BND), *pylh);
    scale_deriv.Clear();

    // check dimension of space
    shared_ptr<L2HighOrderFESpace> fes_check = dynamic_pointer_cast<L2HighOrderFESpace>(fes);
    if(fes_check && (fes->GetDimension() != COMP) )
      throw Exception("set dimension of L2 finite element space to "
		      +ToString(COMP)+" by adding the argument 'dim="+ToString(COMP)+"'");

    u = gfu->GetVectorPtr();
    uinit = u->CreateVector();

    if (ECOMP > 0)
      {
        // Scalar L2 finite element space for entropy residual
        shared_ptr<FESpace> fes_scal =
          CreateFESpace("l2ho", ma,
			Flags().SetFlag("order",order).SetFlag("all_dofs_together"));
        fes_scal->Update();
        fes_scal->FinalizeUpdate();
        gfres = CreateGridFunction(fes_scal,"res",Flags());
    	gfres->Update();

	// residual proxy
	proxy_res = fes_scal->GetProxyFunction(false);
	// graddelta proxy
	proxy_graddelta = make_shared<ProxyFunction>(fes_scal, false, false,
						     fes_scal->GetFluxEvaluator(VOL),
						     nullptr, nullptr, nullptr, nullptr, nullptr);

        // Zero order L2 finite element space for viscosity
    	shared_ptr<FESpace> fes_lo = CreateFESpace("l2ho", ma,
                                                   Flags().SetFlag("order",0));
    	fes_lo->Update();
    	fes_lo->FinalizeUpdate();
    	gfnu = CreateGridFunction(fes_lo,"nu",Flags());
    	gfnu->Update();
	nu.AssignMemory(gfnu->GetVector().FVDouble().Size(),
                        &gfnu->GetVector().FVDouble()(0));
    	nu = 0.0;
      }

    // first order H1 space for the advancing front
    shared_ptr<FESpace> fesh1 = CreateFESpace("h1ho", ma,
                                              Flags().SetFlag("order",1));
    fesh1->Update();
    fesh1->FinalizeUpdate();
    gftau = CreateGridFunction(fesh1,"tau",Flags().SetFlag("novisual"));
    gftau->Update();
    gftau->GetVector() = 0.0;
    cftau = make_shared<GridFunctionCoefficientFunction>(gftau);
  }

  virtual ~T_ConservationLaw() { ; }

  // set boundary condition numbers for given region
  void SetBC(int bc, const BitArray & region)
  {
    def_bcnr = true;
    for(auto i : Range(ma->GetNSE()))
      {
        auto sel = ElementId(BND,i);
        auto fnums = ma->GetElFacets(sel);
        if(region.Test(ma->GetElIndex(sel)))
          bcnr[fnums[0]] = bc;
      }
  }

  // Set old style boundary condition numbers from the mesh boundary elements indices.
  // These indices are 0-based here and 1-based in Python. 
  // 0: outflow, 1: wall, 2: inflow, 3: transparent 
  void CheckBC()
  {
    if(!def_bcnr)
      for(auto i : Range(ma->GetNSE()))
        {
          auto sel = ElementId(BND,i);
          auto fnums = ma->GetElFacets(sel);
          bcnr[fnums[0]] = ma->GetElIndex(sel);
        }
  }

  void SetBoundaryCF(int bcnr, shared_ptr<CoefficientFunction> cf)
  {
    if(cf_bnd.Size()==0)
      cf_bnd.Append(cf);
    else
      cout << "Resetting boundary values (discarding prior set values)"
	   << endl;
  }


  // derive boundary coefficient functions
  // fill table with 1/j! * cf^(j)
  void DeriveBoundaryCF(int stages)
  {
    if(cf_bnd_deriv)
      return;

    if(cf_bnd.Size() == 0)
      return;

    // check if cf_bnd is time-dependent
    auto dtau_cf = cf_bnd[0]->Diff(cftau.get(),
				   make_shared<ConstantCoefficientFunction>(1.0));
    auto components = dtau_cf->InputCoefficientFunctions();
    
    for( auto i : Range(components.Size()))
      if(components[i])
	if( components[i]->GetDescription() != "ZeroCF")
	  scale_deriv.SetBit(i);

    for(size_t i : Range(stages-1))
      {
	auto dtau =
	  cf_bnd[i]->Diff(cftau.get(),
			  make_shared<ConstantCoefficientFunction>(1.0/(i+1)));
	auto temp = cf_bnd[i]->Diff(proxy_u.get(), proxy_u); // identity
	cf_bnd.Append(dtau + temp);
      }

    cf_bnd_deriv = true;
  }

  virtual void SetVectorField(shared_ptr<CoefficientFunction> cf)
  {
    throw Exception("SetVectorField just available for Advection equation");
  }

  virtual void SetMaterialParameters(shared_ptr<CoefficientFunction> cf_mu,
                                     shared_ptr<CoefficientFunction> cf_eps)
  {
    throw Exception("SetMaterialParameters just available for Wave equation");
  }

  virtual void SetViscosityCoefficient(shared_ptr<CoefficientFunction> cf_visc)
  {
    throw Exception("SetViscosityCoefficient just available for SymbolicConsLaw");
  }

  virtual void SetNumEntropyFlux(shared_ptr<CoefficientFunction> cf_numentropyflux)
  {
    throw Exception("SetNumEntropyFlux just available for SymbolicConsLaw");
  }

  template <int W>
  void SolveM (const Tent & tent, int loci, FlatMatrixFixWidth<W> mat,
               LocalHeap & lh) const
  {
    auto fedata = tent.fedata;
    if (!fedata)
        throw Exception ("Expected tent.fedata to be set!");

    HeapReset hr(lh);
    auto & fel = static_cast<const BaseScalarFiniteElement&> (*fedata->fei[loci]);

    bool curved = ma->GetElement(ElementId(VOL,tent.els[loci])).is_curved;
    if (curved)
      {
	FlatVector<> diagmass(mat.Height(),lh);
	fel.GetDiagMassMatrix(diagmass);

	SIMD_IntegrationRule & ir = *fedata->iri[loci];
	SIMD_BaseMappedIntegrationRule & mir = *fedata->miri[loci];
	FlatMatrix<SIMD<double>> pntvals(W, ir.Size(), lh);

	for (int i : Range(mat.Height()))
	  mat.Row(i) /= diagmass(i);
        fel.Evaluate(ir, mat, pntvals);
	for (int comp : Range(W))
	  {
            for(int i : Range(ir))
              pntvals(comp,i) *= ir[i].Weight() / mir[i].GetMeasure();
	  }
        mat = 0.0;
        fel.AddTrans(ir,pntvals, mat);
	for (int i : Range(mat.Height()))
	  mat.Row(i) /= diagmass(i);
      }

    else
      {
        FlatVector<> diagmass(mat.Height(),lh);

        double measure = (*fedata->miri[loci])[0].GetMeasure()[0];
        fel.GetDiagMassMatrix(diagmass);
        for (size_t j = 0; j < diagmass.Size(); j++)
          diagmass(j) = 1.0 / (diagmass(j) * measure);
        for(size_t j = 0; j < mat.Height(); j++)
          mat.Row (j) *= diagmass(j);
      }
  }


  template <int W>
  void SolveM (const Tent & tent, int loci,
               FlatVector<SIMD<double>> delta,
               FlatMatrixFixWidth<W> mat, LocalHeap & lh) const
  {
    auto fedata = tent.fedata;
    if (!fedata)
        throw Exception ("Expected tent.fedata to be set!");

    HeapReset hr(lh);
    auto & fel = static_cast<const BaseScalarFiniteElement&> (*fedata->fei[loci]);

    FlatVector<> diagmass(mat.Height(),lh);
    fel.GetDiagMassMatrix(diagmass);

    SIMD_IntegrationRule & ir = *fedata->iri[loci];
    SIMD_BaseMappedIntegrationRule & mir = *fedata->miri[loci];
    FlatMatrix<SIMD<double>> pntvals(W, ir.Size(), lh);

    for (int i : Range(mat.Height()))
      mat.Row(i) /= diagmass(i);
    fel.Evaluate(ir, mat, pntvals);
    for (int comp : Range(W))
      {
        for(int i : Range(ir))
          pntvals(comp,i) *= ir[i].Weight() *
            delta(i) / mir[i].GetMeasure(); //scale with 1/delta
      }
    mat = 0.0;
    fel.AddTrans(ir,pntvals, mat);
    for (int i : Range(mat.Height()))
      mat.Row(i) /= diagmass(i);
  }

  template <typename SCAL>
  Mat<COMP,DIM,SCAL> Flux (const BaseMappedIntegrationPoint & mip,
                           const FlatVec<COMP,SCAL> & u) const
  {
    throw Exception ("flux not implemented");
  }

  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    throw Exception ("flux for FlatMatrix<SIMD> not implemented");
  }
  
  void NumFlux(const SIMD_BaseMappedIntegrationRule & mir,
	       FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
	       FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    throw Exception ("numerical flux for FlatMatrix<SIMD> not implemented");
  }

  Vec<COMP> NumFlux (const BaseMappedIntegrationPoint & mip,
		     const FlatVec<COMP> & ul, const FlatVec<COMP> & ur,
		     const Vec<DIM> & nv) const
  {
    throw Exception ("numerical flux not implemented");
  }

  void u_reflect(const SIMD_BaseMappedIntegrationRule & mir,
		 FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const
  {
    throw Exception ("reflecting boundary conditions not implemented for "
		     + ToString(this->equation) + " equation!");
  }

  void u_transparent(const SIMD_BaseMappedIntegrationRule & mir,
                     FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                     FlatMatrix<SIMD<double>> u_transp) const
  {
    throw Exception ("Transparent boundary just available for wave equation!");
  }

  void CalcFluxTent(const Tent & tent, const FlatMatrixFixWidth<COMP> u,
		    FlatMatrixFixWidth<COMP> u0, FlatMatrixFixWidth<COMP> flux,
		    double tstar, int derive_cf_bnd, LocalHeap & lh);

  ////////////////////////////////////////////////////////////////
  // entropy viscosity for nonlinear conservation laws
  ////////////////////////////////////////////////////////////////

  // evaluation of the temporal derivative of the entropy E
  // and evaluation the entropy flux F
  void CalcEntropy(FlatMatrix<AutoDiff<1,SIMD<double>>> adu,
                   FlatMatrix<AutoDiff<1,SIMD<double>>> grad,
		   FlatMatrix<SIMD<double>> dEdt,
                   FlatMatrix<SIMD<double>> F) const
  {
    cout << "no overload for CalcEntropy for tent pitching" << endl;
  }

  // numerical flux for the entropy flux
  void NumEntropyFlux (FlatMatrix<SIMD<double>> ml, FlatMatrix<SIMD<double>> mr,
		       FlatMatrix<SIMD<double>> n,
		       FlatMatrix<SIMD<double>> flux) const
  {
    cout << "no overload for NumEntropyFlux for FlatMatrix<SIMD>" << endl;
  }

  // apply viscosity
  void CalcViscosityTent (const Tent & tent, FlatMatrixFixWidth<COMP> u,
                          FlatMatrixFixWidth<COMP> ubnd, FlatVector<double> nu,
                          FlatMatrixFixWidth<COMP> visc, LocalHeap & lh);

  // calculate entropy residual on a tent
  void CalcEntropyResidualTent (const Tent & tent, FlatMatrixFixWidth<COMP> u,
                                FlatMatrixFixWidth<COMP> ut,
                                FlatMatrixFixWidth<ECOMP> res,
                                FlatMatrixFixWidth<COMP> u0, double tstar,
                                LocalHeap & lh);

  // calculate viscosity coefficient based on the entropy residual on a tent
  double CalcViscosityCoefficientTent (const Tent & tent, FlatMatrixFixWidth<COMP> u,
                                       FlatMatrixFixWidth<ECOMP> hres,
				       double tstar, LocalHeap & lh);

  // calculate viscosity coefficient based on the entropy residual on an element
  void CalcViscCoeffEl(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<SIMD<double>> elu_ipts,
                       FlatMatrix<SIMD<double>> res_ipts,
                       const double hi, double & coeff) const
  {
    cout << "no overload for CalcViscCoeffEl with FlatMatrix<SIMD>" << endl;
  }

  ////////////////////////////////////////////////////////////////
  // maps 
  ////////////////////////////////////////////////////////////////

  template <typename T = SIMD<double>>
  void InverseMap(const SIMD_BaseMappedIntegrationRule & mir,
		  FlatMatrix<T> grad, FlatMatrix<T> u) const
  {
    throw Exception ("TransformBack for FlatMatrix<SIMD> not available");
  }

  void Cyl2Tent (const Tent & tent, double tstar,
		 const FlatMatrixFixWidth<COMP> uhat, FlatMatrixFixWidth<COMP> u,
		 LocalHeap & lh);

  void ApplyM1 (const Tent & tent, double tstar,
		FlatMatrixFixWidth<COMP> u, FlatMatrixFixWidth<COMP> res,
		LocalHeap & lh);

  void Tent2Cyl (const Tent & tent, double tstar,
		 FlatMatrixFixWidth<COMP> u, FlatMatrixFixWidth<COMP> uhat,
                 bool solvemass, LocalHeap & lh);
  
  ////////////////////////////////////////////////////////////////
  // time stepping methods 
  ////////////////////////////////////////////////////////////////

  void SetTentSolver(string method, int stages, int substeps)
  {
    if(method == "SAT")
      tentsolver = make_shared<SAT<T_ConservationLaw<EQUATION,DIM,COMP,ECOMP,SYMBOLIC>>>
	(this->shared_from_this(), stages, substeps);
    else if(method == "SARK")
      tentsolver = make_shared<SARK<T_ConservationLaw<EQUATION,DIM,COMP,ECOMP,SYMBOLIC>>>
	(this->shared_from_this(), stages, substeps);
    else
      throw Exception("unknown TentSolver "+method);
  }
  
  void Propagate(LocalHeap & lh, shared_ptr<GridFunction> hdgf);

};


#endif // CONSERVATIONLAW_HPP
