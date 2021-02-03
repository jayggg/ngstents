#ifndef CONSERVATIONLAW_HPP
#define CONSERVATIONLAW_HPP

#ifdef USE_NETGEN_GUI
#include "myvisual.hpp"
#endif
#include "tents.hpp"
#include "tentsolver.hpp"
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
  shared_ptr<GridFunction> gfuorig = nullptr;
  shared_ptr<GridFunction> gfnu = nullptr;

  shared_ptr<LocalHeap> pylh = nullptr;
  shared_ptr<BaseVector> u = nullptr;     // u(n)
  shared_ptr<BaseVector> uinit = nullptr; // initial data, also used for bc

  shared_ptr<TentSolver> tentsolver;
public:
  ConservationLaw (const shared_ptr<GridFunction> & agfu,
		   const shared_ptr<TentPitchedSlab> & atps,
		   const string & eqn)
    : gfu{agfu}, fes{agfu->GetFESpace()}, tps {atps}, ma {atps->ma},
      equation {eqn}, order{agfu->GetFESpace()->GetOrder()}
  { };
  
  virtual ~ConservationLaw() { ; }
  
  virtual void SetBC() = 0;

  virtual void SetFluxField(shared_ptr<CoefficientFunction> cf) = 0;

  virtual void SetMaterialParameters(shared_ptr<CoefficientFunction> cf_mu,
                                     shared_ptr<CoefficientFunction> cf_eps) = 0;

  virtual void SetTentSolver(string method, int stages, int substeps) = 0;

  virtual void Propagate(LocalHeap & lh) = 0;

  virtual void PropagateSAT(int stages, int substeps,
			    BaseVector & hu, BaseVector & hu_init,
			    LocalHeap & lh) = 0;

  virtual void PropagateSARK(int stages, int substeps,
			     BaseVector & hu, BaseVector & hu_init,
			     LocalHeap & lh) = 0;
  
};


/*Template for conservation law classes.
Any conservation law to be used in the MTP scheme should derive
from this class.
The template parameters are as follows:
EQUATION: the class representing the equation
DIM: spatial dimension in which the conservation law is prescribed upon
COMP: number of state variables/equations
ECOMP: number of state variables/equations for entropy residual (non-linear eqs)
*/
template <typename EQUATION, int DIM, int COMP, int ECOMP>
class T_ConservationLaw : public ConservationLaw,
			  public enable_shared_from_this<T_ConservationLaw<EQUATION,DIM,COMP,ECOMP>>
{
protected:
  FlatVector<> nu;  // viscosity coefficient

  bool def_bcnr = false; // check if the array below is properly set
  int maxbcnr = 4;
  Array<int> bcnr; // array of boundary condition numbers

  // collection of tents in timeslab
  Table<int> & tent_dependency = tps->tent_dependency;

  const EQUATION & Cast() const {return static_cast<const EQUATION&> (*this);}

public:

  // advancing front (used for time-dependent bc)
  shared_ptr<GridFunction> gftau = nullptr;

  T_ConservationLaw (const shared_ptr<GridFunction> & gfu,
		     const shared_ptr<TentPitchedSlab> & tps,
		     const string & eqn)
    : ConservationLaw(gfu, tps, eqn)
  {
    // TODO: do I need that
    size_t heapsize = 10*1000000;
    pylh = make_shared<LocalHeap>(heapsize,"ConsLaw - py main heap",true);

    // TODO: set boundaries later
    // store boundary condition numbers
    bcnr = FlatArray<int>(ma->GetNFacets(),*pylh);
    bcnr = -1;

    // Main L2 finite element space based on spatial mesh
    // Flags fesflags = Flags();
    // fesflags.SetFlag("order",order);
    // fesflags.SetFlag("dim",COMP);
    // fesflags.SetFlag("all_dofs_together");
    // fes = dynamic_pointer_cast<L2HighOrderFESpace>(
    //     CreateFESpace("l2ho", ma, fesflags));
    // fes->Update();
    // fes->FinalizeUpdate();
    // 
    // gfu = CreateGridFunction(fes,"u",Flags());
    // gfu->Update();

    // if (ECOMP > 0)
    //   {
    //     // Scalar L2 finite element space for entropy residual
    //     shared_ptr<FESpace> fes_scal =
    //       CreateFESpace("l2ho", ma,
    //           Flags().SetFlag("order",order).SetFlag("all_dofs_together"));
    //     fes_scal->Update();
    //     fes_scal->FinalizeUpdate();
    //     gfres = CreateGridFunction(fes_scal,"res",Flags());
    // 	gfres->Update();
    // 
    //     // Zero order L2 finite element space for viscosity
    // 	shared_ptr<FESpace> fes_lo = CreateFESpace("l2ho", ma,
    //                                                Flags().SetFlag("order",0));
    // 	fes_lo->Update();
    // 	fes_lo->FinalizeUpdate();
    // 	gfnu = CreateGridFunction(fes_lo,"nu",Flags());
    // 	gfnu->Update();
    //   }

    // first order H1 space for the advancing front
    shared_ptr<FESpace> fesh1 = CreateFESpace("h1ho", ma,
                                              Flags().SetFlag("order",1));
    fesh1->Update();
    fesh1->FinalizeUpdate();
    gftau = CreateGridFunction(fesh1,"tau",Flags().SetFlag("novisual"));
    gftau->Update();
    gftau->GetVector() = 0.0;

    AllocateVectors();
  }

  void AllocateVectors()
  {
    u = gfu->GetVectorPtr(); // TODO: should not be needed
    uinit = u->CreateVector(); // TODO: how to set inflow boundary?
    // if(gfnu != NULL)
    //   {
    // 	gfnu->Update();
    // 	nu.AssignMemory(gfnu->GetVector().FVDouble().Size(),
    //                     &gfnu->GetVector().FVDouble()(0));
    // 	nu = 0.0;
    //   }
  }

  virtual ~T_ConservationLaw() { ; }
  
  // Set the boundary condition numbers from the mesh boundary elements indices
  // These indices are 0-based here and 1-based in Python. 
  //  0: outflow, 1: wall, 2: inflow, 3: transparent 
  void SetBC()
  {
    if(!def_bcnr)
      for(int i : Range(ma->GetNSE()))
        {
          auto sel = ElementId(BND,i);
          auto fnums = ma->GetElFacets(sel);
          bcnr[fnums[0]] = ma->GetElIndex(sel);
        }
  }

  virtual void SetFluxField(shared_ptr<CoefficientFunction> cf)
  {
    throw Exception("SetFluxField just available for Advection equation");
  }

  virtual void SetMaterialParameters(shared_ptr<CoefficientFunction> cf_mu,
                                     shared_ptr<CoefficientFunction> cf_eps)
  {
    throw Exception("SetMaterialParameters just available for Wave equation");
  }
  
  template <int W>
  void SolveM (const Tent & tent, int loci, FlatMatrixFixWidth<W> mat,
               LocalHeap & lh) const
  {
    // TODO: check space
    auto fedata = tent.fedata;
    if (!fedata)
        throw Exception ("Expected tent.fedata to be set!");

    HeapReset hr(lh);
    const DGFiniteElement<DIM> & fel =
      static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[loci]);

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
    // TODO: check space
    auto fedata = tent.fedata;
    if (!fedata)
        throw Exception ("Expected tent.fedata to be set!");

    HeapReset hr(lh);
    const DGFiniteElement<DIM> & fel =
      static_cast<const DGFiniteElement<DIM>&> (*fedata->fei[loci]);

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
		     + ToString(this->equation) + "equation!");
  }

  void u_transparent(const SIMD_BaseMappedIntegrationRule & mir,
                     FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                     FlatMatrix<SIMD<double>> u_transp) const
  {
    throw Exception ("Transparent boundary just available for wave equation!");
  }

  void CalcFluxTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                     FlatMatrixFixWidth<COMP> u0, FlatMatrixFixWidth<COMP> flux,
                     double tstar, LocalHeap & lh);

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

  [[deprecated]]
  void EntropyFlux (const Vec<DIM+2> & ml, const Vec<DIM+2> & mr,
                    const Vec<DIM> & n, double & flux) const
  {
    cout << "no overload for EntropyFlux" << endl;
  }

  // numerical flux for the entropy flux
  void EntropyFlux (FlatMatrix<SIMD<double>> ml, FlatMatrix<SIMD<double>> mr,
                    FlatMatrix<SIMD<double>> n,
                    FlatMatrix<SIMD<double>> flux) const
  {
    cout << "no overload for EntropyFlux for FlatMatrix<SIMD>" << endl;
  }

  // apply viscosity
  void CalcViscosityTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                          FlatMatrixFixWidth<COMP> ubnd, FlatVector<double> nu,
                          FlatMatrixFixWidth<COMP> visc, LocalHeap & lh);

  // calculate entropy residual on a tent
  void CalcEntropyResidualTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                                FlatMatrixFixWidth<COMP> ut,
                                FlatMatrixFixWidth<ECOMP> res,
                                FlatMatrixFixWidth<COMP> u0, double tstar,
                                LocalHeap & lh);

  // calculate viscosity coefficient based on the entropy residual on a tent
  double CalcViscosityCoefficientTent (int tentnr, FlatMatrixFixWidth<COMP> u,
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

  void Cyl2Tent (int tentnr, double tstar,
		 FlatMatrixFixWidth<COMP> uhat, FlatMatrixFixWidth<COMP> u,
		 LocalHeap & lh);

  void ApplyM1 (int tentnr, double tstar, FlatMatrixFixWidth<COMP> u,
                FlatMatrixFixWidth<COMP> res, LocalHeap & lh);

  void Tent2Cyl (int tentnr, double tstar,
		 FlatMatrixFixWidth<COMP> u, FlatMatrixFixWidth<COMP> uhat,
                 bool solvemass, LocalHeap & lh);
  
  ////////////////////////////////////////////////////////////////
  // time stepping methods 
  ////////////////////////////////////////////////////////////////

  void SetTentSolver(string method, int stages, int substeps)
  {
    tentsolver = make_shared<SAT<EQUATION,DIM,COMP,ECOMP>>(this->shared_from_this(), stages, substeps);
  }
  
  void Propagate(LocalHeap & lh);

  void PropagateSAT(int stages, int substeps,
		    BaseVector & hu, BaseVector & hu_init,
		    LocalHeap & lh);

  void PropagateSARK(int stages, int substeps,
		     BaseVector & hu, BaseVector & hu_init,
		     LocalHeap & lh);

};


#endif // CONSERVATIONLAW_HPP
