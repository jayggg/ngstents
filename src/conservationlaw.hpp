#ifndef CONSERVATIONLAW_HPP
#define CONSERVATIONLAW_HPP

#ifdef USE_NETGEN_GUI
#include "myvisual.hpp"
#endif
#include "tents.hpp"
#include <atomic>
class ConservationLaw
{
public:
  shared_ptr<MeshAccess> ma = nullptr;
  shared_ptr<TentPitchedSlab> tps = nullptr;

  const int order = {};
  const string equation = {};
  shared_ptr<L2HighOrderFESpace> fes = nullptr;
  shared_ptr<GridFunction> gfu = nullptr;
  shared_ptr<GridFunction> gfres = nullptr;
  shared_ptr<GridFunction> gfuorig = nullptr;
  shared_ptr<GridFunction> gfnu = nullptr;

  shared_ptr<LocalHeap> pylh = nullptr;
  shared_ptr<BaseVector> u = nullptr;     // u(n)
  shared_ptr<BaseVector> uinit = nullptr; // initial data, also used for bc
public:
  ConservationLaw (const shared_ptr<TentPitchedSlab> & atps,
		   const string & eqn, int aorder)
    : tps {atps}, ma {atps->ma}, order {aorder}, equation {eqn}
  { };
  
  virtual ~ConservationLaw() { ; }
  
  virtual void SetBC() = 0;
  
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
XDEPENDENT: whether the flux depends on the spatial coordinates
*/
template <typename EQUATION, int DIM, int COMP, int ECOMP, bool XDEPENDENT>
class T_ConservationLaw : public ConservationLaw
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

  T_ConservationLaw (const shared_ptr<TentPitchedSlab> & tps,
		     const string & eqn, int order)
    : ConservationLaw(tps, eqn, order)
  {
    size_t heapsize = 10*1000000;
    pylh = make_shared<LocalHeap>(heapsize,"ConsLaw - py main heap",true);

    // store boundary condition numbers
    bcnr = FlatArray<int>(ma->GetNFacets(),*pylh);
    bcnr = -1;

    // Main L2 finite element space based on spatial mesh
    Flags fesflags = Flags();
    fesflags.SetFlag("order",order);
    fesflags.SetFlag("dim",COMP);
    fesflags.SetFlag("all_dofs_together");
    fes = dynamic_pointer_cast<L2HighOrderFESpace>(
        CreateFESpace("l2ho", ma, fesflags));
    fes->Update();
    fes->FinalizeUpdate();

    gfu = CreateGridFunction(fes,"u",Flags());
    gfu->Update();

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

        // Zero order L2 finite element space for viscosity
	shared_ptr<FESpace> fes_lo = CreateFESpace("l2ho", ma,
                                                   Flags().SetFlag("order",0));
	fes_lo->Update();
	fes_lo->FinalizeUpdate();
	gfnu = CreateGridFunction(fes_lo,"nu",Flags());
	gfnu->Update();
      }

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
    u = gfu->GetVectorPtr();
    uinit = u->CreateVector();
    if(gfnu != NULL)
      {
	gfnu->Update();
	nu.AssignMemory(gfnu->GetVector().FVDouble().Size(),
                        &gfnu->GetVector().FVDouble()(0));
	nu = 0.0;
      }
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

  template <int W>
  void SolveM (const Tent & tent, int loci, FlatMatrixFixWidth<W> mat,
               LocalHeap & lh) const
  {
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
    return Cast().Flux(u);
  }

  void Flux (const SIMD_BaseMappedIntegrationRule & mir,
             FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> flux) const
  {
    throw Exception ("flux for FlatMatrix<SIMD> not implemented");
  }

  void Flux(SIMD_BaseMappedIntegrationRule & mir,
            FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
            FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    if (!XDEPENDENT)
      Cast().Flux (ul, ur, normals, fna);
    else
      throw Exception ("simd-flux not implemented for X-dependent equation");
  }

  void Flux(FlatMatrix<SIMD<double>> ul, FlatMatrix<SIMD<double>> ur,
            FlatMatrix<SIMD<double>> normals, FlatMatrix<SIMD<double>> fna) const
  {
    throw Exception ("flux for FlatMatrix<SIMD> not implemented for boundary");
  }

  Vec<COMP> Flux (const BaseMappedIntegrationPoint & mip,
                  const FlatVec<COMP> & ul, const FlatVec<COMP> & ur,
                  const Vec<DIM> & nv) const
  {
    return Cast().Flux(ul,ur,nv);
  }

  template<typename SCAL=double>
  Vec<COMP,SCAL> Flux (const FlatVec<COMP,SCAL> & ul,
                       const FlatVec<COMP,SCAL> & ur,
                       const Vec<DIM,SCAL> & nv) const
  {
    throw Exception ("flux not implemented for boundary");
  }

  virtual void u_reflect(FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                 FlatMatrix<SIMD<double>> u_refl) const
  {
    Cast().u_reflect(u,normals,u_refl);
  }

  virtual void u_transparent(SIMD_BaseMappedIntegrationRule & mir,
                     FlatMatrix<SIMD<double>> u, FlatMatrix<SIMD<double>> normals,
                     FlatMatrix<SIMD<double>> u_transp) const
  {
    throw Exception ("Transparent boundary just available for wave equation!");
  }

  void EntropyFlux (const Vec<DIM+2> & ml, const Vec<DIM+2> & mr,
                    const Vec<DIM> & n, double & flux) const
  {
    cout << "no overload for EntropyFlux" << endl;
  }

  void EntropyFlux (FlatMatrix<SIMD<double>> ml, FlatMatrix<SIMD<double>> mr,
                    FlatMatrix<SIMD<double>> n,
                    FlatMatrix<SIMD<double>> flux) const
  {
    cout << "no overload for EntropyFlux for FlatMatrix<SIMD>" << endl;
  }

  void CalcFluxTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                     FlatMatrixFixWidth<COMP> u0, FlatMatrixFixWidth<COMP> flux,
                     double tstar, LocalHeap & lh);

  void Cyl2Tent (int tentnr, FlatMatrixFixWidth<COMP> uhat,
                          FlatMatrixFixWidth<COMP> u, double tstar,
                          LocalHeap & lh);

  void CalcViscosityTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                          FlatMatrixFixWidth<COMP> ubnd, FlatVector<double> nu,
                          FlatMatrixFixWidth<COMP> visc, LocalHeap & lh);

  void CalcEntropy(FlatMatrix<AutoDiff<1,SIMD<double>>> adu,
                   FlatMatrix<AutoDiff<1,SIMD<double>>> grad,
		   FlatMatrix<SIMD<double>> dEdt,
                   FlatMatrix<SIMD<double>> F) const
  {
    cout << "no overload for CalcEntropy for tent pitching" << endl;
  }

  void CalcEntropyResidualTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                                FlatMatrixFixWidth<COMP> ut,
                                FlatMatrixFixWidth<ECOMP> res,
                                FlatMatrixFixWidth<COMP> u0, double tstar,
                                LocalHeap & lh);

  double CalcViscosityCoefficientTent (int tentnr, FlatMatrixFixWidth<COMP> u,
                                       FlatMatrixFixWidth<ECOMP> hres,
				       double tstar, LocalHeap & lh);

  void CalcViscCoeffEl(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<SIMD<double>> elu_ipts,
                       FlatMatrix<SIMD<double>> res_ipts,
                       const double hi, double & coeff) const
  {
    cout << "no overload for CalcViscCoeffEl with FlatMatrix<SIMD>" << endl;
  }

  void ApplyM1 (int tentnr, double tstar, FlatMatrixFixWidth<COMP> u,
                FlatMatrixFixWidth<COMP> res, LocalHeap & lh);

  void ApplyM (int tentnr, double tstar,
               FlatMatrixFixWidth<COMP> u, FlatMatrixFixWidth<COMP> res,
               LocalHeap & lh);

  void Tent2Cyl (int tentnr, FlatMatrixFixWidth<COMP> u,
                 FlatMatrixFixWidth<COMP> uhat,
                 double tstar, LocalHeap & lh);

  template <typename T = SIMD<double>>
  void TransformBackIR(const SIMD_BaseMappedIntegrationRule & mir,
                       FlatMatrix<T> grad, FlatMatrix<T> u) const
  {
    throw Exception ("TransformBack for FlatMatrix<SIMD> not available");
  }

  void PropagateSARK(int stages, int substeps,
		     BaseVector & hu, BaseVector & hu_init,
		     LocalHeap & lh);

};


#endif // CONSERVATIONLAW_HPP
