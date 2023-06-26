#ifndef TENTSOLVER_HPP
#define TENTSOLVER_HPP

template <typename EQUATION, int DIM, int COMP, int ECOMP, bool SYMBOLIC> class T_ConservationLaw;

class TentSolver
{
public:
  TentSolver() = default;

  virtual void Setup() { };

  virtual void PropagateTent(const Tent & tent, BaseVector & hu,
			     const BaseVector & hu0, LocalHeap & lh) = 0;
};

template <typename TCONSLAW>
class SAT : public TentSolver
{
protected:
  const int stages;
  const int substeps;

  // pointer to T_ConservationLaw
  shared_ptr<TCONSLAW> tcl;
  static constexpr int COMP = TCONSLAW::NCOMP;
  
public:
  SAT (const shared_ptr<TCONSLAW> & atcl, int astages, int asubsteps)
    : tcl{atcl}, stages{astages}, substeps{asubsteps}
  {
    cout << "set up SAT timestepping with "+
      ToString(stages)+" stages and "+ToString(substeps)+" substeps/tent" << endl;

    shared_ptr<L2HighOrderFESpace> fes_check = dynamic_pointer_cast<L2HighOrderFESpace>(atcl->fes);
    if(!fes_check)
      throw Exception("Structure-aware Taylor time stepping available for L2 spaces only");
  };

  void Setup() override
  {
    tcl->DeriveBoundaryCF(stages);
  }

  void PropagateTent(const Tent & tent, BaseVector & hu,
		     const BaseVector & hu0, LocalHeap & lh) override;
};

template <typename TCONSLAW>
class SARK : public TentSolver
{
protected:
  const int stages;
  const int substeps;

  // pointer to T_ConservationLaw
  shared_ptr<TCONSLAW> tcl;
  static constexpr int COMP = TCONSLAW::NCOMP;
  static constexpr int ECOMP = TCONSLAW::NECOMP;
  
public:

  Matrix<> acoeff;
  Matrix<> dcoeff;
  Vector<> bcoeff;
  Vector<> ccoeff;

  SARK (const shared_ptr<TCONSLAW> & atcl, int astages, int asubsteps)
    : tcl{atcl}, stages{astages}, substeps{asubsteps}
  {
    shared_ptr<L2HighOrderFESpace> fes_check = dynamic_pointer_cast<L2HighOrderFESpace>(atcl->fes);
    if(!fes_check)
      throw Exception("Structure-aware Runge-Kutta time stepping available for L2 spaces only");

    cout << "set up " + ToString(stages) + "-stage ";

    // SARK coefficients below are derived from order conditions in
    // https://doi.org/10.1007/s42985-020-00020-4
    
    switch(stages)
      {
      case 1:
	acoeff = { {0.0} };
	dcoeff = { {0.0} };
	bcoeff = { 1.0 };
	ccoeff = { 0.0 };
	cout << "(first order) ";
	break;
      case 2:
	acoeff = { {0.0, 0.0},
		   {0.5, 0.0} };
	dcoeff = { {0.0, 0.0},
		   {0.5, 0.0} };
	bcoeff = { 0.0, 1.0 };
	ccoeff = { 0.0, 0.5 };
	cout << "(second order) ";
	break;
      case 3:
	acoeff = { {0.0, 0.0, 0.0},
		   {0.5, 0.0, 0.0},
		   {-1.0, 2.0, 0.0} };
	dcoeff = { {0.0, 0.0, 0.0},
		   {0.5, 0.0, 0.0},
		   {-3.0, 4.0, 0.0} };
	bcoeff = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };
	ccoeff = { 0.0, 0.5, 1.0 };
	cout << "(third order) ";
	break;
      case 5:
	
	// Fourth order 5-stage SARK method can be found in the dissertation
	// "Mapped Tent Pitching Schemes for Hyperbolic Systems"
	// by C. Wintersteiger (2020)

	acoeff = { {0.0, 0.0, 0.0, 0.0, 0.0},
		   {1.0/3.0, 0.0, 0.0, 0.0, 0.0},
		   {-1.0/3.0, 1.0, 0.0, 0.0, 0.0},
		   {-2.0/5.0, 7.0/5.0, 0.0, 0.0, 0.0},
		   {3.0/8.0, -1.0/8.0, 1.0/4.0, 0.0, 0.0} };
	dcoeff = { {0.0, 0.0, 0.0, 0.0, 0.0},
		   {1.0/3.0, 0.0, 0.0, 0.0, 0.0},
		   {-4.0/3.0, 2.0, 0.0, 0.0, 0.0},
		   {-9.0/5.0, 14.0/5.0, 0.0, 0.0, 0.0},
		   {11.0/16.0, -13.0/16.0, 5.0/16.0, 5.0/16.0, 0.0} };
	bcoeff = { 5.0/32.0, 3.0/32.0, 3.0/32.0, 5.0/32.0, 1.0/2.0 };
	ccoeff = { 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0/2.0 };
	cout << "(fouth order) ";
	break;
      default:
	throw Exception("no "+ToString(stages)+"-stage SARK method implemented");
      }
    cout << "SARK timestepping with "
      + ToString(substeps) + " substeps/tent" << endl;
  };

  void PropagateTent(const Tent & tent, BaseVector & hu,
		     const BaseVector & hu0, LocalHeap & lh) override;
};
  
#endif //TENTSOLVER_HPP
