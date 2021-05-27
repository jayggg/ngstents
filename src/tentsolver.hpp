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
    cout << "set up structure-aware Taylor time stepping with "+
      ToString(stages)+" stages and "+ToString(substeps)+" substeps within each tent" << endl;

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
    cout << "set up structure-aware Runge-Kutta time stepping with "
      +ToString(substeps)+" substeps within each tent" << endl;

    shared_ptr<L2HighOrderFESpace> fes_check = dynamic_pointer_cast<L2HighOrderFESpace>(atcl->fes);
    if(!fes_check)
      throw Exception("Structure-aware Runge-Kutta time stepping available for L2 spaces only");

    switch(stages)
      {
	// case 1:
	// break;
	// case 2:
	// break;
      case 3:
	acoeff = { {0.0, 0.0, 0.0},
		   {0.5, 0.0, 0.0},
		   {-1.0, 2.0, 0.0} };
	dcoeff = { {0.0, 0.0, 0.0},
		   {0.5, 0.0, 0.0},
		   {-3.0, 4.0, 0.0} };
	bcoeff = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };
	ccoeff = { 0.0, 0.5, 1.0 };
	break;
	// case 4:
	// break;
	// case 5:
	// break;
      default:
	throw Exception("no "+ToString(stages)+"-stage SARK method implemented");
      }
    cout << "a = " << endl << acoeff << endl;
    cout << "d = " << endl << dcoeff << endl;
    cout << "b = " << endl << bcoeff << endl;
    cout << "c = " << endl << ccoeff << endl;
  };

  void PropagateTent(const Tent & tent, BaseVector & hu,
		     const BaseVector & hu0, LocalHeap & lh) override;
};
  
#endif //TENTSOLVER_HPP
