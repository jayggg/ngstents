#ifndef TENTSOLVER_HPP
#define TENTSOLVER_HPP

template <typename EQUATION, int DIM, int COMP, int ECOMP, bool SYMBOLIC> class T_ConservationLaw;

class TentSolver
{
public:
  TentSolver ()
  { };

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
  SARK (const shared_ptr<TCONSLAW> & atcl, int astages, int asubsteps)
    : tcl{atcl}, stages{astages}, substeps{asubsteps}
  {
    cout << "set up structure-aware Runge-Kutta time stepping with "
      +ToString(substeps)+" substeps within each tent" << endl;

    shared_ptr<L2HighOrderFESpace> fes_check = dynamic_pointer_cast<L2HighOrderFESpace>(atcl->fes);
    if(!fes_check)
      throw Exception("Structure-aware Runge-Kutta time stepping available for L2 spaces only");
  };
  
  void PropagateTent(const Tent & tent, BaseVector & hu,
		     const BaseVector & hu0, LocalHeap & lh) override;
};
  
#endif //TENTSOLVER_HPP
