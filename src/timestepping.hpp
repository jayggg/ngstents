#ifndef TIMESTEPPING_HPP
#define TIMESTEPPING_HPP

template <typename EQUATION, int DIM, int COMP, int ECOMP> class T_ConservationLaw;

class TimeStepping
{
public:
  TimeStepping ()
  { };

  // const Tent & tent instead of tentnr ?
  virtual void PropagateTent(const int tentnr, BaseVector & hu,
			     const BaseVector & hu0, LocalHeap & lh) = 0;
};

template <typename EQUATION, int DIM, int COMP, int ECOMP>
class SAT : public TimeStepping
{
protected:
  const int stages;
  const int substeps;
  shared_ptr<T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>> tcl;
  
public:
  SAT (const shared_ptr<T_ConservationLaw<EQUATION, DIM, COMP, ECOMP>> & atcl,
       int astages, int asubsteps)
    : tcl{atcl}, stages{astages}, substeps{asubsteps}
  {
    shared_ptr<L2HighOrderFESpace> fes_check = dynamic_pointer_cast<L2HighOrderFESpace>(atcl->fes);
    if(!fes_check)
      throw Exception("Structure-aware Taylor time stepping available for L2 spaces only");
  };
  
  void PropagateTent(const int tentnr, BaseVector & hu,
		     const BaseVector & hu0, LocalHeap & lh);
};
  

#endif //TIMESTEPPING_HPP
