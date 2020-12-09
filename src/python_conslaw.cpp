#include "conservationlaw.hpp"
#include <python_ngstd.hpp>

shared_ptr<ConservationLaw> CreateBurgers(const shared_ptr<TentPitchedSlab> & tps, const int & order);

shared_ptr<ConservationLaw> CreateWave(const shared_ptr<TentPitchedSlab> & tps, const int & order);

typedef ConservationLaw CL;

shared_ptr<CL> CreateConsLaw(const shared_ptr<TentPitchedSlab> & tps,
			     const string & eqn, const int & order)
{
  int dim = tps->ma->GetDimension();
  shared_ptr<CL> cl = nullptr;
  if (eqn=="burgers")
    cl = CreateBurgers(tps,order);
  else if(eqn=="wave")
    cl = CreateWave(tps,order);
    throw Exception(string("unknown equation '"+eqn+"'"));
  return cl;
}


void ExportConsLaw(py::module & m)
{
  py::class_<CL, shared_ptr<CL>>
    (m,"ConservationLaw", "Conservation Law")
    .def(py::init([](const shared_ptr<TentPitchedSlab> & tps,
  		     const string & eqn, const int & order)
  		  -> shared_ptr<CL>
                  {
                    auto cl = CreateConsLaw(tps, eqn, order);
                    cl->SetBC(); //use old style bc numbers for now
                    return cl;
                  }),
         py::arg("tentslab"),
  	 py::arg("equation"),
         py::arg("order"))
    .def("__str__",[](shared_ptr<CL> self){ return self->equation; })
    .def_property_readonly("space", [](shared_ptr<CL> self)
  			   -> shared_ptr<FESpace>
                           {
  			     return self->fes;
  			   })
    .def_property_readonly("sol", [](shared_ptr<CL> self)
  			   {
  			     return self->gfu;
  			   })
    .def_property_readonly("res", [](shared_ptr<CL> self)
  			   {
  			     return self->gfres;
  			   })
    .def_property_readonly("nu", [](shared_ptr<CL> self)
  			   {
  			     return self->gfnu;
  			   })
    // Set the initial data
    .def("SetInitial",
         [](shared_ptr<CL> self, shared_ptr<CoefficientFunction> cf)
         {
           SetValues(cf,*(self->gfu),VOL,0,*(self->pylh));
           self->uinit = self->u; // set data used for b.c.
         })
    .def("PropagateSARK",
         [](shared_ptr<CL> self, shared_ptr<BaseVector> vecu, int stages, int substeps)
         {
           self->PropagateSARK(stages, substeps,*vecu,*(self->uinit),*(self->pylh));
         }, py::arg("vec"), py::arg("stages") = 2, py::arg("substeps") = 1)
    ;
}

// Just for initial proof of concept (we need to move from DG -> MTP)
extern void ExportSymbolicDG (py::module & m); 

PYBIND11_MODULE(conslaw, m) {
  ExportConsLaw(m);
  ExportSymbolicDG(m);
}
