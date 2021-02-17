#include "conservationlaw.hpp"
#include <python_ngstd.hpp>

shared_ptr<ConservationLaw> CreateBurgers(const shared_ptr<GridFunction> & gfu,
					  const shared_ptr<TentPitchedSlab> & tps);
shared_ptr<ConservationLaw> CreateEuler(const shared_ptr<GridFunction> & gfu,
					const shared_ptr<TentPitchedSlab> & tps);
shared_ptr<ConservationLaw> CreateWave(const shared_ptr<GridFunction> & gfu,
				       const shared_ptr<TentPitchedSlab> & tps);
shared_ptr<ConservationLaw> CreateAdvection(const shared_ptr<GridFunction> & gfu,
					    const shared_ptr<TentPitchedSlab> & tps);
shared_ptr<ConservationLaw> CreateMaxwell(const shared_ptr<GridFunction> & gfu,
					  const shared_ptr<TentPitchedSlab> & tps);

typedef ConservationLaw CL;

shared_ptr<CL> CreateConsLaw(const shared_ptr<GridFunction> & gfu,
			     const shared_ptr<TentPitchedSlab> & tps,
			     const string & eqn)
{
  shared_ptr<CL> cl = nullptr;
  if (eqn=="burgers")
    cl = CreateBurgers(gfu, tps);
  else if(eqn=="euler")
    cl = CreateEuler(gfu, tps);
  else if(eqn=="wave")
    cl = CreateWave(gfu, tps);
  else if(eqn=="advection")
    cl = CreateAdvection(gfu, tps);
  else if(eqn=="maxwell")
    cl = CreateMaxwell(gfu, tps);
  else
    throw Exception(string("unknown equation '"+eqn+"'"));
  return cl;
}


void ExportConsLaw(py::module & m)
{
  py::class_<CL, shared_ptr<CL>>
    (m,"ConservationLaw", "Conservation Law")
    .def(py::init([](const shared_ptr<GridFunction> & gfu,
		     const shared_ptr<TentPitchedSlab> & tps,
  		     const string & eqn)
  		  -> shared_ptr<CL>
                  {
                    auto cl = CreateConsLaw(gfu, tps, eqn);
                    cl->SetBC(); //use old style bc numbers for now
                    return cl;
                  }),
         py::arg("gridfunction"),
  	 py::arg("tentslab"),
         py::arg("equation"))
    .def_property_readonly("tentslab", [](shared_ptr<CL> self)
                           {
  			     return self->tps;
  			   })
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
           self->uinit->Set(1.0,*(self->u)); // set data used for b.c.
         })
    // Set flux field
    .def("SetVectorField",
         [](shared_ptr<CL> self, shared_ptr<CoefficientFunction> cf)
         {
           self->SetVectorField(cf);
         })
    .def("SetMaterialParameters",
         [](shared_ptr<CL> self,
	    shared_ptr<CoefficientFunction> cf_mu,
	    shared_ptr<CoefficientFunction> cf_eps)
	 {
	   self->SetMaterialParameters(cf_mu,cf_eps);
	 }, py::arg("mu"), py::arg("eps"))
    .def("SetTentSolver",
         [](shared_ptr<CL> self, string method, int stages, int substeps)
         {
           self->SetTentSolver(method, stages, substeps);
         }, py::arg("method") = "SAT", py::arg("stages") = 2, py::arg("substeps") = 1)
    .def("Propagate",
         [](shared_ptr<CL> self)
         {
           self->Propagate(*(self->pylh));
         })
    .def("PropagateSAT",
         [](shared_ptr<CL> self, int stages, int substeps)
         {
           self->PropagateSAT(stages, substeps, self->gfu->GetVector(), *(self->uinit), *(self->pylh));
         }, py::arg("stages") = 2, py::arg("substeps") = 1)
    .def("PropagateSARK",
         [](shared_ptr<CL> self, int stages, int substeps)
         {
           self->PropagateSARK(stages, substeps, self->gfu->GetVector(), *(self->uinit), *(self->pylh));
         }, py::arg("stages") = 2, py::arg("substeps") = 1)
    ;
}

// Just for initial proof of concept (we need to move from DG -> MTP)
extern void ExportSymbolicDG (py::module & m); 

PYBIND11_MODULE(conslaw, m) {
  ExportConsLaw(m);
  ExportSymbolicDG(m);
}
