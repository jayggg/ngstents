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
shared_ptr<ConservationLaw> CreateSymbolicConsLaw (const shared_ptr<GridFunction> & gfu,
						   const shared_ptr<TentPitchedSlab> & tps,
						   const shared_ptr<CoefficientFunction> & flux,
						   const shared_ptr<CoefficientFunction> & numflux,
						   const shared_ptr<CoefficientFunction> & invmap,
						   const shared_ptr<CoefficientFunction> & cf_reflect);

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
  		     const string & eqn,
		     optional<Region> outflow, optional<Region> inflow,
		     optional<Region> reflect, optional<Region> transparent)
  		  -> shared_ptr<CL>
                  {
                    auto cl = CreateConsLaw(gfu, tps, eqn);
		    // set boundary data
                    if(outflow.has_value())
                      cl->SetBC(0,outflow.value().Mask());
                    if(reflect.has_value())
                      cl->SetBC(1,reflect.value().Mask());
                    if(inflow.has_value())
                      cl->SetBC(2,inflow.value().Mask());
                    if(transparent.has_value())
                      {
                        if (eqn=="wave")
                          cl->SetBC(3,transparent.value().Mask());
                        else
                          throw Exception ("Transparent boundary just available for wave equation!");
                      }
                    cl->CheckBC(); //use old style bc numbers if no regions set
                    return cl;
                  }),
         py::arg("gridfunction"), py::arg("tentslab"), py::arg("equation"),
	 py::arg("outflow")=nullptr, py::arg("inflow")=nullptr,
         py::arg("reflect")=nullptr, py::arg("transparent")=nullptr)
    .def(py::init([](const shared_ptr<GridFunction> & gfu,
     		     const shared_ptr<TentPitchedSlab> & tps,
     		     py::object Flux,
		     py::object NumFlux,
		     py::object InverseMap,
		     optional<py::object> ReflectBnd,
		     const bool compile,
		     optional<Region> outflow,
		     optional<Region> inflow,
		     optional<Region> reflect)
     		  -> shared_ptr<CL>
		  {
     		    py::object u = py::cast (gfu->GetFESpace()).attr("TrialFunction")();
     		    py::object flux_u = Flux( u );
     		    py::object numflux_u = NumFlux( u, u.attr("Other")() );
		    py::object invmap = InverseMap( u );

     		    shared_ptr<CoefficientFunction> cpp_flux_u =
     		      py::extract<shared_ptr<CoefficientFunction>> (flux_u)();
     		    cpp_flux_u = Compile(cpp_flux_u, compile);

		    shared_ptr<CoefficientFunction> cpp_numflux_u =
		      py::extract<shared_ptr<CoefficientFunction>> (numflux_u)();
		    cpp_numflux_u = Compile(cpp_numflux_u, compile);

		    shared_ptr<CoefficientFunction> cpp_invmap =
		      py::extract<shared_ptr<CoefficientFunction>> (invmap)();
		    cpp_invmap = Compile(cpp_invmap, compile);

		    // CF for reflecting boundary condition
		    shared_ptr<CoefficientFunction> cpp_cf_reflect = nullptr;
		    if(ReflectBnd.has_value())
		      {
			py::object cf_reflect = ReflectBnd.value()( u );
			cpp_cf_reflect = py::extract<shared_ptr<CoefficientFunction>> (cf_reflect)();
			cpp_cf_reflect = Compile(cpp_cf_reflect, compile);
		      }

		    auto cl = CreateSymbolicConsLaw(gfu, tps, cpp_flux_u, cpp_numflux_u, cpp_invmap,
						    cpp_cf_reflect);
		    // set boundary data
                    if(outflow.has_value())
                      cl->SetBC(0,outflow.value().Mask());
                    if(reflect.has_value())
                      cl->SetBC(1,reflect.value().Mask());
                    if(inflow.has_value())
                      cl->SetBC(2,inflow.value().Mask());
		    cl->CheckBC(); //use old style bc numbers if no regions set
		    return cl;
		  }),
	 py::arg("gridfunction"),
	 py::arg("tentslab"),
     	 py::arg("flux"),
	 py::arg("numflux"),
	 py::arg("inversemap"),
	 py::arg("reflectbnd")=nullptr,
	 py::arg("compile")=false,
	 py::arg("outflow")=nullptr,
	 py::arg("inflow")=nullptr,
         py::arg("reflect")=nullptr
	 )
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
    .def_property_readonly("tau", [](shared_ptr<CL> self)
			   {
			     return self->cftau;
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
    .def("SetBoundaryCF",[](shared_ptr<CL> self, Region region, shared_ptr<CoefficientFunction> cf)
         {
	   auto maxbcnr = self->GetMaxBCNr();
	   self->SetBC(maxbcnr, region);
           self->SetBoundaryCF(maxbcnr, cf);
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
    ;
}

// Just for initial proof of concept (we need to move from DG -> MTP)
extern void ExportSymbolicDG (py::module & m); 

PYBIND11_MODULE(_pyconslaw, m) {
  m.attr("__name__") = "ngstents.conslaw";
  m.attr("__package__") = "ngstents";
  ExportConsLaw(m);
  ExportSymbolicDG(m);
}
