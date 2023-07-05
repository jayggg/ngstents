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

typedef CoefficientFunction CF;
shared_ptr<ConservationLaw>
CreateSymbolicConsLaw(const shared_ptr<GridFunction> & gfu,
		      const shared_ptr<TentPitchedSlab> & tps,
		      const shared_ptr<ProxyFunction> & proxy_u,
		      const shared_ptr<ProxyFunction> & proxy_uother,
		      const shared_ptr<CF> & flux,
		      const shared_ptr<CF> & numflux,
		      const shared_ptr<CF> & invmap,
		      const shared_ptr<CF> & entropy,
		      const shared_ptr<CF> & entropyflux,
		      const shared_ptr<CF> & numentropyflux,
		      const bool compile);

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
		     const bool compile,
		     optional<py::object> Entropy,
		     optional<py::object> EntropyFlux,
		     optional<py::object> NumEntropyFlux,
		     optional<py::object> ViscosityCoefficient)
     		  -> shared_ptr<CL>
		  {
		    // proxies for u and u.Other()
     		    py::object u = py::cast (gfu->GetFESpace()).attr("TrialFunction")();
		    py::object uother = u.attr("Other")();

		    // CF for flux
		    py::object flux_u = Flux( u );
     		    shared_ptr<CF> cpp_flux_u =
     		      py::extract<shared_ptr<CF>> (flux_u)();
		    cpp_flux_u = Compile(cpp_flux_u, compile, 0);

		    //  CF for numerical flux
		    py::object numflux_u = NumFlux( u, uother );
		    shared_ptr<CF> cpp_numflux_u =
		      py::extract<shared_ptr<CF>> (numflux_u)();
		    cpp_numflux_u = Compile(cpp_numflux_u, compile, 0);

		    // CF for inverse map
		    py::object invmap = InverseMap( u );
		    shared_ptr<CF> cpp_invmap =
		      py::extract<shared_ptr<CF>> (invmap)();
		    cpp_invmap = Compile(cpp_invmap, compile, 0);

		    // CF for entropy residual
		    shared_ptr<CF> cpp_entropy = nullptr;
		    shared_ptr<CF> cpp_entropyflux = nullptr;
		    shared_ptr<CF> cpp_numentropyflux = nullptr;
		    if(Entropy.has_value())
		      {
			py::object cf_entropy = Entropy.value()( u );
			cpp_entropy = py::extract<shared_ptr<CF>> (cf_entropy)();
			cpp_entropy = Compile(cpp_entropy, compile, 0);
		      }
		    if(EntropyFlux.has_value())
		      {
			py::object cf_entropyflux = EntropyFlux.value()( u );
			cpp_entropyflux =
			  py::extract<shared_ptr<CF>> (cf_entropyflux)();
			cpp_entropyflux = Compile(cpp_entropyflux, compile, 0);
		      }
		    if(NumEntropyFlux.has_value())
		      {
			py::object cf_numentropyflux = NumEntropyFlux.value()( u, uother);
			cpp_numentropyflux =
			  py::extract<shared_ptr<CF>> (cf_numentropyflux)();
			cpp_numentropyflux = Compile(cpp_numentropyflux, compile, 0);
		      }

		    bool entropy_functions = false;
		    if (cpp_entropy && cpp_entropyflux && cpp_numentropyflux)
		      entropy_functions = true;
		    else if (cpp_entropy || cpp_entropyflux || cpp_numentropyflux)
		      throw Exception("at least one entropy function missing");

		    auto proxy_u = py::extract<shared_ptr<ProxyFunction>> (u)();
		    auto proxy_uother = py::extract<shared_ptr<ProxyFunction>> (uother)();

		    auto cl = CreateSymbolicConsLaw(gfu, tps, proxy_u, proxy_uother,
						    cpp_flux_u, cpp_numflux_u, cpp_invmap,
						    cpp_entropy, cpp_entropyflux,
						    cpp_numentropyflux, compile);
		    if(ViscosityCoefficient.has_value())
		      {
			py::object cf_visccoeff =
			  ViscosityCoefficient.value()( u, py::cast(cl->proxy_res) );
			auto cpp_visccoeff =
			  py::extract<shared_ptr<CF>> (cf_visccoeff)();
			cpp_visccoeff = Compile(cpp_visccoeff, compile, 0);
			cl->SetViscosityCoefficient(cpp_visccoeff);
		      }
		    else if (entropy_functions)
		      throw Exception("function for ViscosityCoefficient missing");

		    // use old style bc numbers
		    cl->CheckBC();

		    return cl;
		  }),
	 py::arg("gridfunction"),
	 py::arg("tentslab"),
     	 py::arg("flux"),
	 py::arg("numflux"),
	 py::arg("inversemap"),
	 py::arg("compile")=false,
	 py::arg("entropy")=nullptr,
	 py::arg("entropyflux")=nullptr,
         py::arg("numentropyflux")=nullptr,
	 py::arg("visccoeff")=nullptr
	 )
    .def_property_readonly("tentslab", [](shared_ptr<CL> self)
                           {
  			     return self->tps;
  			   }, "returns the tent pitched time slab")
    .def("__str__",[](shared_ptr<CL> self){ return self->equation; })
    .def_property_readonly("space", [](shared_ptr<CL> self)
  			   -> shared_ptr<FESpace>
                           {
  			     return self->fes;
  			   }, "returns finite element space")
    .def_property_readonly("sol", [](shared_ptr<CL> self)
  			   {
  			     return self->gfu;
  			   }, "returns grid function representing the solution")
    .def_property_readonly("u_minus", [](shared_ptr<CL> self)
                           {
  			     return self->proxy_u;
  			   }, "returns trial function u(x-s*n) for s->0^+ and the normal vector n")
    .def_property_readonly("time", [](shared_ptr<CL> self)
			   {
			     return self->cftau;
			   }, "the time coordinate of a spacetime point in an advancing front")
    .def_property_readonly("res", [](shared_ptr<CL> self)
  			   {
  			     return self->gfres;
  			   }, "entropy residual for nonlinear hyperbolic eq")
    .def_property_readonly("nu", [](shared_ptr<CL> self)
  			   {
  			     return self->gfnu;
  			   }, "entropy viscosity for nonlinear hyperbolic eq")
    // Set the initial data
    .def("SetInitial",
         [](shared_ptr<CL> self, shared_ptr<CF> cf)
         {
           SetValues(cf,*(self->gfu),VOL,0,*(self->pylh));
           self->uinit->Set(1.0,*(self->u)); // set data used for b.c.
         })
    // Set vector field for advection equation
    .def("SetVectorField",
         [](shared_ptr<CL> self, shared_ptr<CF> cf)
         {
           self->SetVectorField(cf);
         })
    .def("SetBoundaryCF",[](shared_ptr<CL> self, Region region, shared_ptr<CF> cf)
         {
	   // bcnr's 0 - 3 used for default boundary conditions
	   self->SetBC(4, region);
           self->SetBoundaryCF(4 , cf);
         })
    .def("SetBoundaryCF",[](shared_ptr<CL> self, shared_ptr<CF> cf)
         {
	   // bcnr's 0 - 3 used for default boundary conditions
	   self->SetBC(4, Region(self->ma,BND,".*"));
	   self->SetBoundaryCF(4, cf);
         })
    .def("SetNumEntropyFlux",[](shared_ptr<CL> self, shared_ptr<CF> cf)
         {
	   self->SetNumEntropyFlux(cf);
         })
    .def("SetMaterialParameters",
         [](shared_ptr<CL> self,
	    shared_ptr<CF> cf_mu,
	    shared_ptr<CF> cf_eps)
	 {
	   self->SetMaterialParameters(cf_mu,cf_eps);
	 }, py::arg("mu"), py::arg("eps"))
    .def("SetTentSolver",
         [](shared_ptr<CL> self, string method, int stages, int substeps)
         {
           self->SetTentSolver(method, stages, substeps);
         },
	 py::arg("method")="SARK",
	 py::arg("stages")=2, py::arg("substeps")=1,
	 R"(
         Parameters:--
           method: SARK (Structure-Aware Runge Kutta), or 
                   SAT  (Structure-Aware Taylor).
           stages: determines the order of time stepper.
           substeps: number of subtents each tent should be divided into
             before applying the tent solver method.
           ----------- )"
	 )
    .def("SetIdx3d",
         [](shared_ptr<CL> self, py::list lst)
         {
           auto idx3d = Array<shared_ptr<Table<int>>>();
           std::vector<std::map<int, int>> vec
             = lst.cast<std::vector<std::map<int, int>>>();
           for (auto vmp : vec) {
              TableCreator<int> create_vmap;
              for ( ; !create_vmap.Done(); create_vmap++)
                {
                  for ( const auto &p : vmp )
                      create_vmap.Add (p.first, p.second);
                }
              auto vmap = create_vmap.MoveTable();
              auto tptr = make_shared<Table<int>>(vmap);
              idx3d.Append(tptr);
           }
           self->vis3d = make_shared<Visualization3D>(idx3d);
         }, "Set index for visualization on a 3D mesh", py::arg("idx3d"))
    .def("Propagate",
         [](shared_ptr<CL> self,
            shared_ptr<GridFunction> hdgf)
         {
            self->Propagate(*(self->pylh), hdgf);
         }, "GridFunction vector for visualization on 3D mesh"
         , py::arg("hdgf")=nullptr)
    ;
}


PYBIND11_MODULE(_pyconslaw, m) {
  m.attr("__name__") = "ngstents.conslaw";
  m.attr("__package__") = "ngstents";
  ExportConsLaw(m);
}
