#include "conservationlaw.hpp"
#include <python_ngstd.hpp>

/////////////////// export to python ////////////////////////////////////////

shared_ptr<ConservationLaw> CreateBurgers(shared_ptr<MeshAccess> ma,
                                          int order, const Flags & flags);
typedef CoefficientFunction CF;
typedef GridFunction GF;
typedef ConservationLaw CL;

shared_ptr<CL> CreateConsLaw(shared_ptr<MeshAccess> ma, const string & eqn,
                             int order, const Flags & flags)
{
  shared_ptr<CL> cl = nullptr;
  if (eqn=="burgers")
    cl = CreateBurgers(ma,order,flags);
  else
    throw Exception(string("unknown equation '"+eqn+"'"));
  cl->equation = eqn;
  return cl;
}

void ExportConsLaw(py::module & m)
{
  py::class_<ConservationLaw, shared_ptr<CL>>
    (m,"ConsLaw", docu_string(R"raw_string(ConsLaw help)raw_string"))
    .def(py::init([](shared_ptr<MeshAccess> ma, const string & eqn,
                     int order, py::dict pyflags)
                  {
                    const Flags flags = py::extract<Flags> (pyflags)();
                    auto cl = CreateConsLaw(ma,eqn,order,flags);
                    cl->CheckBC(); //use old style bc numbers for now
                    return cl;
                  }),
         py::arg("mesh"),
         py::arg("eqn"),
         py::arg("order"),
         py::arg("flags")=py::dict())
    .def_property_readonly("equation", [] (shared_ptr<CL> self)
         { return self->Equation(); })
    .def_property_readonly("sol", [](shared_ptr<CL> self) { return self->gfu; })
    .def_property_readonly("res", [](shared_ptr<CL> self) { return self->gfres; })
    .def_property_readonly("nu", [](shared_ptr<CL> self) { return self->gfnu; })
    .def_property_readonly("space",
                           [](shared_ptr<CL> self) -> shared_ptr<FESpace>
                           { return self->fes; })

    // Set the initial data
    .def("SetInitial",
         [](shared_ptr<CL> self, shared_ptr<CF> cf)
         {
           SetValues(cf,*(self->gfu),VOL,0,*(self->pylh));
           self->uinit = self->u; // set data used for b.c.
         })
    .def("PropagatePicard",
         [](shared_ptr<CL> self, shared_ptr<BaseVector> vecu, int steps)
         {
           if(steps==-1)
             steps = 1;
           self->PropagatePicard(steps,*vecu,self->uinit,*(self->pylh));
         }, py::arg("vec"),py::arg("steps")=-1)

    /////////////////// tent pitching functions ////////////////////

    .def("PitchTents",
         [](shared_ptr<CL> self, double dt, py::object wavespeed)
         {
           if (py::extract<double> (wavespeed).check())
             {
               auto cf_wavespeed =
                 make_shared<ConstantCoefficientFunction>(
                     py::extract<double> (wavespeed)());
               self->PitchTents(dt,cf_wavespeed);
             }
           else if (py::extract<shared_ptr<CF>>(wavespeed).check())
             self->PitchTents(dt,py::extract<shared_ptr<CF>>(
                   wavespeed)());
           else
             throw Exception("wrong argument type for wavespeed in CreateTents");
         })
    .def("MaxSlope",[](shared_ptr<CL> self) { return self->MaxSlope(); })
      
    /////////////////// visualization functions ////////////////////
    
    .def("DrawPitchedTentsVTK",
         [](shared_ptr<CL> self, string vtkfilename)
         { 
           self->DrawPitchedTentsVTK(vtkfilename); 
         }, py::arg("vtkfilename")="output.vtk")
    .def("DrawPitchedTentsGL",
         [](shared_ptr<CL> self)
         {
           int nlevels;
           Array<int> tentdata;
           Array<double> tenttimes;
           self->DrawPitchedTentsGL(tentdata, tenttimes, nlevels);
           py::list data, times;
           for(auto i : Range(tentdata))
             {
               data.append(tentdata[i]);
               // note: time values make sense only in 2D case.
               // times are not used in 3D case, i.e. they are
               // ignored by tents_visualization (ngsgui) and webgui.
               times.append(tenttimes[i]);
             }
           return py::make_tuple(data,times,self->GetNTents(),nlevels);
         })
    .def("Tau",[](shared_ptr<CL> self) { return self->gftau; })

    ; // please keep me on my own line
}
