#include "tents.hpp"
#include <python_ngstd.hpp>


// python export of tent mesh
auto ExportTimeSlab(py::module &m)
{
  auto pyname = "TentPitchedSlab";
  auto pydocu = "Tent pitched slab in  space + 1 time dimensions";
  auto pyslab = py::class_<TentPitchedSlab, shared_ptr<TentPitchedSlab>>
    (m, pyname, pydocu)
    .def(py::init([](shared_ptr<MeshAccess> ma, string method_name, int heapsize)
    {
      ngstents::PitchingMethod method = [method_name]{
	if(method_name == "edge") return ngstents::EEdgeGrad;
	else if(method_name == "vol") return ngstents::EVolGrad;
	else//just for static analyzers. the code should not reach this case
	  {
	    cout << "Invalid method! Setting edge algorithm as default..." << endl;
	    return ngstents::EEdgeGrad;
	  }
      }();
      auto tps = TentPitchedSlab(ma,  heapsize);
      tps.SetPitchingMethod(method);
      return tps;
    }),
	 py::arg("mesh"), py::arg("method"), py::arg("heapsize") = 1000000
	 )
    ;

  pyslab
    .def_readonly("mesh", &TentPitchedSlab::ma)
    .def("SetWavespeed", static_cast<void (TentPitchedSlab::*)(const double)>(&TentPitchedSlab::SetWavespeed))
    .def("PitchTents",[](shared_ptr<TentPitchedSlab> self,
			 const double dt, const bool local_ct, const double global_ct)
	 {
	   int dim = self->ma->GetDimension();
	   bool success = false;
	   switch(dim){
	   case 1:
	     success = self->PitchTents<1>(dt,local_ct,global_ct);
	     break;
	   case 2:
	     success = self->PitchTents<2>(dt,local_ct,global_ct);
	     break;
	   case 3:
	     success = self->PitchTents<3>(dt,local_ct,global_ct);
	     break;
	   default:
	     throw Exception("TentPitchedSlab not avaiable for dimension "+ToString(dim));
	   }
	   return success;
	 },
	 py::arg("dt"), py::arg("local_ct"), py::arg("global_ct")= 1.0)
    .def("GetNTents", &TentPitchedSlab::GetNTents)
    .def("GetNLayers", &TentPitchedSlab::GetNLayers)
    .def("GetSlabHeight", &TentPitchedSlab::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab::MaxSlope)
    .def("GetTent", &TentPitchedSlab::GetTent, pybind11::return_value_policy::reference_internal)    ;
  return pyslab;
}


void ExportTents(py::module & m) {

  py::class_<Tent, shared_ptr<Tent>>(m, "Tent", "Tent structure")
    .def_readonly("vertex", &Tent::vertex)
    .def_readonly("ttop", &Tent::ttop)
    .def_readonly("tbot", &Tent::tbot)
    .def_readonly("nbv", &Tent::nbv)
    .def_readonly("nbtime", &Tent::nbtime)
    .def_readonly("els", &Tent::els)
    .def_readonly("internal_facets", &Tent::internal_facets)
    .def("MaxSlope", &Tent::MaxSlope);

  ExportTimeSlab(m);
  // .def("DrawPitchedTentsPlt",[](shared_ptr<TentPitchedSlab> self)
  //      {
  // 	 if(self->ma->GetDimension() != 1)
  // 	   throw Exception("Only supported for 1D spatial mesh");
  // 	 else
  // 	   {
  // 	     py::list ret;
  // 	     for(int i = 0; i < self->GetNTents(); i++)
  // 	       {
  // 		 const Tent & tent = self->GetTent(i);
  // 		 py::list reti;
  // 		 reti.append(py::make_tuple(tent.vertex, tent.ttop,
  // 					    tent.tbot, tent.level));
  // 		 for(int j = 0; j< tent.nbv.Size(); j++)
  // 		   reti.append(py::make_tuple(tent.nbv[j],tent.nbtime[j]));
  // 		 ret.append(reti);
  // 	       }
  // 	     return ret;
  // 	   }
  //      })
  //    ;
  //ExportTimeSlab<1>(m)
  //  .def("DrawPitchedTentsPlt",[](shared_ptr<TentPitchedSlab<1>> self)
  //	 {
  //	   py::list ret;
  //	   for(int i = 0; i < self->GetNTents(); i++)
  //	     {
  //	       const Tent & tent = self->GetTent(i);
  //	       py::list reti;
  //	       reti.append(py::make_tuple(tent.vertex, tent.ttop,
  //					  tent.tbot, tent.level));
  //	       for(int j = 0; j< tent.nbv.Size(); j++)
  //		 reti.append(py::make_tuple(tent.nbv[j],tent.nbtime[j]));
  //	       ret.append(reti);
  //	     }
  //	   return ret;
  //	 })
  //  ;
  //
  //ExportTimeSlab<2>(m)
  //  .def("DrawPitchedTentsVTK",
  //	 [](shared_ptr<TentPitchedSlab<2>> self, string vtkfilename)
  //	 {
  //	   self->DrawPitchedTentsVTK(vtkfilename);
  //	 }, py::arg("vtkfilename")="output")
  //  .def("DrawPitchedTentsGL",
  //	 [](shared_ptr<TentPitchedSlab<2>> self)
  //	 {
  //	   int nlevels;
  //	   Array<int> tentdata;
  //	   Array<double> tenttimes;
  //	   self->DrawPitchedTentsGL(tentdata, tenttimes, nlevels);
  //	   py::list data, times;
  //	   for(auto i : Range(tentdata))
  //	     {
  //	       data.append(tentdata[i]);
  //	       // note: time values make sense only in 2D case.
  //	       // They are not used in 3D case, i.e. they are
  //	       // ignored by tents_visualization (ngsgui) and webgui.
  //	       times.append(tenttimes[i]);
  //	     }
  //	   return py::make_tuple(data,times,self->GetNTents(),nlevels);
  //	 })
  //  ;
  //
  //ExportTimeSlab<3>(m)
  //  .def("DrawPitchedTentsGL",
  //	 [](shared_ptr<TentPitchedSlab<3>> self)
  //	 {
  //	   int nlevels;
  //	   Array<int> tentdata;
  //	   Array<double> tenttimes;
  //	   self->DrawPitchedTentsGL(tentdata, tenttimes, nlevels);
  //	   py::list data;
  //	   for(auto i : Range(tentdata))
  //	     {
  //	       data.append(tentdata[i]);
  //	     }
  //	   return py::make_tuple(data, self->GetNTents(), nlevels);
  //	 })
  //  ;
}


PYBIND11_MODULE(tents, m) {
  ExportTents(m);
}
