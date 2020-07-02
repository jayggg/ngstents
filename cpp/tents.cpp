#include "tents.hpp"
#include <python_ngstd.hpp>


//////////////// For handling periodicity //////////////////////////////////

// Get the slave vertex elements for a master vertex in periodic case 3D
void GetVertexElements(shared_ptr<MeshAccess> ma, int vnr_master,
                       const FlatArray<int> vnr_slaves, Array<int> & elems)
{
  ma->GetVertexElements(vnr_master,elems);
  for(auto slave : vnr_slaves)
    for(auto elnr : ma->GetVertexElements(slave))
      elems.Append(elnr);
}

Array<int> Tent::vmap;


/////////////////// Tent meshing ///////////////////////////////////////////

template <int DIM> void
TentPitchedSlab<DIM>::PitchTents(double dt, double wavespeed, LocalHeap & lh) 
{				 
  auto cf = make_shared<ConstantCoefficientFunction>(wavespeed);
  PitchTents(dt, cf, lh);
}


template <int DIM> void
TentPitchedSlab <DIM>::PitchTents(double dt,
				  shared_ptr<CoefficientFunction> wavespeed,
				  LocalHeap & lh)
{

  // maps regular vertices to themselves, periodic slave vertices to masters
  auto & vmap = Tent::vmap;
  vmap.SetSize(ma->GetNV());

  for (int i : Range(ma->GetNV()))
    vmap[i] = i;

  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_nodes = ma->GetPeriodicNodes(NT_VERTEX, idnr);
      for (const auto& per_verts : periodic_nodes)
        vmap[per_verts[1]] = vmap[per_verts[0]];
    }

  // element-wise maximal wave-speeds
  Array<double> cmax (ma->GetNE());

  // compute edge-based max time-differences
  Array<double> edge_refdt(ma->GetNEdges());

  edge_refdt = 1e99;  // max-double ??
  // check if edge is contained in mesh
  BitArray fine_edges(ma->GetNEdges());
  fine_edges.Clear();

  // Compute a reference dt for each edge based on edge length and
  // the wavespeed on each element
  for (Ngs_Element el : ma->Elements(VOL))
    {
      ElementId ei = ElementId(el);
      ELEMENT_TYPE eltype = ma->GetElType(ei);
      IntegrationRule ir (eltype, 0);
      ElementTransformation & trafo = ma->GetTrafo (ei, lh);
      MappedIntegrationPoint<DIM,DIM> mip(ir[0], trafo);
      cmax[el.Nr()] = wavespeed->Evaluate(mip);

      for (int e : el.Edges())
	{
          auto pnts = ma->GetEdgePNums(e);
          auto v1 = pnts[0], v2 = pnts[1];
	  double len = L2Norm (ma-> template GetPoint<DIM>(v1)
                             - ma-> template GetPoint<DIM>(v2));
	  edge_refdt[e] = min (edge_refdt[e], len/cmax[el.Nr()]);
          fine_edges.SetBit(e);
	}
    }
  // remove periodic edges
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_edges = ma->GetPeriodicNodes(NT_EDGE, idnr);
      for (const auto& per_edges : periodic_edges)
        fine_edges.Clear(per_edges[1]);
    }

  // Set the reference dt for each vertex to be the min of the reference dt
  // values for its adjacent edges.
  Array<double> vertex_refdt(ma->GetNV());
  vertex_refdt = 1e99;
  for (int e : IntRange (0, ma->GetNEdges()))
    if(fine_edges.Test(e))
      {
        auto vts = ma->GetEdgePNums (e);
        int v1 = vts[0], v2 = vts[1];
        vertex_refdt[v1] = min (vertex_refdt[v1], edge_refdt[e]);
        vertex_refdt[v2] = min (vertex_refdt[v2], edge_refdt[e]);
      }

  Array<double> tau(ma->GetNV());  // advancing front values at vertices
  tau = 0.0;
  // max time increase allowed at vertex, depends on tau of neighbors
  Array<double> ktilde(ma->GetNV());
  ktilde = vertex_refdt;

  // array of vertices ready for pitching a tent 
  Array<int> ready_vertices;
  Array<bool> vertex_ready(ma->GetNV());
  vertex_ready = false;
  for (int i = 0; i < ma->GetNV(); i++)
    if(vmap[i]==i) // non-periodic
      {
        ready_vertices.Append (i);
        vertex_ready[i] = true;
      }

  // build vertex2edge and vertex2vertex tables
  TableCreator<int> create_v2e, create_v2v;
  for ( ; !create_v2e.Done(); create_v2e++, create_v2v++)
    {
      for (int e : IntRange (0, ma->GetNEdges()))
        if(fine_edges.Test(e))
          {
            auto vts = ma->GetEdgePNums (e);
            int v1 = vts[0], v2 = vts[1];
            if(v1==vmap[v1]) // non-periodic
              {
                create_v2v.Add (v1, v2);
                create_v2e.Add (v1, e);
              }
            else
              {
                create_v2v.Add (vmap[v1], v2);
                create_v2e.Add (vmap[v1], e);
              }
            if(v2==vmap[v2])
              {
                create_v2v.Add (v2, v1);
                create_v2e.Add (v2, e);
              }
            else
              {
                create_v2v.Add (vmap[v2], v1);
                create_v2e.Add (vmap[v2], e);
              }
          }
    }


  Table<int> v2v = create_v2v.MoveTable();
  Table<int> v2e = create_v2e.MoveTable();

  // added for periodic tents
  TableCreator<int> create_slave_verts(ma->GetNV());
  for ( ; !create_slave_verts.Done(); create_slave_verts++)
    {
      for(auto i : Range(vmap))
        if(vmap[i]!=i)
          create_slave_verts.Add(vmap[i],i);
    }
  Table<int> slave_verts = create_slave_verts.MoveTable();

  Array<int> latest_tent(ma->GetNV()), vertices_level(ma->GetNV());
  latest_tent = -1;
  vertices_level = 0;

  // ---------------------------------------------
  // Main loop: constructs one tent each iteration
  // ---------------------------------------------
  while (ready_vertices.Size())
    {
      int minlevel = 1000;
      int posmin = 0;
      // Choose tent pole vertex vi and remove it from vertex_ready
      for(int i = 0; i < ready_vertices.Size(); i++)
	if(vertices_level[ready_vertices[i]] < minlevel)
	  {
	    minlevel = vertices_level[ready_vertices[i]];
	    posmin = i;
	  }
      int vi = ready_vertices[posmin];
      ready_vertices.DeleteElement(posmin);
      vertex_ready[vi] = false;

      // advance by ktilde:
      Tent * tent = new Tent;
      tent->vertex = vi;
      tent->tbot = tau[vi];
      tent->ttop = min (dt, tau[vi]+ktilde[vi]);
      tent->level = vertices_level[vi]; // 0;
      tau[vi] = tent->ttop;

      // update level of neighbor vertices
      for (int nb : v2v[vi])
	{
          nb = vmap[nb]; // only update master if periodic
	  tent->nbv.Append (nb);
	  tent->nbtime.Append (tau[nb]);
	  if(vertices_level[nb] < tent->level +1)
	    vertices_level[nb] = tent->level + 1;
          // tent number is just array index in tents
	  if (latest_tent[nb] != -1)
	      tents[latest_tent[nb]]->dependent_tents.Append (tents.Size());
	}
      latest_tent[vi] = tents.Size();
      vertices_level[vi]++;

      // Set tent internal facets
      if(DIM==1)
        // vertex itself represents the only internal edge/facet
       	tent->internal_facets.Append (vi);
      else if (DIM == 2)
        for (int e : v2e[vi]) tent->internal_facets.Append (e);
      else
        {
	  // DIM == 3 => internal facets are faces
          ArrayMem<int,4> fpnts;
          for (auto elnr : ma->GetVertexElements(vi))
            for (auto f : ma->GetElement(ElementId(VOL,elnr)).Faces())
              {
                ma->GetFacetPNums(f, fpnts);
                if (fpnts.Contains(vi) && !tent->internal_facets.Contains(f))
                  tent->internal_facets.Append(f);
              }
        }

      if(slave_verts[vi].Size()==0)
        ma->GetVertexElements (vi, tent->els);
      else
        GetVertexElements(ma,vi,slave_verts[vi],tent->els);

      // update max step ktilde for neighbors, and append them
      // to ready_vertices (if not there) once ktilde is large enough
      for (int nb : v2v[vi])
	{
          nb = vmap[nb]; // map periodic vertices
          if (tau[nb] >= dt) continue;
	  double kt = 1e99;
	  for (int nb2_index : v2v[nb].Range())
	    {
	      int nb2 = vmap[v2v[nb][nb2_index]];
	      double kt1 = tau[nb2]-tau[nb]+edge_refdt[v2e[nb][nb2_index]];
	      kt = min (kt, kt1);
	    }

	  ktilde[nb] = kt;
	  if (kt > 0.5 * vertex_refdt[nb])
            if (!vertex_ready[nb])
              {
                ready_vertices.Append (nb);
                vertex_ready[nb] = true;
              }
	}
      tents.Append (tent);
    }

  // set lists of internal facets of each element of each tent
  ParallelFor
    (Range(tents),
     [&] (int i) 
     {
       Array<int> dnums;
       Tent & tent = *tents[i];
       
       TableCreator<int> elfnums_creator(tent.els.Size());
       
       for ( ; !elfnums_creator.Done(); elfnums_creator++)  {
	 for(int j : Range(tent.els)) {
	   
	   auto fnums = ma->GetElFacets (tent.els[j]);
	   for(int fnum : fnums)
	     if (tent.internal_facets.Pos(fnum) != -1)
	       elfnums_creator.Add(j,fnum);
	 }
       }
       tent.elfnums = elfnums_creator.MoveTable();
     });

  // build dependency graph (used by RunParallelDependency)
  TableCreator<int> create_dag(tents.Size());
  for ( ; !create_dag.Done(); create_dag++)
    {
      for (int i : tents.Range())
	for (int d : tents[i]->dependent_tents)
	  create_dag.Add(i, d);
    }
  tent_dependency = create_dag.MoveTable();

  // set advancing front data members
  SetFrontData();

}


template <int DIM> void
TentPitchedSlab<DIM>::SetFrontData() {


  // The plan is to set these here:
  
  // gradphi_bot, gradphi_top;
  // Array<Vector<double>> delta; // phi_top - phi_bot
  // Array<Vector<>> graddelta;
  // Table<Matrix<>> gradphi_facet_bot, gradphi_facet_top;
  // Table<Vector<double>> delta_facet;

}


///////////////////// Output routines //////////////////////////////////////

void VTKOutputTents(shared_ptr<MeshAccess> maptr, Array<Tent*> & tents,
                    string filename)
{
  const MeshAccess & ma = *maptr;
  ofstream out(filename+".vtk");
  Array<Vec<3>> points;
  Array<INT<4>> cells;
  Array<int> level, tentnr;
  int ptcnt = 0;

  for(int i : Range(tents))
    {
      int firstpt = ptcnt;
      Tent & tent = *tents[i];
      Vec<2> pxy = ma.GetPoint<2> (tent.vertex);
      points.Append (Vec<3> (pxy(0), pxy(1), tent.tbot));
      points.Append (Vec<3> (pxy(0), pxy(1), tent.ttop));
      INT<4> tet(ptcnt,ptcnt+1,0,0);
      ptcnt+=2;

      for (int elnr : tent.els)
	{
	  Ngs_Element el = ma.GetElement(ElementId(VOL,elnr));

	  for (int v : el.Vertices())
	    if (v != tent.vertex)
	      {
		pxy = ma.GetPoint<2> (v);
		points.Append (Vec<3> (pxy(0), pxy(1),
                               tent.nbtime[tent.nbv.Pos(v)]));
	      }
	  for (int j = 2; j < 4; j++)
	    tet[j] = ptcnt++;

	  cells.Append(tet);
	}
      for(int j = firstpt; j < ptcnt; j++)
	{
	  level.Append(tent.level);
	  tentnr.Append(i);
	}
    }

  // header
  out << "# vtk DataFile Version 3.0" << endl;
  out << "vtk output" << endl;
  out << "ASCII" << endl;
  out << "DATASET UNSTRUCTURED_GRID" << endl;


  out << "POINTS " << points.Size() << " float" << endl;
  for (auto p : points)
    out << p << endl;


  out << "CELLS " << cells.Size() << " " << 5 * cells.Size() << endl;
  for (auto c : cells)
    out << 4 <<" " << c << endl;


  out << "CELL_TYPES " << cells.Size() << endl;
  for (auto c : cells)
    out << "10 " << endl;

  out << "CELL_DATA " << cells.Size() << endl;
  out << "POINT_DATA " << points.Size() << endl;

  out << "FIELD FieldData " << 2 << endl;

  out << "tentlevel" << " 1 " << level.Size() << " float" << endl;
  for (auto i : level)
    out << i << " ";
  out << endl;

  out << "tentnumber" << " 1 " << tentnr.Size() << " float" << endl;
  for (auto i : tentnr)
    out << i << " ";
  out << endl;
}



ostream & operator<< (ostream & ost, const Tent & tent) 
{
  ost << "vertex: " << tent.vertex << ", tbot = " << tent.tbot
      << ", ttop = " << tent.ttop << endl;
  ost << "neighbour vertices: " << endl;
  for (int k = 0; k < tent.nbv.Size(); k++)
    ost << k << ": " << tent.nbv[k] << " " << tent.nbtime[k] << endl;
  ost << "elements: " << endl << tent.els << endl;
  ost << "internal_facets: " << endl << tent.internal_facets << endl;
  ost << "elfnums: " << endl << tent.elfnums << endl;
  return ost;
}


///////////// Instantiate in expected dimensions ///////////////////////////

template class TentPitchedSlab<1>;
template class TentPitchedSlab<2>;
template class TentPitchedSlab<3>;



///////////// For python export ////////////////////////////////////////////

void ExportTents(py::module & m) {

  py::class_<TentPitchedSlab<2>, shared_ptr<TentPitchedSlab<2>>>
    (m, "TPS2", "Tent pitched slab in 2 space + 1 time dimensions")
    .def(py::init([](shared_ptr<MeshAccess> ma, double dt, double c, int heapsize)
		  {
		    auto tlh = LocalHeap(heapsize, "Tents heap");
		    auto tps = TentPitchedSlab<2>(ma);
		    tps.PitchTents(dt, c, tlh);
		    return tps;
		  }), 
      py::arg("mesh"), py::arg("dt"), py::arg("c"),
      py::arg("heapsize") = 1000000
      )

    .def("GetNTents",&TentPitchedSlab<2>::GetNTents);
    
}
