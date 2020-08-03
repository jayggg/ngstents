#include "tents.hpp"
#include <python_ngstd.hpp>
#include <limits>


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
TentPitchedSlab<DIM>::PitchTents(double dt, double wavespeed)
{
  auto cf = make_shared<ConstantCoefficientFunction>(wavespeed);
  PitchTents(dt, cf);
}


template <int DIM> void
TentPitchedSlab <DIM>::PitchTents(double dt,
				  shared_ptr<CoefficientFunction> wavespeed)
{
  this->dt = dt; // set it so that GetSlabHeight can return it

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

  edge_refdt = std::numeric_limits<double>::max();
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
  vertex_refdt =  std::numeric_limits<double>::max();
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
      for(size_t i = 0; i < ready_vertices.Size(); i++)
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
	  if(vertices_level[nb] < tent->level + 1)
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
	  double kt = std::numeric_limits<double>::max();
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
       Tent & tent = *tents[i];
       TableCreator<int> elfnums_creator(tent.els.Size());

       for ( ; !elfnums_creator.Done(); elfnums_creator++)  {
	 for(int j : Range(tent.els)) {

	   auto fnums = ma->GetElFacets (tent.els[j]);
	   for(int fnum : fnums)
	     if (tent.internal_facets.Pos(fnum) !=
                 tent.internal_facets.ILLEGAL_POSITION)
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

  // set advancing front gradients
  ParallelFor
    (Range(tents),
     [&] (int i)     {

       Tent & tent = *tents[i];
       int nels = tent.els.Size();
       tent.gradphi_bot.SetSize(nels);
       tent.gradphi_top.SetSize(nels);

       for (int j : Range(nels)) { //  loop over elements in a tent

	 ElementId ej (VOL, tent.els[j]);
	 ELEMENT_TYPE eltype = ma->GetElType(ej);
	 BaseScalarFiniteElement *fe;
	 switch(DIM)
	   {
	   case 1: fe = new (lh) ScalarFE<ET_SEGM,1>(); break;
	   case 2: fe = new (lh) ScalarFE<ET_TRIG,1>(); break;
	   default: fe = new (lh) ScalarFE<ET_TET,1>();
	   }
	 Vector<> shape_nodal(fe->GetNDof());
	 Matrix<> dshape_nodal(fe->GetNDof(), DIM);
	 Vector<> coef_bot, coef_top; // coefficient of tau (top & bot)
	 coef_bot.SetSize(fe->GetNDof());
	 coef_top.SetSize(fe->GetNDof());
	 auto vnums = ma->GetElVertices (ej);
	 for (size_t k = 0; k < vnums.Size(); k++) {
	   if (vnums[k] == tent.vertex)  { // central vertex
	     coef_bot(k) = tent.tbot;
	     coef_top(k) = tent.ttop;
	   }
	   else
	     for (size_t l = 0; l < tent.nbv.Size(); l++)
	       if (tent.nbv[l] == vnums[k])
		 coef_bot(k) = coef_top(k) = tent.nbtime[l];
	 }

	 IntegrationRule ir(eltype, 0);
	 ElementTransformation & trafo = ma->GetTrafo (ej, lh);
	 MappedIntegrationPoint<DIM, DIM> mip(ir[0], trafo);
	 tent.gradphi_bot[j].SetSize(DIM);
	 tent.gradphi_top[j].SetSize(DIM);
	 fe->CalcMappedDShape(mip, dshape_nodal);
	 tent.gradphi_bot[j] = Trans(dshape_nodal) * coef_bot;
	 tent.gradphi_top[j] = Trans(dshape_nodal) * coef_top;
       }
     });
}

template <int DIM> double TentPitchedSlab<DIM>::GetPoleHeight(const int vi, const Array<double> & tau, const Array<double> & cmax, LocalHeap & lh) const{

  constexpr ELEMENT_TYPE  el_type = [&]() -> ELEMENT_TYPE {
    if constexpr (DIM == 1) return ET_SEGM;
    else if constexpr (DIM == 2) return ET_TRIG;
    else if constexpr (DIM == 3) return ET_TET;
    else return ET_POINT;//
  }();//this assumes that there is only one type of element per mesh

  //number of vertices of the current element (always the simplex associated to DIM)
  constexpr int n_vertices = DIM+1;
  //finite element created for calculating the barycentric coordinates
  BaseScalarFiniteElement *my_fel = new (lh) ScalarFE<el_type,1>();
  // array of all elements containing vertex vi
  Array<int>
      els;
  ma->GetVertexElements(vi, els);

  double pole_height = std::numeric_limits<double>::max();
  //vector containing the advancing front time for each vertex (except vi)
  Vector<double> coeff_vec(n_vertices);
  //gradient of basis functions onthe current element
  Matrix<double> gradphi(n_vertices,DIM);
  for (int el : els)
    {
      ElementId ei(VOL,el);
      //c_max^2
      const double c_max_sq = cmax[ei.Nr()] * cmax[ei.Nr()]; 
      //mapping of the current el
      ElementTransformation &trafo = ma->GetTrafo(ei, lh);
      //vertices of current el
      auto v_indices = ma->GetElVertices(ei);
      //vi position in current el
      const auto local_vi = v_indices.Pos(vi);
      //integration rule for reference el
      IntegrationRule ir(el_type,1);
      //integration point on deformed element
      MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
      my_fel->CalcMappedDShape(mip,gradphi);

      //sets the coefficient vec
      for (auto k : IntRange(0, v_indices.Size()))
        coeff_vec[k] = tau[v_indices[k]];
      coeff_vec[local_vi] = 0;
      cout<<"el = "<<ei.Nr()<<"\tvi = "<<vi;
      cout<<"\tlocal vi = "<<local_vi<<endl;
      for (auto k : IntRange(0, v_indices.Size()))
        {
          cout << "\tk = "<<k<<"\tv_n = "<<v_indices[k];
          cout <<"\talpha_k = "<<coeff_vec[k]<<"\ttau_k = "<<tau[v_indices[k]]<<endl;
        }
      /*writing the quadratic eq for tau_v
       \tau_{v}^2 ( \nabla\,\phi_{vi} \cdot \nabla\, \phi_{vi})
                   ^^^^^^^^^^^^^^^^^^alpha^^^^^^^^^^^^^^^^^^^^^^
       +\tau_{v} (\sum_{i \neq v} (\tau_i \nabla,\phi_i \cdot \nabla\,\phi_v))
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^beta^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       + (\sum_{i\neq v}\sum_{j\neq v}\left(\tau_i\tau_j \nabla \phi_i 
``        \cdot \nabla \phi_j\right)-\frac{1}{c}^2)
         ^^^^^^^^^^^^^^^^^^gamma^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

      
      // square of the norm of gradphi_vi
      const double alpha = InnerProduct(gradphi.Row(local_vi),gradphi.Row(local_vi));
      Vec<DIM> tau_i_grad_i = Trans(gradphi) * coeff_vec;
      const double beta = 2 * InnerProduct(tau_i_grad_i,gradphi.Row(local_vi));
      const double gamma = InnerProduct(tau_i_grad_i, tau_i_grad_i) - 1.0/c_max_sq;
      const double sq_delta = sqrt(beta * beta - 4 * alpha * gamma);
      try
        {
          if(beta * beta - 4 * alpha * gamma < 0) throw std::domain_error("Error solving quadratic equation!");
        }
      catch (const std::domain_error &error)
        {
          Matrix<double> jac(DIM,DIM);
          trafo.CalcJacobian(ir[0],jac);
          cout << error.what() << endl;
          cout<<"transform =\n"<<jac<<endl;
          cout<<"gradphi =\n"<<gradphi<<endl;
          cout << "a = "<<alpha<<"\tb = "<<beta<<"\tc = "<<gamma<<endl;
          cout << "delta = " << beta * beta - 4 * alpha * gamma << endl;
          exit(-1);
        }

      //since alpha will always be positive (and so will sq_delta),
      //we dont need to consider the solution (-beta-sq_delta)/(2*alpha)
      const double sol1 = (-beta+sq_delta)/(2*alpha);
      const double sol2 = (-beta-sq_delta)/(2*alpha);
      
      if(sol1 > sol2) pole_height = min(pole_height,sol1);
      else pole_height = min(pole_height,sol2);
      cout<<"a = "<<alpha<<"\tb = "<<beta<<"\tc = "<<gamma<<endl;
      cout<<"sol 1:"<< sol1<<"\tsol 2:"<<sol2<<endl;
    }

  //the return value is actually the ADVANCE in the current vi
  pole_height -= tau[vi];
  try
    {
      if(pole_height < 0) throw std::runtime_error("Error in pole height!");
    }
  catch (const std::runtime_error &error)
    {
      cout<< error.what() << endl;
      cout<<"pole_height(old) = "<<pole_height+tau[vi]<<endl;
      cout<<"pole_height = "<<pole_height<<endl;
      cout<<"tau(old) = "<<tau[vi]<<endl;
      exit(-1);
    }
  return pole_height;
 }

template <int DIM> void
TentPitchedSlab<DIM>::PitchTentsGradient(double dt, double wavespeed)
{
  auto cf = make_shared<ConstantCoefficientFunction>(wavespeed);
  PitchTentsGradient(dt, cf);
}

template <int DIM> void
TentPitchedSlab <DIM>::PitchTentsGradient(double dt,
				  shared_ptr<CoefficientFunction> wavespeed)
{
  this->dt = dt; // set it so that GetSlabHeight can return it

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
  // check if edge is contained in mesh
  BitArray fine_edges(ma->GetNEdges());
  fine_edges.Clear();  
  
  //identify fine edges and evaluate maximum wavespeed for each element
  for (Ngs_Element el : ma->Elements(VOL))
    {
      //set all edges belonging to the mesh
      for (int e : el.Edges())	
        fine_edges.SetBit(e);

      auto ei = ElementId(el);
      auto eltype = ma->GetElType(ei);
      ElementTransformation & trafo = ma->GetTrafo (ei, lh);
      IntegrationRule ir(eltype, 0);
      MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
      cmax[el.Nr()] = wavespeed->Evaluate(mip);
    }
  // remove periodic edges
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_edges = ma->GetPeriodicNodes(NT_EDGE, idnr);
      for (const auto& per_edges : periodic_edges)
        fine_edges.Clear(per_edges[1]);
    }
  Array<double> tau(ma->GetNV());  // advancing front values at vertices
  tau = 0.0;

  /*The reference dt would previously be computed
    analysing the minimum length of the edges around
    a given vertex. However, this does not imply
    in causality. The inverse of the gradient of the
    barycentric coordinates (H1 vertex functions)
    are now used for this purpose, as one can obtain the
    minimum distance to the opposite facet in relation
    to the vertex patch. */

  // compute vertex-based max time-differences
  Array<double> vertex_refdt(ma->GetNV());
  vertex_refdt = std::numeric_limits<double>::max();

  // array of vertices ready for pitching a tent
  Array<int> ready_vertices;
  // array for checking if a given vertex is ready
  Array<bool> vertex_ready(ma->GetNV());
  vertex_ready = false;
  for (auto i = 0; i < ma->GetNV(); i++)
    if(vmap[i]==i) // non-periodic
      {
        ready_vertices.Append (i);
        vertex_ready[i] = true;
        vertex_refdt[i] = GetPoleHeight(i, tau, cmax, lh);
      }
  // max time increase allowed at vertex, depends on tau of neighbors
  Array<double> ktilde(ma->GetNV());
  //at the beginning the advancing front is at a constant t=0
  //so ktilde can be set as vertex_refdt
  ktilde = vertex_refdt;

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

  //array containing the latest tent in which the vertex was included
  Array<int> latest_tent(ma->GetNV());
  //array containing the current level of the tent to be pitched at vertex vi
  Array<int> vertices_level(ma->GetNV());
  latest_tent = -1;
  vertices_level = 0;

  // ---------------------------------------------
  // Main loop: constructs one tent each iteration
  // ---------------------------------------------
  while (ready_vertices.Size())
    {
      //minimum vertex level among the ready_vertices
      int minlevel = 1000;
      //position of tent pole vertex vi in ready_vertices
      int posmin = 0;
      // Choose tent pole vertex vi and remove it from vertex_ready
      for(auto i = 0; i < ready_vertices.Size(); i++)
	if(vertices_level[ready_vertices[i]] < minlevel)
	  {
	    minlevel = vertices_level[ready_vertices[i]];
	    posmin = i;
	  }
      const int vi = ready_vertices[posmin];
      ready_vertices.DeleteElement(posmin);
      vertex_ready[vi] = false;

      // advance by ktilde:
      Tent * tent = new Tent;
      tent->vertex = vi;
      tent->tbot = tau[vi];
      tent->ttop = min (dt, tau[vi]+ktilde[vi]);
      tent->level = vertices_level[vi]; // 0;
      tau[vi] = tent->ttop;

      //add neighboring vertices and update their level
      for (int nb : v2v[vi])
	{
          nb = vmap[nb]; // only update master if periodic
	  tent->nbv.Append (nb);
	  tent->nbtime.Append (tau[nb]);
          //update level of vertices if needed
	  if(vertices_level[nb] < tent->level + 1)
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
          //points contained in a given facet
          ArrayMem<int,4> fpnts;
          for (auto elnr : ma->GetVertexElements(vi))
            for (auto f : ma->GetElement(ElementId(VOL,elnr)).Faces())
              {
                //get facet vertices
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
	  double kt = std::numeric_limits<double>::max();
	  for (int nb2_index : v2v[nb].Range())
	    {
	      const int nb2 = vmap[v2v[nb][nb2_index]];
	      const double kttemp = GetPoleHeight(nb2,tau,cmax,lh);
	      kt = min (kt, kttemp);
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
       Tent & tent = *tents[i];
       TableCreator<int> elfnums_creator(tent.els.Size());

       for ( ; !elfnums_creator.Done(); elfnums_creator++)  {
	 for(int j : Range(tent.els)) {

	   auto fnums = ma->GetElFacets (tent.els[j]);
	   for(int fnum : fnums)
	     if (tent.internal_facets.Pos(fnum) !=
                 tent.internal_facets.ILLEGAL_POSITION)
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

  // set advancing front gradients
  ParallelFor
    (Range(tents),
     [&] (int i)     {

       Tent & tent = *tents[i];
       int nels = tent.els.Size();
       tent.gradphi_bot.SetSize(nels);
       tent.gradphi_top.SetSize(nels);

       for (int j : Range(nels)) { //  loop over elements in a tent

	 ElementId ej (VOL, tent.els[j]);
	 ELEMENT_TYPE eltype = ma->GetElType(ej);
	 BaseScalarFiniteElement *fe;
	 switch(DIM)
	   {
	   case 1: fe = new (lh) ScalarFE<ET_SEGM,1>(); break;
	   case 2: fe = new (lh) ScalarFE<ET_TRIG,1>(); break;
	   default: fe = new (lh) ScalarFE<ET_TET,1>();
	   }
	 Vector<> shape_nodal(fe->GetNDof());
	 Matrix<> dshape_nodal(fe->GetNDof(), DIM);
	 Vector<> coef_bot, coef_top; // coefficient of tau (top & bot)
	 coef_bot.SetSize(fe->GetNDof());
	 coef_top.SetSize(fe->GetNDof());
	 auto vnums = ma->GetElVertices (ej);
	 for (size_t k = 0; k < vnums.Size(); k++) {
	   if (vnums[k] == tent.vertex)  { // central vertex
	     coef_bot(k) = tent.tbot;
	     coef_top(k) = tent.ttop;
	   }
	   else
	     for (size_t l = 0; l < tent.nbv.Size(); l++)
	       if (tent.nbv[l] == vnums[k])
		 coef_bot(k) = coef_top(k) = tent.nbtime[l];
	 }

	 IntegrationRule ir(eltype, 0);
	 ElementTransformation & trafo = ma->GetTrafo (ej, lh);
	 MappedIntegrationPoint<DIM, DIM> mip(ir[0], trafo);
	 tent.gradphi_bot[j].SetSize(DIM);
	 tent.gradphi_top[j].SetSize(DIM);
	 fe->CalcMappedDShape(mip, dshape_nodal);
	 tent.gradphi_bot[j] = Trans(dshape_nodal) * coef_bot;
	 tent.gradphi_top[j] = Trans(dshape_nodal) * coef_top;
       }
     });
}


template <int DIM>
double TentPitchedSlab <DIM>::MaxSlope() {

  // Return  max(|| gradphi_top||, ||gradphi_bot||)

  double maxgrad = 0.0;
  ParallelFor
    (Range(tents),
     [&] (int i)
     {
       Tent & tent = *tents[i];
       for (int j : Range(tent.els.Size())) {
	 auto norm = L2Norm(tent.gradphi_top[j]);
	 AtomicMax(maxgrad, norm);
       }
     });
  return maxgrad;
}



///////////////////// Output routines //////////////////////////////////////

template <int DIM> void
TentPitchedSlab <DIM>::DrawPitchedTentsVTK(string filename)
{
  ofstream out(filename+".vtk");
  Array<Vec<3>> points;
  Array<INT<4>> cells;
  Array<int> level, tentnr;
  int ptcnt = 0;

  for(int i : Range(GetNTents()))
    {
      int firstpt = ptcnt;
      const Tent & tent = GetTent(i);
      Vec<2> pxy = ma->GetPoint<2> (tent.vertex);
      points.Append (Vec<3> (pxy(0), pxy(1), tent.tbot));
      points.Append (Vec<3> (pxy(0), pxy(1), tent.ttop));
      INT<4> tet(ptcnt,ptcnt+1,0,0);
      ptcnt+=2;

      for (int elnr : tent.els)
	{
	  Ngs_Element el = ma->GetElement(ElementId(VOL,elnr));

	  for (int v : el.Vertices())
	    if (v != tent.vertex)
	      {
		pxy = ma->GetPoint<2> (v);
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


// Used with OpenGL (ngsgui/tents_visualization) and WebGL (webgui).
template <int DIM> void
TentPitchedSlab <DIM>::DrawPitchedTentsGL(
    Array<int> & tentdata, Array<double> & tenttimes, int & nlevels)
{
  nlevels = 0;
  tentdata.SetAllocSize(4*tents.Size());
  tenttimes.SetAllocSize(4*tents.Size());

  for(int i : Range(GetNTents()))
    {
      const Tent & tent = GetTent(i);
      for(int el : Range(tent.els))
        {
          tentdata.Append(i);
          tentdata.Append(tent.level);
          tentdata.Append(tent.vertex);
          tentdata.Append(tent.els[el]);
          if(tent.level > nlevels)
            nlevels = tent.level;

          if constexpr (DIM == 2)
          {
            auto verts = ma->GetElVertices(ElementId(VOL,tent.els[el]));
            for(auto v : verts)
              {
                auto pos = tent.nbv.Pos(v);
                if (pos != tent.nbv.ILLEGAL_POSITION)
                  tenttimes.Append(tent.nbtime[pos]);
                else
                  tenttimes.Append(tent.tbot);
              }
            tenttimes.Append(tent.ttop);
          }
        }
    }
  nlevels+=1;
}


ostream & operator<< (ostream & ost, const Tent & tent)
{
  ost << "vertex: " << tent.vertex << ", tbot = " << tent.tbot
      << ", ttop = " << tent.ttop << endl;
  ost << "neighbour vertices: " << endl;
  for (size_t k = 0; k < tent.nbv.Size(); k++)
    ost << k << ": " << tent.nbv[k] << " " << tent.nbtime[k] << endl;
  ost << "elements: " << endl << tent.els << endl;
  ost << "internal_facets: " << endl << tent.internal_facets << endl;
  ost << "elfnums: " << endl << tent.elfnums << endl;
  return ost;
}


///////////// TentDataFE ///////////////////////////////////////////////////


TentDataFE::TentDataFE(const Tent & tent, const FESpace & fes,
		       const MeshAccess & ma, LocalHeap & lh)
  : fei(tent.els.Size(), lh),
    iri(tent.els.Size(), lh),
    miri(tent.els.Size(), lh),
    trafoi(tent.els.Size(), lh),
    mesh_size(tent.els.Size(), lh),
    agradphi_bot(tent.els.Size(), lh),
    agradphi_top(tent.els.Size(), lh),
    adelta(tent.els.Size(), lh),
    felpos(tent.internal_facets.Size(), lh),
    firi(tent.internal_facets.Size(), lh),
    mfiri1(tent.internal_facets.Size(), lh),
    mfiri2(tent.internal_facets.Size(), lh),
    agradphi_botf1(tent.internal_facets.Size(), lh),
    agradphi_topf1(tent.internal_facets.Size(), lh),
    agradphi_botf2(tent.internal_facets.Size(), lh),
    agradphi_topf2(tent.internal_facets.Size(), lh),
    anormals(tent.internal_facets.Size(), lh),
    adelta_facet(tent.internal_facets.Size(), lh)
{
  int dim = ma.GetDimension();

  FlatArray<BaseScalarFiniteElement*> fe_nodal(tent.els.Size(),lh);
  FlatArray<FlatVector<double>> coef_delta(tent.els.Size(),lh);
  FlatArray<FlatVector<double>> coef_top(tent.els.Size(),lh);
  FlatArray<FlatVector<double>> coef_bot(tent.els.Size(),lh);

  for (size_t i = 0; i < tent.els.Size(); i++)
    {
      ElementId ei(VOL, tent.els[i]);
      // ranges and dofs were previously members of tent
      Array<int> dnums;
      fes.GetDofNrs (ei, dnums);
      ranges.Append(IntRange(dnums.Size()) + dofs.Size());
      dofs += dnums;

      fei[i] = &fes.GetFE (ei, lh);
      iri[i] = new (lh) SIMD_IntegrationRule(fei[i]->ElementType(),
					     2*fei[i]->Order());
      trafoi[i] = &ma.GetTrafo (ei, lh);
      miri[i] =  &(*trafoi[i]) (*iri[i], lh);

      mesh_size[i] = pow(fabs((*miri[i])[0].GetJacobiDet()[0]),
                         1.0/miri[i]->DimElement());

      int nipt = miri[i]->Size();
      agradphi_bot[i].AssignMemory(dim, nipt, lh);
      agradphi_top[i].AssignMemory(dim, nipt, lh);
      adelta[i].AssignMemory(nipt, lh);

      switch(dim)
        {
        case 1: fe_nodal[i] = new (lh) ScalarFE<ET_SEGM,1>(); break;
        case 2: fe_nodal[i] = new (lh) ScalarFE<ET_TRIG,1>(); break;
        default: fe_nodal[i] = new (lh) ScalarFE<ET_TET,1>();
        }

      coef_top[i].AssignMemory(fe_nodal[i]->GetNDof(), lh);
      coef_bot[i].AssignMemory(fe_nodal[i]->GetNDof(), lh);
      auto vnums = ma.GetElVertices(ei);
      for (size_t k = 0; k < vnums.Size(); k++)
        {
          auto mapped_vnum = tent.vmap[vnums[k]]; // map periodic vertices
          auto pos = tent.nbv.Pos(mapped_vnum);
          if (pos != tent.nbv.ILLEGAL_POSITION)

            coef_bot[i](k) = coef_top[i](k) = tent.nbtime[pos];

          else {

	    coef_bot[i](k) = tent.tbot;
	    coef_top[i](k) = tent.ttop;
	  }
        }
      fe_nodal[i]->EvaluateGrad(*miri[i], coef_top[i], agradphi_top[i]);
      fe_nodal[i]->EvaluateGrad(*miri[i], coef_bot[i], agradphi_bot[i]);

      coef_delta[i].AssignMemory(fe_nodal[i]->GetNDof(), lh);
      coef_delta[i] = coef_top[i]-coef_bot[i];
      fe_nodal[i]->Evaluate(*iri[i], coef_delta[i], adelta[i]);
    }
    nd = dofs.Size();

  for (size_t i = 0; i < tent.internal_facets.Size(); i++)
    {
      int order = 0;
      INT<2> loc_facetnr;

      ArrayMem<int,2> elnums;
      ArrayMem<int,2> elnums_per;
      ma.GetFacetElements(tent.internal_facets[i],elnums);

      bool periodic_facet = false;
      int facet2;
      if(elnums.Size() < 2)
        {
          facet2 = ma.GetPeriodicFacet(tent.internal_facets[i]);
          if(facet2 != tent.internal_facets[i])
            {
              ma.GetFacetElements (facet2, elnums_per);
              if (elnums_per.Size())
                {
                  periodic_facet = true;
                  elnums.Append(elnums_per[0]);
                }
            }
        }
      SIMD_IntegrationRule * simd_ir_facet;

      // Note: size_t(-1) = 18446744073709551615 = ar.ILLEGAL_POSITION
      // Is it OK to rely on this cast/conversion?  Is there a better way
      // to check for the condition that an element is not in the array?
      felpos[i] = INT<2,size_t>(size_t(-1));
      for(int j : Range(elnums.Size()))
        {
          felpos[i][j] = tent.els.Pos(elnums[j]);
          if(felpos[i][j] != size_t(-1))
            {
              int elorder = fei[felpos[i][j]]->Order();
              if(elorder > order) order = elorder;

              auto fnums = ma.GetElFacets (elnums[j]);
              int fnr = tent.internal_facets[i];
              if(periodic_facet)
                {
                  auto pos = fnums.Pos(facet2);
                  if(pos != size_t(-1))
                    fnr = facet2; // change facet nr to slave
                }
              for (int k : Range(fnums.Size()))
                if (fnums[k] == fnr) loc_facetnr[j] = k;

              auto & trafo = *trafoi[felpos[i][j]];

              auto vnums = ma.GetElVertices (elnums[j]);
              Facet2ElementTrafo transform(trafo.GetElementType(), vnums);

              auto etfacet = ElementTopology::GetFacetType (
                  trafo.GetElementType(), loc_facetnr[j]);
              if(j == 0)
                {
                  simd_ir_facet = new (lh)
		    SIMD_IntegrationRule (etfacet, 2*order+1);

                  // quick fix to avoid usage of TP elements (slows down)
                  simd_ir_facet->SetIRX(nullptr);
                }

              firi[i][j] = &transform(loc_facetnr[j],*simd_ir_facet, lh);
              if(j == 0)
                {
                  mfiri1[i] = &trafo(*firi[i][j], lh);
                  mfiri1[i]->ComputeNormalsAndMeasure(trafo.GetElementType(),
                                                      loc_facetnr[j]);

                  anormals[i].AssignMemory(dim,mfiri1[i]->Size(),lh);
                  adelta_facet[i].AssignMemory(mfiri1[i]->Size(),lh);
                  agradphi_botf1[i].AssignMemory(dim, mfiri1[i]->Size(), lh);
                  agradphi_topf1[i].AssignMemory(dim, mfiri1[i]->Size(), lh);

                  for (size_t k : Range(mfiri1[i]->Size()))
                    {
                      switch(dim)
                        {
                        case 1:
                          anormals[i].Col(k) =
                            static_cast<const SIMD<DimMappedIntegrationPoint<1>>&>
                              ((*mfiri1[i])[k]).GetNV(); break;
                        case 2:
                          anormals[i].Col(k) =
                            static_cast<const SIMD<DimMappedIntegrationPoint<2>>&>
                              ((*mfiri1[i])[k]).GetNV(); break;
                        default:
                          anormals[i].Col(k) =
                            static_cast<const SIMD<DimMappedIntegrationPoint<3>>&>
                              ((*mfiri1[i])[k]).GetNV();
                        }
                    }
                  size_t elpos = felpos[i][j];
                  fe_nodal[elpos]->Evaluate(*firi[i][j],coef_delta[elpos],
                                            adelta_facet[i]);
                  fe_nodal[elpos]->EvaluateGrad(*mfiri1[i],coef_bot[elpos],
                                            agradphi_botf1[i]);
                  fe_nodal[elpos]->EvaluateGrad(*mfiri1[i],coef_top[elpos],
                                            agradphi_topf1[i]);
                }
              else
                {
                  mfiri2[i] = &trafo(*firi[i][j], lh);
                  mfiri2[i]->ComputeNormalsAndMeasure(trafo.GetElementType(),
                                                      loc_facetnr[j]);
                  agradphi_botf2[i].AssignMemory(dim, mfiri2[i]->Size(), lh);
                  agradphi_topf2[i].AssignMemory(dim, mfiri2[i]->Size(), lh);
                  size_t elpos = felpos[i][j];
                  fe_nodal[elpos]->EvaluateGrad(*mfiri2[i],coef_bot[elpos],
                                                agradphi_botf2[i]);
                  fe_nodal[elpos]->EvaluateGrad(*mfiri2[i],coef_top[elpos],
                                                agradphi_topf2[i]);
                }
            }
        }
    }
}


///////////// Instantiate in expected dimensions ///////////////////////////

template class TentPitchedSlab<1>;
template class TentPitchedSlab<2>;
template class TentPitchedSlab<3>;


///////////// For python export ////////////////////////////////////////////

void ExportTents(py::module & m) {

  py::class_<Tent, shared_ptr<Tent>>(m, "Tent", "Tent structure")
    .def_readonly("vertex", &Tent::vertex)
    .def_readonly("ttop", &Tent::ttop)
    .def_readonly("tbot", &Tent::tbot)
    .def_readonly("nbv", &Tent::nbv)
    .def_readonly("nbtime", &Tent::nbtime)
    .def_readonly("els", &Tent::els)
    .def_readonly("internal_facets", &Tent::internal_facets);

  //
  // 1D spatial mesh
  //
  py::class_<TentPitchedSlab<1>, shared_ptr<TentPitchedSlab<1>>>
    (m, "TentPitchedSlab1", "Tent pitched slab in 1 space + 1 time dimensions")
    .def(py::init([](shared_ptr<MeshAccess> ma, double dt, double c,
                     bool exact, int heapsize)
		  {
		    auto tps = TentPitchedSlab<1>(ma, heapsize);
                    if(exact) tps.PitchTentsGradient(dt,c);
		    else tps.PitchTents(dt, c);
		    return tps;
		  }),
      py::arg("mesh"), py::arg("dt"), py::arg("c"),
      py::arg("exact") = false, py::arg("heapsize") = 1000000
      )

    .def_readonly("mesh", &TentPitchedSlab<1>::ma)
    .def("GetNTents", &TentPitchedSlab<1>::GetNTents)
    .def("GetSlabHeight", &TentPitchedSlab<1>::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab<1>::MaxSlope)
    .def("GetTent", &TentPitchedSlab<1>::GetTent)
    .def("DrawPitchedTentsPlt",[](shared_ptr<TentPitchedSlab<1>> self)
     {
       py::list ret;
       for(int i = 0; i < self->GetNTents(); i++)
         {
           const Tent & tent = self->GetTent(i);
           py::list reti;
           reti.append(py::make_tuple(tent.vertex, tent.ttop,
                                      tent.tbot, tent.level));
           for(int j = 0; j< tent.nbv.Size(); j++)
             reti.append(py::make_tuple(tent.nbv[j],tent.nbtime[j]));
           ret.append(reti);
         }
       return ret;
     })


     ; // please leave me on my own line

  py::class_<TentPitchedSlab<2>, shared_ptr<TentPitchedSlab<2>>>
    (m, "TentPitchedSlab2", "Tent pitched slab in 2 space + 1 time dimensions")
    .def(py::init([](shared_ptr<MeshAccess> ma, double dt, double c, bool exact, int heapsize)
		  {
		    auto tps = TentPitchedSlab<2>(ma, heapsize);
		    if(exact) tps.PitchTentsGradient(dt,c);
		    else tps.PitchTents(dt, c);
		    return tps;
		  }),
      py::arg("mesh"), py::arg("dt"), py::arg("c"),
      py::arg("exact") = false, py::arg("heapsize") = 1000000
      )

    .def_readonly("mesh", &TentPitchedSlab<2>::ma)
    .def("GetNTents", &TentPitchedSlab<2>::GetNTents)
    .def("GetSlabHeight", &TentPitchedSlab<2>::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab<2>::MaxSlope)
    .def("GetTent", &TentPitchedSlab<2>::GetTent)
    .def("DrawPitchedTentsVTK",
         [](shared_ptr<TentPitchedSlab<2>> self, string vtkfilename)
         {
           self->DrawPitchedTentsVTK(vtkfilename);
         }, py::arg("vtkfilename")="output")
    .def("DrawPitchedTentsGL",
         [](shared_ptr<TentPitchedSlab<2>> self)
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
               // They are not used in 3D case, i.e. they are
               // ignored by tents_visualization (ngsgui) and webgui.
               times.append(tenttimes[i]);
             }
           return py::make_tuple(data,times,self->GetNTents(),nlevels);
         })

     ; // please leave me on my own line

  py::class_<TentPitchedSlab<3>, shared_ptr<TentPitchedSlab<3>>>
    (m, "TentPitchedSlab3", "Tent pitched slab in 3 space + 1 time dimensions")
    .def(py::init([](shared_ptr<MeshAccess> ma, double dt, double c, bool exact,  int heapsize)
		  {
		    auto tps = TentPitchedSlab<3>(ma, heapsize);
		    if(exact) tps.PitchTentsGradient(dt,c);
		    else tps.PitchTents(dt, c);
		    return tps;
		  }),
      py::arg("mesh"), py::arg("dt"), py::arg("c"),
      py::arg("exact") = false, py::arg("heapsize") = 1000000
      )

    .def_readonly("mesh", &TentPitchedSlab<3>::ma)
    .def("GetNTents", &TentPitchedSlab<3>::GetNTents)
    .def("GetSlabHeight", &TentPitchedSlab<3>::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab<3>::MaxSlope)
    .def("GetTent", &TentPitchedSlab<3>::GetTent)
    .def("DrawPitchedTentsGL",
         [](shared_ptr<TentPitchedSlab<3>> self)
         {
           int nlevels;
           Array<int> tentdata;
           Array<double> tenttimes;
           self->DrawPitchedTentsGL(tentdata, tenttimes, nlevels);
           py::list data;
           for(auto i : Range(tentdata))
             {
               data.append(tentdata[i]);
             }
           return py::make_tuple(data, self->GetNTents(), nlevels);
         })

     ; // please leave me on my own line

}
