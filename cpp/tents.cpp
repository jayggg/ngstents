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

void MapPeriodicVertices(shared_ptr<MeshAccess> ma)
{
  auto &vmap = Tent::vmap;
  vmap.SetSize(ma->GetNV());
  for (int i : Range(ma->GetNV()))
    vmap[i] = i;
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_nodes = ma->GetPeriodicNodes(NT_VERTEX, idnr);
      for (const auto& per_verts : periodic_nodes)
        vmap[per_verts[1]] = vmap[per_verts[0]];
    }
}


void RemovePeriodicEdges(shared_ptr<MeshAccess> ma, BitArray &fine_edges)
{
  for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
    {
      const auto & periodic_edges = ma->GetPeriodicNodes(NT_EDGE, idnr);
      for (const auto& per_edges : periodic_edges)
        fine_edges.Clear(per_edges[1]);
    }
}

Array<int> Tent::vmap;

/////////////////// Tent methods ////////////////////////////////////////////
double Tent::MaxSlope() const{
  double maxgrad = 0.0;
  for (int j : Range(this->els.Size()))
    {
      const auto norm = L2Norm(this->gradphi_top[j]);
      maxgrad =  max(maxgrad, norm);
    }
  return maxgrad;
}

/////////////////// Tent meshing ///////////////////////////////////////////

constexpr ELEMENT_TYPE EL_TYPE(int DIM)
{
  return DIM == 1 ? ET_SEGM : DIM == 2 ? ET_TRIG : ET_TET;
}//this assumes that there is only one type of element per mesh

template <int DIM>
bool TentPitchedSlab <DIM>::PitchTents(double dt, bool calc_local_ct, const double global_ct)
{
  if(has_been_pitched)
    {
      tents.DeleteAll();
    }
  if(cmax == nullptr)
    {
      throw std::logic_error("Wavespeed has not been set!");
    }
  this->dt = dt; // set it so that GetSlabHeight can return it
  auto &vmap = Tent::vmap;
  TentSlabPitcher * slabpitcher = [this]() ->TentSlabPitcher* {
    switch (this->method)
      {
      case ngstents::EVolGrad:
        return new VolumeGradientPitcher<DIM>(this->ma);
        break;
      case ngstents::EEdgeGrad:
        return new EdgeGradientPitcher<DIM>(this->ma);
      default:
        cout << "Trying to pitch tent without setting a pitching method." << endl;
        return nullptr;
        break;
      }
  }();
  if(!slabpitcher) return false;
  cout << "Created slab pitcher"<<endl;
  
  //map periodic vertices
  MapPeriodicVertices(ma);
  cout << "Mapped periodic vertices" << endl;
  //calc wavespeed for each element and perhaps other stuff (i..e, calculating edge gradients, checking fine edges, etc)
  BitArray fine_edges(ma->GetNEdges());
  slabpitcher->InitializeMeshData<DIM>(lh,fine_edges,cmax, calc_local_ct, global_ct);
  cout << "Initialised mesh data" << endl;
  //remove periodic edges
  RemovePeriodicEdges(ma,fine_edges);
  
  //compute neighbouring data
  TableCreator<int> create_v2e, create_v2v;
  for ( ; !create_v2e.Done(); create_v2e++, create_v2v++)
    {
      for (int e : IntRange (0, ma->GetNEdges()))
        if(fine_edges.Test(e))
          {
            auto vts = ma->GetEdgePNums (e);
            int v1 = vts[0], v2 = vts[1];
            //if v1 (or v2) is not periodic, vmap[v1] == v1
            create_v2v.Add (vmap[v1], v2);
            create_v2e.Add (vmap[v1], e);
            create_v2v.Add (vmap[v2], v1);
            create_v2e.Add (vmap[v2], e);
          }
    }
  auto v2v = create_v2v.MoveTable();
  auto v2e = create_v2e.MoveTable();

  cout << "Computed neighbouring data" << endl;
  Array<double> tau(ma->GetNV());  // advancing front values at vertices
  tau = 0.0;

  
  slabpitcher->ComputeVerticesReferenceHeight(v2v, v2e, tau, lh);
  cout << "Computed reference heights" << endl;
  // max time increase allowed at vertex, depends on tau of neighbors
  //at the beginning the advancing front is at a constant t=0
  //so ktilde can be set as vertex_refdt
  Array<double> ktilde = slabpitcher->GetVerticesReferenceHeight();
  
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
  //every vertex start at level 0
  vertices_level = 0;

  //for an advance to be considered good, dt >= factor * refdt
  double adv_factor{0.5};
  //whether to reset the adv_factor to its initial value after populating ready_vertices
  constexpr bool reset_adv_factor = true;
  // array of vertices ready for pitching a tent
  Array<int> ready_vertices;
  // array for checking if a given vertex is ready
  Array<bool> vertex_ready(ma->GetNV());
  vertex_ready = false;
  bool slab_complete{false};
  //array for checking if a given vertex is complete (tau[vi] = dt)
  Array<bool> complete_vertices(ma->GetNV());
  complete_vertices = false;
  //numerical tolerance
  const double num_tol = std::numeric_limits<double>::epsilon() * dt;
  while ( !slab_complete )
    {
      cout << "Setting ready vertices" << endl;
      const bool found_vertices =
        slabpitcher->GetReadyVertices(adv_factor,reset_adv_factor,ktilde,complete_vertices, vertex_ready,ready_vertices);
      //no possible vertex in which a tent could be pitched was found
      if(!found_vertices) break;
      // ---------------------------------------------
      // Main loop: constructs one tent each iteration
      // ---------------------------------------------
      cout << "Pitching tents..." << endl;
      while (ready_vertices.Size())
        {
          int minlevel, posmin;
          std::tie(minlevel,posmin) =
            slabpitcher->PickNextVertexForPitching(ready_vertices,ktilde,vertices_level);
          nlayers = max(minlevel,nlayers);
          //vertex index at which the current tent is being pitched
          const int vi = ready_vertices[posmin];
          ready_vertices.DeleteElement(posmin);
          vertex_ready[vi] = false;

          //current tent
          Tent * tent = new Tent;
          tent->vertex = vi;
          tent->tbot = tau[vi];

          const auto new_ttop = tau[vi] + ktilde[vi];
          if(dt - new_ttop > num_tol)
            {//not close to the end of the time slab
              tent->ttop = new_ttop;
            }
          else
            {//vertex is complete
              tent->ttop = dt;
              complete_vertices[vi] = true;
            }
          //let us ignore this for now
          // else if(new_ttop >= dt)
          //   {//vertex is complete
          //     tent->ttop = dt;
          //     complete_vertices[vi] = true;
          //   }
          // else
          //   {//vertex is really close to the end of time slab.
          //     //in this scenario, we might want to pitch a lower
          //     //tent to avoid numerical issues with degenerate tents
          //     tent->ttop = ktilde[vi] * 0.75 + tau[vi];
          //   }
          
          tent->level = vertices_level[vi]; // 0;
          tau[vi] = tent->ttop;
          ktilde[vi] = 0;//assuming that ktilde[vi] was the maximum advance

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

          slabpitcher->UpdateNeighbours(vi,adv_factor,v2v,v2e,tau,complete_vertices,
                                        ktilde,vertex_ready,ready_vertices,lh);
          tents.Append (tent);
        }
      //check if slab is complete
      slab_complete = true;
      for(int i = 0; i < ma->GetNV(); i++)
        if(vmap[i] == i)
          if(complete_vertices[i] == false)
            {
              slab_complete = false;
              break;
            }
    }

 
  if(!slab_complete)
    {
      cout << "Error: the algorithm could not pitch the whole slab" << endl;
      int iv;
      for(iv = 0; iv < ma->GetNV(); iv++)
        if(vmap[iv] == iv && !complete_vertices[iv]) break;
      if(iv == ma->GetNV())//just as a precaution, let us check that it really didnt pitch.
        {
          cout << "Inconsistent data structure. Aborting..." << endl;
          exit(-1);
        }
    }
  delete slabpitcher;
  

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
         constexpr auto el_type = EL_TYPE(DIM);
         //number of vertices of the current element (always the simplex associated to DIM)
         constexpr int n_vertices = DIM+1;
         //finite element created for calculating the barycentric coordinates
         ScalarFE<el_type,1> fe;
	 Vector<> shape_nodal(n_vertices);
	 Matrix<> dshape_nodal(n_vertices, DIM);
	 Vector<> coef_bot, coef_top; // coefficient of tau (top & bot)
	 coef_bot.SetSize(n_vertices);
	 coef_top.SetSize(n_vertices);
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
	 fe.CalcMappedDShape(mip, dshape_nodal);
	 tent.gradphi_bot[j] = Trans(dshape_nodal) * coef_bot;
	 tent.gradphi_top[j] = Trans(dshape_nodal) * coef_top;
       }
     });
  has_been_pitched = slab_complete;
  return has_been_pitched;
}


template <int DIM>
double TentPitchedSlab <DIM>::MaxSlope() const{

  // Return  max(|| gradphi_top||, ||gradphi_bot||)

  double maxgrad = 0.0;
  ParallelFor
    (Range(tents),
     [&] (int i)
     {
       Tent & tent = *tents[i];
       AtomicMax(maxgrad , tent.MaxSlope() );
     });
  return maxgrad;
}


///////////////////// Pitching Algo Routines ///////////////////////////////
TentSlabPitcher::TentSlabPitcher(shared_ptr<MeshAccess> ama) : ma(ama), cmax(ama->GetNE()), vertex_refdt(ama->GetNV()), edge_len(ama->GetNEdges()), local_ctau([](const int, const int){return 1.;}) {
}


bool TentSlabPitcher::GetReadyVertices(double &adv_factor, bool reset_adv_factor,
                                       const Array<double> &ktilde, const Array<bool> &complete_vertices,
                                       Array<bool> &vertex_ready, Array<int> &ready_vertices){
  auto &vmap = Tent::vmap;
  bool found{false};
  //how many times the adv_factor will be relaxed looking for new vertices
  constexpr int n_attempts = 5;
  vertex_ready = false;
  const double initial_adv_factor = adv_factor;
  for(auto ia = 0; ia < n_attempts; ia++)
    {
      for (auto iv = 0; iv < ma->GetNV(); iv++)
        if(vmap[iv] == iv && !complete_vertices[iv] )
          {
            if (ktilde[iv] > adv_factor * vertex_refdt[iv])
              if (!vertex_ready[iv])
                {
                  ready_vertices.Append (iv);
                  vertex_ready[iv] = true;
                }
          }
      if(ready_vertices.Size())
        {
          found = true;
          break;
        }
      adv_factor /= 2;
    }
  if(reset_adv_factor)
    adv_factor = initial_adv_factor;
  //the algorithm is most likely stuck
  else if (adv_factor < 0.05) return false;
  return found;
}

void TentSlabPitcher::ComputeVerticesReferenceHeight(const Table<int> &v2v, const Table<int> &v2e, const Array<double> &tau, LocalHeap &lh)
{
  auto &vmap = Tent::vmap;
  this->vertex_refdt = std::numeric_limits<double>::max();
  for (auto i = 0; i < this->ma->GetNV(); i++)
    if(vmap[i]==i) // non-periodic
      {
        this->vertex_refdt[i] = this->GetPoleHeight(i, tau, this->cmax, v2v[i],v2e[i],lh);
      }
  
}

std::tuple<int,int> TentSlabPitcher::PickNextVertexForPitching(const Array<int> &ready_vertices,
                                                               const Array<double> &ktilde,
                                                               const Array<int> &vertices_level){  
  int minlevel = std::numeric_limits<int>::max();
  int posmin = -1;
  for(auto i = 0; i < ready_vertices.Size(); i++)
    if(vertices_level[ready_vertices[i]] < minlevel)
      {
        minlevel = vertices_level[ready_vertices[i]];
        posmin = i;
      }
  return std::make_tuple(minlevel,posmin);
}

void TentSlabPitcher::UpdateNeighbours(const int vi, const double adv_factor, const Table<int> &v2v,
                                       const Table<int> &v2e, const Array<double> &tau,
                                       const Array<bool> &complete_vertices, Array<double> &ktilde,
                                       Array<bool> &vertex_ready, Array<int> &ready_vertices,
                                       LocalHeap &lh){
  auto &vmap = Tent::vmap;
  for (int nb : v2v[vi])
    {
      nb = vmap[nb]; // map periodic vertices
      if (complete_vertices[nb]) continue;
      const double kt = GetPoleHeight(nb, tau, cmax, v2v[nb], v2e[nb],lh);
      ktilde[nb] = kt;
      if (kt > adv_factor * vertex_refdt[nb])
        {
          if (!vertex_ready[nb])
            {
              ready_vertices.Append (nb);
              vertex_ready[nb] = true;
            }
        }
      else
        {
          vertex_ready[nb] = false;
          const auto pos_nb = ready_vertices.Pos(nb);
          if(pos_nb != ready_vertices.ILLEGAL_POSITION)
            {
              ready_vertices.RemoveElement(pos_nb);
            }
        }
    } 
}

template<int DIM>
void TentSlabPitcher::InitializeMeshData(LocalHeap &lh, BitArray &fine_edges, shared_ptr<CoefficientFunction>wavespeed, bool calc_local_ct, const double global_ct)
{
  constexpr auto el_type = EL_TYPE(DIM);//simplex of dimension dim
  constexpr auto n_el_vertices = DIM + 1;//number of vertices of that simplex
  //sets global constant
  this->global_ctau = global_ct;
  HeapReset hr(lh);
  fine_edges.Clear();

  const auto n_vol_els = ma->Elements(VOL).Size();
  //this table will contain the local mesh-dependent constant
  //local_ctau is the ratio between the distance to the opposite facet
  //and the biggest edge to the opposite facet, so local_ctau <1
  TableCreator<double> create_local_ctau;
  create_local_ctau.SetSize(n_vol_els);
  create_local_ctau.SetMode(2);
  ArrayMem<int,n_el_vertices> el_vertices(n_el_vertices);
  //just calculating the size of the table
  for(auto el : IntRange(0,n_vol_els))
    {
      for(auto v : IntRange(0, n_el_vertices))
        create_local_ctau.Add(el,v);
    }
  create_local_ctau++;// it is in insert mode
   
  //gradient of basis functions onthe current element  
  Matrix<> gradphi(n_el_vertices, DIM);
  //used to calculate distance to opposite facet
  ScalarFE<el_type,1> my_fel;
  //minimum length of the adjacent edges for each element's vertices
  ArrayMem<double, n_el_vertices> max_edge(n_el_vertices);
  //the mesh contains only simplices so only one integration rule is needed
  IntegrationRule ir(el_type, 0);
  for (Ngs_Element el : this->ma->Elements(VOL))
    {
      max_edge = -1;
      auto ei = ElementId(el);
      ElementTransformation & trafo = this->ma->GetTrafo (ei, lh);
      MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
      this->cmax[el.Nr()] = wavespeed->Evaluate(mip);


      auto v_indices = ma->GetElVertices(ei);
      //set all edges belonging to the mesh
      for (int e : el.Edges())
        {
          auto pnts = ma->GetEdgePNums(e);
          auto v1 = pnts[0], v2 = pnts[1];
          if(!fine_edges[e])
            {
              fine_edges.SetBit(e);
              double len = L2Norm (ma-> template GetPoint<DIM>(v1)
                               - ma-> template GetPoint<DIM>(v2));
              edge_len[e] = len;
            }
          const auto v1_local = v_indices.Pos(v1);
          const auto v2_local = v_indices.Pos(v2);
          max_edge[v1_local] = max(max_edge[v1_local],edge_len[e]);
          max_edge[v2_local] = max(max_edge[v2_local],edge_len[e]);
        }
      my_fel.CalcMappedDShape(mip,gradphi);
      const auto el_num = ei.Nr();
      const auto detjac_inv = 1./mip.GetJacobiDet();
      for(int vi_local = 0; vi_local < v_indices.Size(); vi_local++) 
        {
          const auto dist_opposite_facet = 1./L2Norm(gradphi.Row(vi_local));
          const auto val = dist_opposite_facet / max_edge[vi_local];
          create_local_ctau.Add(el_num,val);
        }
    }

  if(calc_local_ct && DIM > 1)
    {
      local_ctau_table = create_local_ctau.MoveTable();
      this->local_ctau = [this](const int el, const int v){return local_ctau_table[el][v];};
    }
  else
    {
      this->local_ctau = [](const int el, const int v){return 1;};
    }
  
}

template <int DIM> double VolumeGradientPitcher<DIM>::GetPoleHeight(const int vi, const Array<double> & tau, const Array<double> & cmax, FlatArray<int> nbv, FlatArray<int> nbe, LocalHeap & lh) const{
  HeapReset hr(lh);
  constexpr auto el_type = EL_TYPE(DIM);
  //number of vertices of the current element (always the simplex associated to DIM)
  constexpr int n_vertices = DIM+1;
  //finite element created for calculating the barycentric coordinates
  ScalarFE<el_type,1> my_fel;
  // array of all elements containing vertex vi
  ArrayMem<int,30>  els;
  els.SetSize(0);
  this->ma->GetVertexElements(vi, els);

  constexpr double init_pole_height = std::numeric_limits<double>::max();
  double pole_height = init_pole_height;
  //vector containing the advancing front time for each vertex (except vi)
  Vector<double> coeff_vec(n_vertices);
  //gradient of basis functions onthe current element
  Matrix<double> gradphi(n_vertices,DIM);
  //numerical tolerance (NOT YET SCALED)
  double num_tol = std::numeric_limits<double>::epsilon();

  double det_jac_inv = -1; 
  for (int el : els)
    {
      ElementId ei(VOL,el);
      //c_max^2
      const double c_max_sq = cmax[ei.Nr()] * cmax[ei.Nr()]; 
      //mapping of the current el
      ElementTransformation &trafo = this->ma->GetTrafo(ei, lh);
      //vertices of current el
      auto v_indices = this->ma->GetElVertices(ei);
      //vi position in current el
      const auto local_vi = v_indices.Pos(vi);
      //integration rule for reference el
      IntegrationRule ir(el_type,1);
      //integration point on deformed element
      MappedIntegrationPoint<DIM,DIM> mip(ir[0],trafo);
      det_jac_inv = max (det_jac_inv, 1/mip.GetJacobiDet());
      my_fel.CalcMappedDShape(mip,gradphi);

      //sets the coefficient vec
      for (auto k : IntRange(0, v_indices.Size()))
        coeff_vec(k) = tau[v_indices[k]];
      coeff_vec[local_vi] = 0;
      
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
      //since alpha>0 we can scale the equation by alpha
      Vec<DIM> tau_i_grad_i = Trans(gradphi) * coeff_vec;
      const double beta = 2 * InnerProduct(tau_i_grad_i,gradphi.Row(local_vi))/alpha;
      const double gamma = (InnerProduct(tau_i_grad_i, tau_i_grad_i) - 1.0/c_max_sq)/alpha;
      const double delta = beta * beta - 4 * gamma;

      //since sq_delta will always be positive
      //we dont need to consider the solution (-beta-sq_delta)/(2*alpha)
      const double sol = [alpha,beta, gamma, delta,init_pole_height,num_tol](){
        if(delta > num_tol * alpha)//positive delta
          {
            if(beta <= num_tol * alpha)
              return (sqrt(delta)-beta)/2.0;
            else
              return - (2.0 * gamma) / (beta + sqrt(delta));
          }
        else if (delta > -num_tol*alpha)//zero delta
          {
            return -beta/2.0;
          }
        return init_pole_height;//negative delta
      }();
      pole_height = min(pole_height,this->global_ctau * sol);
    }

  //scaling of numerical tolerance
  num_tol *= det_jac_inv;
  //check if a real solution to the quadratic equation was found
  if(fabs(pole_height - init_pole_height) < num_tol) return 0.0;
  //the return value is actually the ADVANCE in the current vi
  pole_height -= tau[vi];
  if( pole_height < num_tol ) return 0.0;
  return pole_height - num_tol;//just to enforce causality
 }

template <int DIM>
double EdgeGradientPitcher<DIM>::GetPoleHeight(const int vi, const Array<double> & tau, const Array<double> & cmax, FlatArray<int> nbv, FlatArray<int> nbe, LocalHeap & lh) const{
  auto &vmap = Tent::vmap;
  double num_tol = std::numeric_limits<double>::epsilon();
  double kt = std::numeric_limits<double>::max();

  // array of all elements containing vertex vi
  ArrayMem<int,30> els;
  els.SetSize(0);
  ma->GetVertexElements(vi,els);
  for (int nb_index : nbv.Range())
    {
      int nb = vmap[nbv[nb_index]];
      const double length = edge_len[nbe[nb_index]];
      for(int el : els)
        {
          ElementId ei(VOL,el);
          const auto el_num = ei.Nr();
          const auto vi_local = ma->GetElVertices(ei).Pos(vi);
          const double c_max = cmax[el_num];
          const double local_ct = this->local_ctau(el_num,vi_local);
          const double kt1 = tau[nb]-tau[vi]+ global_ctau * local_ct * length/c_max;
          if (kt1 > 0)
            {          
              kt = min (kt, kt1);
            }
          else continue;
        }
      
      
    }
  //could not advance
  if (fabs(kt-std::numeric_limits<double>::max()) < 1) return 0.0;
  //ensuring causality (inequality)
  kt -= num_tol;
  if(kt > num_tol) return kt;
  else return 0.0;
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
    .def_readonly("internal_facets", &Tent::internal_facets)
    .def("MaxSlope", &Tent::MaxSlope);

  //
  // 1D spatial mesh
  //
  py::class_<TentPitchedSlab<1>, shared_ptr<TentPitchedSlab<1>>>
    (m, "TentPitchedSlab1", "Tent pitched slab in 1 space + 1 time dimensions")
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
        auto tps = TentPitchedSlab<1>(ma,  heapsize);
        tps.SetPitchingMethod(method);
        return tps;
      }),
      py::arg("mesh"), py::arg("method"), py::arg("heapsize") = 1000000
      )
    
    .def_readonly("mesh", &TentPitchedSlab<1>::ma)
    .def("SetWavespeed", static_cast<void (TentPitchedSlab<1>::*)(const double)>(&TentPitchedSlab<1>::SetWavespeed))
    .def("PitchTents", &TentPitchedSlab<1>::PitchTents, py::arg("dt"), py::arg("local_ct"), py::arg("global_ct")= 1.0)
    .def("GetNTents", &TentPitchedSlab<1>::GetNTents)
    .def("GetNLayers", &TentPitchedSlab<1>::GetNLayers)
    .def("GetSlabHeight", &TentPitchedSlab<1>::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab<1>::MaxSlope)
    .def("GetTent", &TentPitchedSlab<1>::GetTent, pybind11::return_value_policy::reference_internal)
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
        auto tps = TentPitchedSlab<2>(ma,  heapsize);
        tps.SetPitchingMethod(method);
        return tps;
      }),
      py::arg("mesh"), py::arg("method"), py::arg("heapsize") = 1000000
      )

    .def_readonly("mesh", &TentPitchedSlab<2>::ma)
    .def("SetWavespeed", static_cast<void (TentPitchedSlab<2>::*)(const double)>(&TentPitchedSlab<2>::SetWavespeed))
    .def("PitchTents", &TentPitchedSlab<2>::PitchTents, py::arg("dt"), py::arg("local_ct"), py::arg("global_ct")= 1.0)
    .def("GetNTents", &TentPitchedSlab<2>::GetNTents)
    .def("GetNLayers", &TentPitchedSlab<2>::GetNLayers)
    .def("GetSlabHeight", &TentPitchedSlab<2>::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab<2>::MaxSlope)
    .def("GetTent", &TentPitchedSlab<2>::GetTent, pybind11::return_value_policy::reference_internal)
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
        auto tps = TentPitchedSlab<3>(ma,  heapsize);
        tps.SetPitchingMethod(method);
        return tps;
      }),
      py::arg("mesh"), py::arg("method"), py::arg("heapsize") = 1000000
      )

    .def_readonly("mesh", &TentPitchedSlab<3>::ma)
    .def("SetWavespeed", static_cast<void (TentPitchedSlab<3>::*)(const double)>(&TentPitchedSlab<3>::SetWavespeed))
    .def("PitchTents", &TentPitchedSlab<3>::PitchTents, py::arg("dt"), py::arg("local_ct"), py::arg("global_ct")= 1.0)
    .def("GetNTents", &TentPitchedSlab<3>::GetNTents)
    .def("GetNLayers", &TentPitchedSlab<3>::GetNLayers)
    .def("GetSlabHeight", &TentPitchedSlab<3>::GetSlabHeight)
    .def("MaxSlope", &TentPitchedSlab<3>::MaxSlope)
    .def("GetTent", &TentPitchedSlab<3>::GetTent, pybind11::return_value_policy::reference_internal)
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
