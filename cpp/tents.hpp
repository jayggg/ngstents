#ifndef TENTSHEADER
#define TENTSHEADER

#include <solve.hpp>
using namespace ngsolve;
using namespace std;

// A spacetime tent is a macroelement consisting of a tentpole erected at
// a central vertex in space and all the space-time tetrahedra with
// the tentpole as an edge.

// We represent the tent by its projection on space (a vertex patch),
// the central vertex, and the heights (times) of its neighboring
// vertices.

////////////////////////////////////////////////////////////////////////////

// Class to describe one spacetime tent in a mesh of tents

class Tent {

public:

  int vertex;                 // central vertex
  double tbot, ttop;          // bottom and top times of central vertex
  Array<int> nbv;             // neighbour vertices
  Array<double> nbtime;       // height/time of neighbouring vertices
  Array<int> els;             // all elements in the tent's vertex patch
  Array<int> internal_facets; // all internal facets in the tent's vertex patch
  Table<int> elfnums;         /* elfnums[k] lists all internal facets of
				 the k-th element of tent */

  // Tent top and bottom are graphs of phi_top, phi_bot, which are
  // p.w.linear functions on non-curved elements (with p.w.constant gradients).
  Array<Vector<>> gradphi_bot; // gradphi_bot[l], gradphi_top[l] =
  Array<Vector<>> gradphi_top; /* gradients of phi_bot/top at some point in the
				  l-th simplex of the tent */

  // access to global periodicity identifications
  static Array<int> vmap;      // vertex map for periodic spaces

  // access to the finite element & dofs
  mutable class TentDataFE * fedata = nullptr;

  // other global details from a mesh of tents
  int level;                   // parallel layer number
  Array<int> dependent_tents;  // these tents depend on me

  double MaxSlope() const;
};

ostream & operator<< (ostream & ost, const Tent & tent);


////////////////////////////////////////////////////////////////////////////

// Class with dofs, finite element & integration info for a tent:

class TentDataFE
{
public:
  // moved from Tent
  int nd;            // total # interior and interface dofs in space
  Array<int> dofs;   // all interior and interface dof nums, size(dofs)=nd.
  // ranges[k] IntRange (of dof numbers) of k-th element of local matrix
  Array<IntRange> ranges;

  // finite elements for all elements in the tent
  Array<FiniteElement*> fei;
  // integration rules for all elements in the tent
  Array<SIMD_IntegrationRule*> iri;
  // mapped integration rules for all elements in the tent
  Array<SIMD_BaseMappedIntegrationRule*> miri;
  // element transformations for all elements in the tent
  Array<ElementTransformation*> trafoi;
  // mesh size for each element
  Array<double> mesh_size;
  // gradients of tent bottom at integration points (in possibly curved elements)
  Array<FlatMatrix<SIMD<double>>> agradphi_bot;
  // gradient of (tent top) the new advancing front in the IP's
  Array<FlatMatrix<SIMD<double>>> agradphi_top;
  // height of the tent in the IP's
  Array<FlatVector<SIMD<double>>> adelta;
  // local numbers of the neighbors
  Array<INT<2,size_t>> felpos;
  // facet integration rules for all facets in the tent
  // transformed to local coordinated of the neighboring elements
  Array<Vec<2,const SIMD_IntegrationRule*>> firi;
  // mapped facet integration rules for all facets
  Array<SIMD_BaseMappedIntegrationRule*> mfiri1;
  // mapped facet integration rules for all facets
  Array<SIMD_BaseMappedIntegrationRule*> mfiri2;
  // gradient phi face first and second element
  Array<FlatMatrix<SIMD<double>>> agradphi_botf1;
  Array<FlatMatrix<SIMD<double>>> agradphi_topf1;
  Array<FlatMatrix<SIMD<double>>> agradphi_botf2;
  Array<FlatMatrix<SIMD<double>>> agradphi_topf2;
  // normal vectors in the IP's
  Array<FlatMatrix<SIMD<double>>> anormals;
  // height of the tent in the IP's
  Array<FlatVector<SIMD<double>>> adelta_facet;

  TentDataFE(int n, LocalHeap & lh)
    : fei(n, lh), iri(n, lh), miri(n, lh), trafoi(n, lh) { ; }

  TentDataFE(const Tent & tent, const FESpace & fes,
	     const MeshAccess & ma, LocalHeap & lh);
};


////////////////////////////////////////////////////////////////////////////
namespace ngstents{
  enum PitchingMethod {EVolGrad =1, EEdgeGrad};
}

template <int DIM>
class TentPitchedSlab {
public:
  Array<Tent*> tents;         // tents between two time slices
  double dt;                  // time step between two time slices
  //wavespeed
  shared_ptr<CoefficientFunction> cmax;
  int nlayers;//number of layers in the time slab
  //whether the slab has been already pitched
  bool has_been_pitched;
  LocalHeap lh;
  ngstents::PitchingMethod method;

public:
  // access to base spatial mesh (public for export to Python visualization)
  shared_ptr<MeshAccess> ma;
  // Constructor and initializers
  TentPitchedSlab(shared_ptr<MeshAccess> ama, int heapsize) :
    dt(0), ma(ama), cmax(nullptr), nlayers(0),
    has_been_pitched(false), lh(heapsize, "Tents heap") { ; };
  
  //uses a gradient based method for pitching the tent
  //calc_local_ct will indicate wether to use a local mesh-dependent
  //constant for the algorithm
  //global_ct is a globalwise constant that can be independently used
  //its return value will indicate whether the slab was successfully pitched.
  bool PitchTents(const double dt, const bool calc_local_ct, const double global_ct = 1.0);
  
  // Get object features
  int GetNTents() { return tents.Size(); }
  int GetNLayers() { return nlayers; }


  void SetWavespeed(const double c){cmax =  make_shared<ConstantCoefficientFunction>(c);}
  void SetWavespeed(shared_ptr<CoefficientFunction> c){ cmax = c;}
  
  double GetSlabHeight() { return dt; }
  const Tent & GetTent(int i) { return *tents[i];}

  // Return  max(|| gradphi_top||, ||gradphi_bot||)
  double MaxSlope() const;

  // Drawing
  void DrawPitchedTents(int level=1) ;
  void DrawPitchedTentsVTK(string vtkfilename);
  void DrawPitchedTentsGL(Array<int> & tentdata,
                          Array<double> & tenttimes, int & nlevels);

  void SetPitchingMethod(ngstents::PitchingMethod amethod) {this->method = amethod;}

  // Propagate methods need to access this somehow
  Table<int> tent_dependency; // DAG of tent dependencies
};

//Abstract class with the interface of methods used for pitching a tent
class TentSlabPitcher{
protected:
  //access to base spatial mesh
  shared_ptr<MeshAccess> ma;
  //element-wise maximal wave-speeds
  Array<double> cmax;
  //reference heights for each vertex
  Array<double> vertex_refdt;
  //array containing the length of each edge
  Array<double> edge_len;
  //returns the constant associated with a given vertex of a given element
  std::function<double(const int, const int)> local_ctau;
  //table for storing local geometric constants
  Table<double> local_ctau_table;
  //global constant (defaulted to 1)
  double global_ctau;
public:
  //constructor
  TentSlabPitcher(shared_ptr<MeshAccess> ama);
  //destructor
  virtual ~TentSlabPitcher(){;}
  //calculates the wavespeed for each element and the edge length
  template<int DIM>
  void InitializeMeshData(LocalHeap &lh, BitArray &fine_edges,
                          shared_ptr<CoefficientFunction> wavespeed, bool calc_local_ctau, const double global_ct );

  //compute the vertex based max time-differences assumint tau=0
  //corresponding to a non-periodic vertex
  void ComputeVerticesReferenceHeight(const Table<int> &v2v, const Table<int> &v2e, const FlatArray<double> &tau,
                                      LocalHeap &lh);

  void UpdateNeighbours(const int vi, const double adv_factor, const Table<int> &v2v,const Table<int> &v2e,
                        const FlatArray<double> &tau, const FlatArray<bool> &complete_vertices,
                        Array<double> &ktilde, Array<bool> &vertex_ready,
                        Array<int> &ready_vertices, LocalHeap &lh);
  
  //it does NOT compute, only returns a copy of vertex_refdt
  Array<double> GetVerticesReferenceHeight(){ return Array<double>(vertex_refdt);}

  //Populate the set of ready vertices with vertices satisfying ktilde > adv_factor * refdt. returns false if
  //no such vertex was found
  [[nodiscard]] bool GetReadyVertices(double &adv_factor, bool reset_adv_factor,
                                      const FlatArray<double> &ktilde, const FlatArray<bool> &complete_vertices,
                                      Array<bool> &vertex_ready, Array<int> &ready_vertices);

  //Given the current advancing (time) front, calculates the
  //maximum advance on a tent centered on vi that will still
  //guarantee causality
  virtual double GetPoleHeight(const int vi, const FlatArray<double> & tau,  FlatArray<int> nbv, FlatArray<int> nbe, LocalHeap & lh) const = 0;

  //Returns the position in ready_vertices containing the vertex in which a tent will be pitched (and its level)
  [[nodiscard]] std::tuple<int,int> PickNextVertexForPitching(const FlatArray<int> &ready_vertices, const FlatArray<double> &ktilde, const FlatArray<int> &vertices_level);
};

template <int DIM>
class VolumeGradientPitcher : public TentSlabPitcher{
public:
  
  VolumeGradientPitcher(shared_ptr<MeshAccess> ama) : TentSlabPitcher(ama){;}

  double GetPoleHeight(const int vi, const FlatArray<double> & tau, FlatArray<int> nbv,
                       FlatArray<int> nbe, LocalHeap & lh) const override;
};

template <int DIM>
class EdgeGradientPitcher : public TentSlabPitcher{
public:
  
  EdgeGradientPitcher(shared_ptr<MeshAccess> ama) : TentSlabPitcher(ama) {;}

  double GetPoleHeight(const int vi, const FlatArray<double> & tau, FlatArray<int> nbv,
                       FlatArray<int> nbe, LocalHeap & lh) const override;
};
#endif
