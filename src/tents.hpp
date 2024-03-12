#ifndef TENTSHEADER
#define TENTSHEADER

#ifdef WIN32
        #define NGTENT_API_EXPORT __declspec(dllexport)
        #define NGTENT_API_IMPORT __declspec(dllimport)
#else
        #define NGTENT_API_EXPORT
        #define NGTENT_API_IMPORT
#endif

#ifdef NGSTENT_EXPORTS
        #define NGSTENT_API NGTENT_API_EXPORT
#else
        #define NGSTENT_API NGTENT_API_IMPORT
#endif

#include <solve.hpp>
#include <h1lofe.hpp> // seems needed for ScalarFE (post 2021-06-22 NGSolve update)
using namespace ngsolve;
using namespace std;


////////////////////////////////////////////////////////////////////////////
///
/// Class to describe one spacetime tent in a mesh of tents
///
/// A spacetime tent is a macroelement consisting of a tentpole erected at
/// a central vertex in space and all the space-time tetrahedra with
/// the tentpole as an edge.
///
/// We represent the tent by its projection on space (a vertex patch),
/// the central vertex, and the heights (times) of its neighboring
/// vertices.
///


class NGSTENT_API Tent {

public:
  Tent(const Array<int> &avmap) : vmap(avmap){}
  Tent() = delete;
  int vertex;                 ///< central vertex
  double tbot, ttop;          ///< bottom and top times of central vertex
  Array<int> nbv;             ///< neighbouring vertices of central vertex
  Array<double> nbtime;       ///< height/time of neighbouring vertices
  Array<int> els;             ///< all elements in the tent's vertex patch
  Array<int> internal_facets; ///< all internal facets in the tent's vertex patch

  /// elfnums[k] lists all internal facets of the k-th element of tent
  Table<int> elfnums;         

  const Array<int> &vmap;     ///< vertex map for any periodicity identification

  /// access to the finite element & dofs
  mutable class TentDataFE * fedata = nullptr;

  int level;                  ///< my parallel layer number in a mesh of tents
  Array<int> dependent_tents; ///< these tents depend on me

  double maxslope = 0.0;      ///< maximal slope of the top advancing front
  double MaxSlope() const { return maxslope; }

  /// global physical time at vertex (stored in ConservationLaw::gftau)
  mutable double * time;     
  mutable double timebot;     ///< global physical bottom time at vertex

  void InitTent(shared_ptr<GridFunction> gftau) const
  {
    time = &(gftau->GetVector().FVDouble()(vertex));
    timebot = *time;
  }

  void SetFinalTime() const { *time = timebot + (ttop - tbot); }
};

ostream & operator<< (ostream & ost, const Tent & tent);


////////////////////////////////////////////////////////////////////////////
///
/// Class with dofs, finite element & integration info for a tent:
///
class NGSTENT_API TentDataFE
{
public:
  int nd;           ///< total # interior and interface dofs in space
  Array<int> dofs;  ///< all interior and interface dof nums, size(dofs)=nd.
  /// ranges[k] = IntRange (of dof numbers) of k-th element of local matrix
  Array<IntRange> ranges;
  /// finite elements for all elements in the tent
  Array<FiniteElement*> fei;
  /// integration rules for all elements in the tent
  Array<SIMD_IntegrationRule*> iri;
  /// mapped integration rules for all elements in the tent
  Array<SIMD_BaseMappedIntegrationRule*> miri;
  /// element transformations for all elements in the tent
  Array<ElementTransformation*> trafoi;
  /// mesh size for each element
  Array<double> mesh_size;
  //// gradients of tent bottom at integration points (in possibly curved elements)
  Array<FlatMatrix<SIMD<double>>> agradphi_bot;
  /// gradient of (tent top) the new advancing front in the IP's
  Array<FlatMatrix<SIMD<double>>> agradphi_top;
  /// height of the tent in the IP's
  Array<FlatVector<SIMD<double>>> adelta;
  /// local numbers of the neighbors
  Array<ngcore::IVec<2,size_t>> felpos;
  /// facet integration rules for all facets in the tent
  Array<SIMD_IntegrationRule*> fir;
  /// facet integration rules for all internal facets in the tent
  /// transformed to local coordinates of the two neighboring elements
  Array<Vec<2,const SIMD_IntegrationRule*>> firi;
  /// mapped facet integration rules for all facets
  Array<SIMD_BaseMappedIntegrationRule*> mfiri1;
  /// mapped facet integration rules for all facets
  Array<SIMD_BaseMappedIntegrationRule*> mfiri2;
  /// gradient phi face first and second element
  Array<FlatMatrix<SIMD<double>>> agradphi_botf1;
  Array<FlatMatrix<SIMD<double>>> agradphi_topf1;
  Array<FlatMatrix<SIMD<double>>> agradphi_botf2;
  Array<FlatMatrix<SIMD<double>>> agradphi_topf2;
  /// normal vectors in the IP's
  Array<FlatMatrix<SIMD<double>>> anormals;
  /// height of the tent in the IP's
  Array<FlatVector<SIMD<double>>> adelta_facet;

  TentDataFE(const Tent & tent, const FESpace & fes, LocalHeap & lh);
};

////////////////////////////////////////////////////////////////////////////
///
/// Class representing the spatial gradient of the advancing front
/// φ(x, τ) = (1-τ) φbot(x) + τ φtop(x)
///

class GradPhiCoefficientFunction : public CoefficientFunction
{
public:
  GradPhiCoefficientFunction (int adim)
    : CoefficientFunction(adim)
  { }

  double Evaluate(const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("Evaluate not implemented for BaseMappedIntegrationPoint!");
  }

  /// Return "values" of grad(φ) at points of a mapped integration rule "mir".
  /// Since grad(φ) depends on x and pseudotime τ, this evaluation makes sense 
  /// only when τ is known. Whoever needs grad(φ) at some known τ should compute
  /// it and send it down as ProxyUserData of mir's trafo. This member function
  /// simply copies and returns these precomputed values – no actual computation
  /// is done here.
  void Evaluate(const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> values) const  {
    ProxyUserData & ud =
      *static_cast<ProxyUserData*>(mir.GetTransformation().userdata);
    values.AddSize(Dimension(), mir.Size()) =
      BareSliceMatrix<SIMD<double>> (ud.GetAMemory (this));
  }

  void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;

  shared_ptr<CoefficientFunction>
  Diff(const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const
  {
    if(var == this)
      return dir;
    else
      return ZeroCF(Dimensions());
  }
};

////////////////////////////////////////////////////////////////////////////
namespace ngstents{
  enum PitchingMethod {EVolGrad =1, EEdgeGrad};
}

class NGSTENT_API TentPitchedSlab {
protected:
  double dt;                              // time step between two time slices
  shared_ptr<CoefficientFunction> cmax;   // wavespeed
  ngstents::PitchingMethod method;
  bool has_been_pitched;                  // whether the slab has been already pitched
  Array<Tent*> tents;                     // tents between two time slices
  int nlayers;                            // number of layers in the time slab

  Array<int> vmap;                        // vertex map for periodic boundaries
  LocalHeap lh;

public:
  // access to base spatial mesh (public for export to Python visualization)
  shared_ptr<MeshAccess> ma;
  // Propagate methods need access to DAG of tent dependencies
  Table<int> tent_dependency;
  // access to grad(phi) coefficient function
  shared_ptr<CoefficientFunction> cfgradphi = nullptr;

  // Constructor and initializers
  TentPitchedSlab(shared_ptr<MeshAccess> ama, int heapsize) :
    dt(0), ma(ama), cmax(nullptr), nlayers(0),
    has_been_pitched(false), lh(heapsize, "Tents heap")
  {
    cfgradphi = make_shared<GradPhiCoefficientFunction>(ma->GetDimension());
  };
  
  //uses a gradient based method for pitching the tent
  //calc_local_ct will indicate whether to use a local mesh-dependent
  //constant for the algorithm
  //global_ct is a globalwise constant that can be independently used
  
  template <int DIM>
  bool PitchTents(const double dt, const bool calc_local_ct, const double global_ct = 1.0);
  
  // Get object features
  int GetNTents() { return tents.Size(); }
  int GetNLayers() { return nlayers + 1; }

  void SetMaxWavespeed(const double c){cmax =  make_shared<ConstantCoefficientFunction>(c);}
  void SetMaxWavespeed(shared_ptr<CoefficientFunction> c){ cmax = c;}
  
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
};

//Abstract class with the interface of methods used for pitching a tent
class NGSTENT_API TentSlabPitcher{
protected:
  //access to base spatial mesh
  shared_ptr<MeshAccess> ma;
  //element-wise (vol grad algo) or edge-wise (edge grad algo)  maximal wave-speeds
  Array<double> cmax;
  //reference heights for each vertex
  Array<double> vertex_refdt;
  //array containing the length of each edge
  Array<double> edge_len;
  //returns the mesh dependent local constants. the first parameter is the vertex,
  //and the second parameter is the edge (edge algo) or element (vol algo)
  std::function<double(const int, const int)> local_ctau;
  //table for storing local geometric constants
  Table<double> local_ctau_table;
  //global constant (defaulted to 1)
  double global_ctau;
  const ngstents::PitchingMethod method;
  
  // Calculate c_tau to ensure causality (edge algo) / prevent locks (vol algo)
  virtual Table<double> CalcLocalCTau(LocalHeap& lh, const Table<int> &v2e) = 0;



public:
  
  TentSlabPitcher(shared_ptr<MeshAccess> ama,
		  ngstents::PitchingMethod m, Array<int> &avmap);
  
  virtual ~TentSlabPitcher(){;}

  // Precompute mesh-dependent data, including the wavespeed (per
  // element) and neighbouring data. It returns the table v2v
  // (neighbouring vertices), v2e(edges adjacent to a given vertex)
  // and per_verts (used for periodicity).
  
  template<int DIM> std::tuple<Table<int>,Table<int>>
  InitializeMeshData(LocalHeap &lh,
		     shared_ptr<CoefficientFunction> wavespeed,
		     bool calc_local_ctau, const double global_ct );

  // Compute vertex based max time-differences assumint tau=0
  // corresponding to a non-periodic vertex

  void ComputeVerticesReferenceHeight(const Table<int> &v2v,
				      const Table<int> &v2e,
				      const FlatArray<double> &tau, LocalHeap &lh);
  
  void UpdateNeighbours(const int vi, const double adv_factor,
			const Table<int> &v2v,const Table<int> &v2e,
                        const FlatArray<double> &tau,
			const BitArray &complete_vertices,
                        Array<double> &ktilde, BitArray &vertex_ready,
                        Array<int> &ready_vertices, LocalHeap &lh);
  
  // Return a copy of vertex_refdt  (without any computation)
  Array<double> GetVerticesReferenceHeight(){ return Array<double>(vertex_refdt);}

  // Populate the set of ready vertices with vertices satisfying
  //   ktilde > adv_factor * refdt. Returns false if no such vertex was found.
  [[nodiscard]] bool
  GetReadyVertices(double &adv_factor, bool reset_adv_factor,
                                      const FlatArray<double> &ktilde,
				      const BitArray &complete_vertices,
                                      BitArray &vertex_ready,
				      Array<int> &ready_vertices);

  // Given the current advancing (time) front, calculates the maximum
  // advance on a tent centered on vi that will still guarantee causality

  virtual double
  GetPoleHeight(const int vi, const FlatArray<double> & tau,
		FlatArray<int> nbv, FlatArray<int> nbe, LocalHeap & lh) const = 0;

  //Returns the position in ready_vertices containing the vertex in which a tent will be pitched (and its level)
  [[nodiscard]] std::tuple<int,int> PickNextVertexForPitching(const FlatArray<int> &ready_vertices, const FlatArray<double> &ktilde, const FlatArray<int> &vertices_level);

  //////////////// For handling periodicity //////////////////////////////////

  // Get all elements connected to a given vertex (contemplating periodicity)
  void GetVertexElements(int vnr_master, Array<int> & elems) const;

  // Get all elements connected to a given edge (contemplating periodicity)
  void GetEdgeElements(int edge, Array<int> & elems) const;
  
  void MapPeriodicVertices();


  void RemovePeriodicEdges(BitArray &fine_edges);

  // access to global periodicity identifications
  Array<int> &vmap;      // vertex map for periodic spaces
  //per_verts[v] will contain all periodic
  // vertices associated with v
  Table<int> per_verts;

};

template <int DIM>
class VolumeGradientPitcher : public TentSlabPitcher{
public:
  
  VolumeGradientPitcher(shared_ptr<MeshAccess> ama, Array<int> &avmap)
    : TentSlabPitcher(ama, ngstents::PitchingMethod::EVolGrad, avmap){;}

  double GetPoleHeight(const int vi, const FlatArray<double> & tau,
		       FlatArray<int> nbv,
                       FlatArray<int> nbe, LocalHeap & lh) const override;

  // Calculate c_tau to prevent locks 

  Table<double> CalcLocalCTau(LocalHeap& lh, const Table<int> &v2e) override;
};

template <int DIM>
class EdgeGradientPitcher : public TentSlabPitcher{
public:
  
  EdgeGradientPitcher(shared_ptr<MeshAccess> ama, Array<int> &avmap)
    : TentSlabPitcher(ama, ngstents::PitchingMethod::EEdgeGrad, avmap) {;}

  double GetPoleHeight(const int vi, const FlatArray<double> & tau,
		       FlatArray<int> nbv,
                       FlatArray<int> nbe, LocalHeap & lh) const override;

  // Calculate c_tau to ensure causality

  Table<double> CalcLocalCTau(LocalHeap& lh, const Table<int> &v2e) override;
  
};
#endif
