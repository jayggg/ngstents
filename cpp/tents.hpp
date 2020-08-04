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

template <int DIM>
class TentPitchedSlab {

  Array<Tent*> tents;         // tents between two time slices
  double dt;                  // time step between two time slices
  LocalHeap lh;

public:
  // access to base spatial mesh (public for export to Python visualization)
  shared_ptr<MeshAccess> ma;
  // Constructor and initializers
  TentPitchedSlab(shared_ptr<MeshAccess> ama, int heapsize) :
    dt(0), ma(ama), lh(heapsize, "Tents heap") { ; };
  void PitchTents(double dt, double cmax);
  void PitchTents(double dt, shared_ptr<CoefficientFunction> cmax);

  //Given the current advancing (time) front, calculates the
  //maximum advance on a tent centered on vi that will still
  //guarantee causality
  double GetPoleHeight(const int vi, const Array<double> & tau, const Array<double> & cmax, LocalHeap & lh) const;
  
  //uses the exact calculation of the gradient for pitching the tent
  void PitchTentsGradient(double dt, double cmax);
  void PitchTentsGradient(double dt, shared_ptr<CoefficientFunction> cmax);
  // Get object features
  int GetNTents() { return tents.Size(); }
  double GetSlabHeight() { return dt; }
  const Tent & GetTent(int i) { return *tents[i];}

  // Return  max(|| gradphi_top||, ||gradphi_bot||)
  double MaxSlope() const;

  // Drawing
  void DrawPitchedTents(int level=1) ;
  void DrawPitchedTentsVTK(string vtkfilename);
  void DrawPitchedTentsGL(Array<int> & tentdata,
                          Array<double> & tenttimes, int & nlevels);


  // Propagate methods need to access this somehow
  Table<int> tent_dependency; // DAG of tent dependencies
};


#endif
