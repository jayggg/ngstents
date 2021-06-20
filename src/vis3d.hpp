#ifndef VIS3DHEADER
#define VIS3DHEADER

#include <solve.hpp>
#include "tents.hpp"

class Visualization3D {

public:
  Visualization3D(const ngsolve::Array<shared_ptr<ngsolve::Table<int>>> &aidx3d) : idx3d(aidx3d) {
    tmph1 = nullptr;
    vtmp = nullptr;
  }

  // Set the initial data (layer 0) for the timestep/slab
  // into the GridFunction vector for the 3D H1 space.
  void SetInitialHd(shared_ptr<GridFunction> gfu,
                    shared_ptr<GridFunction> hdgf, LocalHeap & lh);

  // Interpolate the solution on elements of a tent into a temp H1 space
  // Then transfer the tent vertex value to the 3D H1 space
  void SetForTent(Tent &tent, shared_ptr<GridFunction> gfu,
                  shared_ptr<GridFunction> hdgf, LocalHeap & lh);

private:
  // idx3d[layer][2D vertex nr] --> vertex nr in 3D mesh
  Array<shared_ptr<Table<int>>> idx3d;
  // temporary H1 space used by SetForTent only for tmph1->TransformVec()
  shared_ptr<FESpace> tmph1;
  // GridFunction vector associated with temp H1 space
  // When used by SetForTent, only vector components associated with tent
  // element dofs are affected; so no race condition is present.
  shared_ptr<BaseVector> vtmp;
};
#endif
