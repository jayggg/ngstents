#include "tents.hpp"
#include <python_ngstd.hpp>



// Just for initial proof of concept (we need to move from DG -> MTP)
extern void ExportSymbolicDG (py::module & m); 

// Export tent mesh
extern void ExportTents(py::module & m);


PYBIND11_MODULE(tents, m) {

  ExportSymbolicDG(m);
  ExportTents(m);
}

  
