#include "tconservationlaw_tp_impl.hpp"
#include <python_ngstd.hpp>

extern void ExportBurgers (py::module & m);

void ExportConsLaw(py::module & m)
{
  ExportBurgers(m);
}

// Just for initial proof of concept (we need to move from DG -> MTP)
extern void ExportSymbolicDG (py::module & m); 

PYBIND11_MODULE(conslaw, m) {
  ExportConsLaw(m);
  ExportSymbolicDG(m);
}
