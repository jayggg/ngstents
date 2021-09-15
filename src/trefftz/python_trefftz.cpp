#include <python_ngstd.hpp>
#include <solve.hpp>
#include <fem.hpp>
using namespace ngsolve;
#include "python_fem.hpp"
#include "trefftzfespace.hpp"
#include "specialcoefficientfunction.hpp"
#include <tents.hpp>
#include "twavetents.hpp"

PYBIND11_MODULE(_trefftz,m) {
    py::module::import("ngsolve");
    //py::module::import("ngstents");
    m.attr("__name__") = "ngstents.trefftz";
    m.attr("__package__") = "ngstents";

    ExportTrefftzFESpace(m);
    ExportTWaveTents(m);
    ExportSpecialCoefficientFunction(m);
}

