#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Sequence/VariantMatrix.hpp>

namespace py = pybind11;

py::object
create_VariantMatrix(py::list g, py::list p)
{
    auto __l = py::module::import("libsequence");
    return __l.attr("create_VariantMatrix")(g,p);
}

void
init_create_VariantMatrix(py::module& m)
{
    m.def("create_VariantMatrix", &create_VariantMatrix);
}

