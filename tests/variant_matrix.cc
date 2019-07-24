#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Sequence/VariantMatrix.hpp>

namespace py = pybind11;

Sequence::VariantMatrix
create_VariantMatrix(py::list g, py::list p)
{
    return Sequence::VariantMatrix(g.cast<std::vector<std::int8_t>>(),
                                   p.cast<std::vector<double>>());
}

void
init_create_VariantMatrix(py::module& m)
{
    m.def("create_VariantMatrix", &create_VariantMatrix);
}

