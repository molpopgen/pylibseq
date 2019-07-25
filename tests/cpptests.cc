#include <pybind11/pybind11.h>

namespace py = pybind11;
void init_create_AlleleCountMatrix(py::module&);
void init_create_VariantMatrix(py::module&);

PYBIND11_MODULE(cpptests, m)
{
    init_create_AlleleCountMatrix(m);
    init_create_VariantMatrix(m);
}

