#include <pybind11/pybind11.h>
#include <Sequence/AlleleCountMatrix.hpp>

namespace py = pybind11;

Sequence::AlleleCountMatrix
create_AlleleCountMatrix()
{
    std::vector<std::int32_t> counts;
    return Sequence::AlleleCountMatrix(std::move(counts), 0, 0, 0);
}

void
init_create_AlleleCountMatrix(py::module& m)
{
    m.def("create_AlleleCountMatrix", &create_AlleleCountMatrix);
}
