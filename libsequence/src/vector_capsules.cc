#include <pybind11/pybind11.h>
#include <Sequence/VectorCapsules.hpp>

namespace py = pybind11;

void
init_vector_capsules(py::module& m)
{
    py::class_<Sequence::GenotypeCapsule>(m, "__GenotypeCapsule");
    py::class_<Sequence::PositionCapsule>(m, "__PositionCapsule")
        .def("nsites", &Sequence::PositionCapsule::nsites)
        .def("__get__", [](const Sequence::PositionCapsule& self,
                           std::size_t i) { return self[i]; });

    py::class_<Sequence::VectorPositionCapsule, Sequence::PositionCapsule>(
        m, "__VectorPositionCapsule")
        .def(py::init<const std::vector<double>&>())
        .def("nsites", &Sequence::VectorPositionCapsule::nsites)
        .def("__get__", [](const Sequence::VectorPositionCapsule& self,
                           std::size_t i) { return self[i]; });
}
