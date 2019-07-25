#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Sequence/VariantMatrix.hpp>

namespace py = pybind11;
using genotype_vector = std::vector<std::int8_t>;
using position_vector = std::vector<double>;

py::object
create_VariantMatrix(py::list g, py::list p)
{
    auto __l = py::module::import("libsequence");
    return __l.attr("create_VariantMatrix")(g, p);
}

py::object
create_VariantMatrix_from_cpp_vectors()
// This is not the only way to go here.
// You can imagine fancier solutions
// where the capsules DO NOT delete
// the pointers.  For example, if the
// vectors are themselves held in classes
// and all object lifetimes are managed
// on the Python side.
{
    // We need pointers to our C++ containers
    genotype_vector* g = new genotype_vector{ 0, 1, 1, 0 };
    position_vector* p = new position_vector{ 0.1, 0.2 };
    // These pointers then go into PyCapsule objects,
    // which are basically smart pointers with custom
    // deleters
    auto gc = py::capsule(
        g, [](void* v) { delete reinterpret_cast<genotype_vector*>(v); });
    auto pc = py::capsule(
        p, [](void* v) { delete reinterpret_cast<position_vector*>(v); });

    //These capsules can be the memory buffers of numpy arrays...
    py::array_t<std::int8_t> ga(g->size(), g->data(), gc);
    py::array_t<double> pa(p->size(), p->data(), pc);
    //...from which we may fill libsequence's VariantMatrix
    //from a new Python module that we've written
    auto __l = py::module::import("libsequence");
    return __l.attr("VariantMatrix")(ga, pa);
}

void
init_create_VariantMatrix(py::module& m)
{
    m.def("create_VariantMatrix", &create_VariantMatrix);
    m.def("create_VariantMatrix_from_cpp_vectors",
          &create_VariantMatrix_from_cpp_vectors);
}

