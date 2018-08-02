#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::int8_t>);

PYBIND11_MODULE(variant_matrix, m)
{
    py::bind_vector<std::vector<double>>(m, "VectorDouble",
                                         py::buffer_protocol());
    py::bind_vector<std::vector<std::int8_t>>(m, "VectorInt8",
                                              py::buffer_protocol());

    py::class_<Sequence::VariantMatrix>(m, "VariantMatrix")
        .def(py::init<const std::vector<std::int8_t> &,
                      const std::vector<double> &>(),
             py::arg("data"), py::arg("positions"))
        .def(py::init([](py::list data, py::list pos) {
                 std::vector<std::int8_t> d;
                 std::vector<double> p;
                 for (auto i : data)
                     {
                         d.push_back(i.cast<std::int8_t>());
                     }
                 for (auto i : pos)
                     {
                         p.push_back(i.cast<double>());
                     }
                 return Sequence::VariantMatrix(std::move(d), std::move(p));
             }),
             py::arg("data"), py::arg("pos"))
        .def_readonly("data", &Sequence::VariantMatrix::data)
        .def_readonly("positions", &Sequence::VariantMatrix::positions)
        .def_readonly("nsites", &Sequence::VariantMatrix::nsites)
        .def_readonly("nsam", &Sequence::VariantMatrix::nsam)
        .def_readonly_static("mask", &Sequence::VariantMatrix::mask)
        .def("site",
             [](const Sequence::VariantMatrix &m, const std::size_t i) {
                 return Sequence::get_ConstRowView(m, i);
             })
        .def("sample",
             [](const Sequence::VariantMatrix &m, const std::size_t i) {
                 return Sequence::get_ConstColView(m, i);
             });

    py::class_<Sequence::ConstColView>(m, "ConstColView")
        .def("__len__",
             [](const Sequence::ConstColView &c) { return c.size(); })
        .def("__iter__",
             [](const Sequence::ConstColView &c) {
                 return py::make_iterator(c.begin(), c.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list", [](const Sequence::ConstColView &c) {
            py::list rv;
            for (auto i : c)
                {
                    rv.append(static_cast<int>(i));
                }
            return rv;
        });

    py::class_<Sequence::ConstRowView>(m, "ConstRowView")
        .def("__len__",
             [](const Sequence::ConstRowView &r) { return r.size(); })
        .def("__iter__",
             [](const Sequence::ConstRowView &r) {
                 return py::make_iterator(r.begin(), r.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list", [](const Sequence::ConstRowView &r) {
            py::list rv;
            for (auto i : r)
                {
                    rv.append(static_cast<int>(i));
                }
            return rv;
        });
}
