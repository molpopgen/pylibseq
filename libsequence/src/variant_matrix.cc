#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/filtering.hpp>
#include <Sequence/StateCounts.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::int8_t>);

PYBIND11_MODULE(variant_matrix, m)
{
    py::bind_vector<std::vector<double>>(m, "VectorDouble",
                                         py::buffer_protocol());
    py::bind_vector<std::vector<std::int8_t>>(m, "VectorInt8",
                                              py::buffer_protocol());

    py::class_<Sequence::VariantMatrix>(m, "VariantMatrix",py::buffer_protocol())
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
        .def(py::init([](py::array_t<std::int8_t> data,
                         py::array_t<double> pos) {
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
             })
        .def_buffer([](Sequence::VariantMatrix &m) -> py::buffer_info {
            return py::buffer_info(
                m.data.data(), /* Pointer to buffer */
                sizeof(std::int8_t), /* Size of one scalar */
                py::format_descriptor<std::int8_t>::
                    format(), /* Python struct-style format descriptor */
                2,            /* Number of dimensions */
                { m.nsites, m.nsam }, /* Buffer dimensions */
                { sizeof(std::int8_t)
                      * m.nsites, /* Strides (in bytes) for each index */
                  sizeof(std::int8_t) });
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

    py::class_<Sequence::ColView>(m, "ColView")
        .def("__len__", [](const Sequence::ColView &c) { return c.size(); })
        .def("__iter__",
             [](const Sequence::ColView &c) {
                 return py::make_iterator(c.begin(), c.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list", [](const Sequence::ColView &c) {
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

    py::class_<Sequence::RowView>(m, "RowView")
        .def("__len__", [](const Sequence::RowView &r) { return r.size(); })
        .def("__iter__",
             [](const Sequence::RowView &r) {
                 return py::make_iterator(r.begin(), r.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list", [](const Sequence::RowView &r) {
            py::list rv;
            for (auto i : r)
                {
                    rv.append(static_cast<int>(i));
                }
            return rv;
        });

    py::class_<Sequence::StateCounts>(m, "StateCounts")
        .def(py::init<const Sequence::ConstRowView &, const std::int8_t>(),
             py::arg("site"), py::arg("refstate"))
        .def(py::init<const Sequence::ConstRowView &>())
        .def_readonly("counts", &Sequence::StateCounts::counts)
        .def_readonly("refstate", &Sequence::StateCounts::counts)
        .def_readonly("n", &Sequence::StateCounts::n)
        .def("__iter__",
             [](const Sequence::StateCounts &sc) {
                 return py::make_iterator(sc.counts.begin(), sc.counts.end());
             },
             py::keep_alive<0, 1>());

    m.def("process_variable_sites",
          [](const Sequence::VariantMatrix &m,
             const std::vector<std::int8_t> &refstates) {
              return Sequence::process_variable_sites(m, refstates);
          },
          py::arg("m"), py::arg("refstates"));
    m.def("process_variable_sites",
          [](const Sequence::VariantMatrix &m, const std::int8_t refstate) {
              return Sequence::process_variable_sites(m, refstate);
          },
          py::arg("m"), py::arg("refstate"));
    m.def("process_variable_sites",
          [](const Sequence::VariantMatrix &m) {
              return Sequence::process_variable_sites(m);
          },
          py::arg("m"));

    m.def("filter_haplotypes", &Sequence::filter_haplotypes);

    m.def("filter_sites", &Sequence::filter_sites);
}
