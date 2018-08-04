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

PYBIND11_MODULE(variant_matrix, m)
{
    py::class_<Sequence::VariantMatrix>(m, "VariantMatrix",
                                        py::buffer_protocol(),
                                        R"delim(
        Representation of variation data in matrix format.

        see :ref:`variantmatrix` for discussion.
        )delim")
        .def(py::init<const std::vector<std::int8_t> &,
                      const std::vector<double> &>(),
             R"delim(
        Construct with a lists of input data.

        :param data: The state data.
        :type data: list
        :param positons: List of mutation positions.
        :type positions: list

        >>> import libsequence.variant_matrix as vm
        >>> m = vm.VariantMatrix([0,1,1,0],[0.1,0.2])
             )delim",
             py::arg("data"), py::arg("positions"))
        .def(py::init([](py::array_t<std::int8_t,
                                     py::array::c_style | py::array::forcecast>
                             data,
                         py::array_t<double> pos) {
                 if (data.ndim() != 2)
                     {
                         throw std::invalid_argument(
                             "data must be a 2d ndarray");
                     }
                 if (pos.ndim() != 1)
                     {
                         throw std::invalid_argument(
                             "pos must be a 1d ndarray");
                     }
                 if (pos.size() != data.shape(0))
                     {
                         throw(std::invalid_argument(
                             "len(pos) must equal data.shape[0]"));
                     }
                 auto data_ptr = data.unchecked<2>();
                 std::vector<std::int8_t> d;
                 d.reserve(data.size());
                 for (decltype(data_ptr.shape(0)) i = 0; i < data_ptr.shape(0);
                      ++i)
                     {
                         for (decltype(data_ptr.shape(1)) j = 0;
                              j < data_ptr.shape(1); ++j)
                             {
                                 d.push_back(*data_ptr.data(i, j));
                             }
                     }
                 auto pos_ptr = pos.unchecked<1>();
                 std::vector<double> p(pos_ptr.data(0),
                                       pos_ptr.data(0) + pos.size());
                 return Sequence::VariantMatrix(std::move(d), std::move(p));
             }),
             R"delim(
             Construct with numpy arrays

            :param data: 2d ndarray with dtype numpy.int8
            :type data: list
            :param positons: 1d array with dtype np.float
            :type positions: list

            >>> import libsequence.variant_matrix as vm
            >>> import numpy as np
            >>> d = np.array([0,1,1,0],dtype=np.int8).reshape((2,2))
            >>> p = np.array([0.1,0.2])
            >>> m = vm.VariantMatrix(d,p)
            )delim",
             py::arg("data"), py::arg("pos"))
        .def_readonly("data", &Sequence::VariantMatrix::data,
                      "Return raw data as list")
        .def_readonly("positions", &Sequence::VariantMatrix::positions,
                      "Return positions as list")
        .def_readonly("nsites", &Sequence::VariantMatrix::nsites,
                      "Number of positions")
        .def_readonly("nsam", &Sequence::VariantMatrix::nsam,
                      "Number of samples")
        .def_readonly_static("mask", &Sequence::VariantMatrix::mask,
                             "Reserved missing data state")
        .def("site",
             [](const Sequence::VariantMatrix &m, const std::size_t i) {
                 return Sequence::get_ConstRowView(m, i);
             },
             R"delim(
             Return a view of the i-th site.
             
             :param i: Index
             :type i: int
             :rtype: :class:`libsequence.variant_matrix.ConstRowView`
             )delim",
             py::arg("i"))
        .def("sample",
             [](const Sequence::VariantMatrix &m, const std::size_t i) {
                 return Sequence::get_ConstColView(m, i);
             },
             R"delim(
             Return a view of the i-th sample.
             
             :param i: Index
             :type i: int
             :rtype: :class:`libsequence.variantmatrix.ConstColView`
             )delim",
             py::arg("i"))
        .def_buffer([](Sequence::VariantMatrix &m) -> py::buffer_info {
            return py::buffer_info(
                m.data.data(),       /* Pointer to buffer */
                sizeof(std::int8_t), /* Size of one scalar */
                py::format_descriptor<std::int8_t>::
                    format(), /* Python struct-style format descriptor */
                2,            /* Number of dimensions */
                { m.nsites, m.nsam }, /* Buffer dimensions */
                { sizeof(std::int8_t)
                      * m.nsam, /* Strides (in bytes) for each index */
                  sizeof(std::int8_t) });
        })
        .def(py::pickle(
            [](const Sequence::VariantMatrix &m) {
                return py::make_tuple(m.data, m.positions);
            },
            [](py::tuple t) {
                if (t.size() != 2)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto d = t[0].cast<std::vector<std::int8_t>>();
                auto p = t[1].cast<std::vector<double>>();
                return Sequence::VariantMatrix(std::move(d), std::move(p));
            }));

    py::class_<Sequence::ConstColView>(m, "ConstColView",
                                       R"delim(
            Immutable view of a VariantMatrix column.

            See :ref:`variantmatrix`.
            )delim")
        .def("__len__",
             [](const Sequence::ConstColView &c) { return c.size(); })
        .def("__iter__",
             [](const Sequence::ConstColView &c) {
                 return py::make_iterator(c.begin(), c.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list",
             [](const Sequence::ConstColView &c) {
                 py::list rv;
                 for (auto i : c)
                     {
                         rv.append(static_cast<int>(i));
                     }
                 return rv;
             },
             "Return contents as a list.");

    py::class_<Sequence::ColView>(m, "ColView",
                                  R"delim(
        View of a VariantMatrix column.

        See :ref:`variantmatrix`
        )delim")
        .def("__len__", [](const Sequence::ColView &c) { return c.size(); })
        .def("__iter__",
             [](const Sequence::ColView &c) {
                 return py::make_iterator(c.begin(), c.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list",
             [](const Sequence::ColView &c) {
                 py::list rv;
                 for (auto i : c)
                     {
                         rv.append(static_cast<int>(i));
                     }
                 return rv;
             },
             "Return contents as a list.");

    py::class_<Sequence::ConstRowView>(m, "ConstRowView",
                                       R"delim(
        Immutable view of a sample.

        See :ref:`variantmatrix`.
        )delim")
        .def("__len__",
             [](const Sequence::ConstRowView &r) { return r.size(); })
        .def("__iter__",
             [](const Sequence::ConstRowView &r) {
                 return py::make_iterator(r.begin(), r.end());
             },
             py::keep_alive<0, 1>())
        .def("as_list",
             [](const Sequence::ConstRowView &r) {
                 py::list rv;
                 for (auto i : r)
                     {
                         rv.append(static_cast<int>(i));
                     }
                 return rv;
             },
             "Return contents as a list.");

    py::class_<Sequence::RowView>(m, "RowView",
                                  R"delim(
        View of a sample in a VariantMatrix.

        See :ref:`variantmatrix`.
        )delim")
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

    py::class_<Sequence::StateCounts>(m, "StateCounts",
                                      R"delim(
            Count the states at a site in a VariantMatrix.

            See :ref:`variantmatrix`
            )delim")
        .def(py::init<>())
        .def(py::init<std::int8_t>())
        .def_readonly("counts", &Sequence::StateCounts::counts)
        .def_readonly("refstate", &Sequence::StateCounts::refstate)
        .def_readonly("n", &Sequence::StateCounts::n)
        .def("__iter__",
             [](const Sequence::StateCounts &sc) {
                 return py::make_iterator(sc.counts.begin(), sc.counts.end());
             },
             py::keep_alive<0, 1>());

    m.def("process_variable_sites",
          [](const Sequence::VariantMatrix &m, py::object refstates) {
              if (refstates.is_none())
                  {
                      return Sequence::process_variable_sites(m);
                  }
              try
                  {
                      py::int_ rs(refstates);
                      return Sequence::process_variable_sites(
                          m, rs.cast<std::int8_t>());
                  }
              catch (...)
                  {
                  }
              return Sequence::process_variable_sites(
                  m, refstates.cast<std::vector<std::int8_t>>());
          },
          py::arg("m"), py::arg("refstates")=nullptr,
          R"delim(
          Obtain state counts for all sites

          :param m: data
          :type m: :class:`libsequence.variant_matrix.VariantMatrix`
          :param refstates: The reference states for each sites.
          :type refstates: object
          :return: list of :class:`libsequence.variant_matrix.StateCounts`
          :type: list

          See :ref:`variantmatrix` for examples.
          )delim");

    //m.def("process_variable_sites",
    //      [](const Sequence::VariantMatrix &m, const std::int8_t refstate) {
    //          return Sequence::process_variable_sites(m, refstate);
    //      },
    //      py::arg("m"), py::arg("refstate"));

    //m.def("process_variable_sites",
    //      [](const Sequence::VariantMatrix &m) {
    //          return Sequence::process_variable_sites(m);
    //      },
    //      py::arg("m"));

    m.def("filter_haplotypes", &Sequence::filter_haplotypes,
          R"delim(
            Remove site data from a VariantMatrix

            :param m: A variant matrix
            :type m: :class:`libsequence.variant_matrix.VariantMatrix`
            :param f: A function
            :type f: callable

            See :ref:`variantmatrix` for details.
            )delim",
          py::arg("m"), py::arg("f"));

    m.def("filter_sites", &Sequence::filter_sites,
          R"delim(
            Remove sample data from a VariantMatrix

            :param m: A variant matrix
            :type m: :class:`libsequence.variant_matrix.VariantMatrix`
            :param f: A function
            :type f: callable

            See :ref:`variantmatrix` for details.
            )delim",
          py::arg("m"), py::arg("f"));
}
