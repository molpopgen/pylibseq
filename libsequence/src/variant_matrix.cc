#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Sequence/AlleleCountMatrix.hpp>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/filtering.hpp>
#include <Sequence/variant_matrix/windows.hpp>
#include <Sequence/variant_matrix/msformat.hpp>
#include <Sequence/StateCounts.hpp>

namespace py = pybind11;

class NumpyGenotypeCapsule : public Sequence::GenotypeCapsule
{
  private:
    py::array_t<std::int8_t> buffer;
    std::size_t nsites_, nsam_;

  public:
    explicit NumpyGenotypeCapsule(
        py::array_t<std::int8_t, py::array::c_style | py::array::forcecast>
            input)
        : buffer(std::move(input)), nsites_(buffer.shape(0)),
          nsam_(buffer.shape(1))
    {
    }

    std::size_t &
    nsites()
    {
        return nsites_;
    }

    std::size_t &
    nsam()
    {
        return nsam_;
    }

    std::size_t
    nsites() const
    {
        return nsites_;
    }

    std::size_t
    nsam() const
    {
        return nsam_;
    }

    std::int8_t &operator[](std::size_t i) { return buffer.mutable_data()[i]; }

    const std::int8_t &operator[](std::size_t i) const
    {
        return buffer.data()[i];
    }

    std::int8_t *
    data() final
    {
        return buffer.mutable_data();
    }

    const std::int8_t *
    data() const final
    {
        return buffer.data();
    }

    const std::int8_t *
    cdata() const final
    {
        return buffer.data();
    }

    std::unique_ptr<GenotypeCapsule>
    clone() const final
    {
        return std::unique_ptr<GenotypeCapsule>(
            new NumpyGenotypeCapsule(this->buffer));
    }

    std::int8_t *
    begin() final
    {
        return buffer.mutable_data();
    }

    const std::int8_t *
    begin() const final
    {
        return buffer.data();
    }

    std::int8_t *
    end() final
    {
        return buffer.mutable_data() + buffer.size();
    }

    const std::int8_t *
    end() const final
    {
        return buffer.data() + buffer.size();
    }

    const std::int8_t *
    cbegin() const final
    {
        return begin();
    }

    const std::int8_t *
    cend() const final
    {
        return end();
    }

    bool
    empty() const final
    {
        return !nsites();
    }

    std::size_t
    size() const final
    {
        return buffer.size();
    }

    std::size_t
    row_offset() const
    {
        return 0;
    }

    std::size_t
    col_offset() const
    {
        return 0;
    }

    std::size_t
    stride() const
    {
        return nsam();
    }

    std::int8_t &
    operator()(std::size_t site, std::size_t sample)  final
    {
        return buffer.mutable_data()[nsam()*site + sample];
    }

    const std::int8_t &
    operator()(std::size_t site, std::size_t sample) const final
    {
        return buffer.data()[nsam()*site + sample];
    }

    bool
    resizable() const final
    {
        return false;
    }
};

class NumpyPositionCapsule : public Sequence::PositionCapsule
{
  private:
    py::array_t<double> buffer;

  public:
    explicit NumpyPositionCapsule(py::array_t<double> input)
        : buffer(std::move(input))
    {
        if (buffer.ndim() != 1)
            {
                throw std::invalid_argument(
                    "positions must be a one-dimensional array");
            }
    }

    double &operator[](std::size_t i) { return buffer.mutable_data()[i]; }

    const double &operator[](std::size_t i) const { return buffer.data()[i]; }

    double *
    data() final
    {
        return buffer.mutable_data();
    }

    const double *
    data() const final
    {
        return buffer.data();
    }

    const double *
    cdata() const final
    {
        return buffer.data();
    }

    std::unique_ptr<PositionCapsule>
    clone() const final
    {
        return std::unique_ptr<PositionCapsule>(
            new NumpyPositionCapsule(*this));
    }

    double *
    begin() final
    {
        return buffer.mutable_data();
    }

    const double *
    begin() const final
    {
        return buffer.data();
    }

    const double *
    cbegin() const final
    {
        return buffer.data();
    }

    double *
    end() final
    {
        return buffer.mutable_data() + buffer.size();
    }

    const double *
    end() const final
    {
        return buffer.data() + buffer.size();
    }

    const double *
    cend() const final
    {
        return buffer.data() + buffer.size();
    }

    bool
    empty() const final
    {
        return buffer.size() == 0;
    }

    std::size_t
    size() const final
    {
        return buffer.size();
    }

    std::size_t
    nsites() const
    {
        return buffer.size();
    }

    bool
    resizable() const final
    {
        return false;
    }
};

class MockVM
{
  private:
    std::unique_ptr<Sequence::GenotypeCapsule> g;
    std::unique_ptr<Sequence::PositionCapsule> p;

  public:
    MockVM(std::unique_ptr<Sequence::GenotypeCapsule> g_,
           std::unique_ptr<Sequence::PositionCapsule> p_)
        : g(std::move(g_)), p(std::move(p_))
    {
    }
};

Sequence::VariantMatrix
mslike_from_numpy(
    py::array_t<std::int8_t, py::array::c_style | py::array::forcecast>
        genotypes,
    py::array_t<double> positions)
{
    std::unique_ptr<Sequence::GenotypeCapsule> pp(
        new NumpyGenotypeCapsule(std::move(genotypes)));
    std::unique_ptr<Sequence::PositionCapsule> gp(
        new NumpyPositionCapsule(std::move(positions)));
    return Sequence::VariantMatrix(std::move(pp), std::move(gp), 1);
}

void
init_VariantMatrix(py::module &m)
{
    py::class_<Sequence::AlleleCountMatrix>(
        m, "AlleleCountMatrix", py::buffer_protocol(),
        "A matrix of allele counts. This object supports the buffer "
        "protocol.")
        .def(py::init<const Sequence::VariantMatrix &>(),
             "Construct from a "
             ":class:`libsequence.variant_matrix.VariantMatrix`")
        .def_readonly("counts", &Sequence::AlleleCountMatrix::counts,
                      "Flattened view of the raw data.")
        .def_readonly("nrow", &Sequence::AlleleCountMatrix::nrow,
                      "Number of rows (sites) in the matrix.")
        .def_readonly("ncol", &Sequence::AlleleCountMatrix::ncol,
                      "Number of columns (allelic states) in the matrix.")
        .def_readonly("nsam", &Sequence::AlleleCountMatrix::nsam,
                      "Sample size of the original VariantMatrix.")
        .def(
            "row",
            [](const Sequence::AlleleCountMatrix &c, const std::size_t i) {
                auto x = c.row(i);
                return py::make_iterator(x.first, x.second);
            },
            py::keep_alive<0, 1>(), py::arg("i"),
            "Return an iterator over the i-th site.")
        .def("__getitem__",
             [](const Sequence::AlleleCountMatrix &am, py::slice slice) {
                 std::size_t start, stop, step, slicelength;
                 if (!slice.compute(am.counts.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();
                 std::vector<int> c;
                 int nrow = 0;
                 for (size_t i = 0; i < slicelength; ++i, ++nrow)
                     {
                         auto r = am.row(start);
                         c.insert(c.end(), r.first, r.second);
                         start += step;
                     }
                 return Sequence::AlleleCountMatrix(std::move(c), am.ncol,
                                                    nrow, am.nsam);
             })
        .def("__getitem__",
             [](const Sequence::AlleleCountMatrix &am,
                py::array_t<std::size_t> x) {
                 auto r = x.unchecked<1>();
                 std::vector<int> c;
                 int nrow = 0;
                 for (std::size_t i = 0; i < r.shape(0); ++i, ++nrow)
                     {
                         auto row = am.row(r(i));
                         c.insert(c.end(), row.first, row.second);
                     }
                 return Sequence::AlleleCountMatrix(std::move(c), am.ncol,
                                                    nrow, am.nsam);
             })
        .def("__len__",
             [](const Sequence::AlleleCountMatrix &self) { return self.nrow; })
        .def_buffer(
            [](const Sequence::AlleleCountMatrix &c) -> py::buffer_info {
                using value_type = Sequence::AlleleCountMatrix::value_type;
                return py::buffer_info(
                    //TODO: fix this const_cast hack
                    //once pybind11 supports read-only
                    //buffers!
                    const_cast<value_type *>(
                        c.counts.data()), /* Pointer to buffer */
                    sizeof(value_type),   /* Size of one scalar */
                    py::format_descriptor<value_type>::
                        format(), /* Python struct-style format descriptor */
                    2,            /* Number of dimensions */
                    { c.nrow, c.ncol }, /* Buffer dimensions */
                    {
                        sizeof(value_type) * c.ncol, sizeof(value_type)
                        /* Strides (in bytes) for each index */
                    });
            });

    py::class_<Sequence::VariantMatrix>(m, "VariantMatrix",
                                        //py::buffer_protocol(),
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
        :param positions: List of mutation positions.
        :type positions: list

        >>> import libsequence
        >>> m = libsequence.VariantMatrix([0,1,1,0],[0.1,0.2])
             )delim",
             py::arg("data"), py::arg("positions"))
        .def(py::init([](py::array_t<std::int8_t,
                                     py::array::c_style | py::array::forcecast>
                             data,
                         py::array_t<double> pos,
                         std::int8_t max_allele_value) {
                 std::unique_ptr<Sequence::GenotypeCapsule> dp(
                     new NumpyGenotypeCapsule(std::move(data)));
                 std::unique_ptr<Sequence::PositionCapsule> pp(
                     new NumpyPositionCapsule(std::move(pos)));
                 return Sequence::VariantMatrix(std::move(dp), std::move(pp),
                                                max_allele_value);
             }),
             R"delim(
             Construct with numpy arrays

            :param data: 2d ndarray with dtype numpy.int8
            :type data: list
            :param positions: 1d array with dtype np.float
            :type positions: list

            >>> import libsequence
            >>> import numpy as np
            >>> d = np.array([0,1,1,0],dtype=np.int8).reshape((2,2))
            >>> p = np.array([0.1,0.2])
            >>> m = libsequence.VariantMatrix(d,p)
            )delim",
             py::arg("data"), py::arg("pos"), py::arg("max_allele_value") = -1)
        .def_static(
            "from_TreeSequence",
            [](py::object ts) -> Sequence::VariantMatrix {
                //If a not-TreeSequence is passed in, "duck typing"
                //fails, and an exception will be raised.
                auto g = ts.attr("genotype_matrix")();
                auto p = ts.attr("tables")
                             .attr("sites")
                             .attr("position")
                             .cast<py::array_t<double>>();
                return mslike_from_numpy(std::move(g), std::move(p));
            },
            py::arg("ts"),
            R"delim(
            Create a VariantMatrix from an msprime.TreeSequence
            
            :param ts: A TreeSequence from msprime :cite:`Kelleher2016-cb`.
            
            A TreeSequence object is the output of `msprime.simulate`,
            or, equivalently, certain forward simulations that use
            that format for storing results.

            This function is a convenience function. Internally,
            the output from msprime are cast from 8-bit unsigned
            integers to 8-bit signed integers.
            )delim")
        .def_property_readonly(
            "data",
            [](const Sequence::VariantMatrix &self) {
                auto rv = pybind11::array_t<std::int8_t>(
                    { self.nsites(), self.nsam() }, self.cdata(),
                    pybind11::cast(self));
                rv.attr("flags").attr("writeable") = false;
                return rv;
            },
            "Return raw data as numpy array")
        .def_property_readonly(
            "positions",
            [](const Sequence::VariantMatrix &self)
                -> pybind11::array_t<double> {
                auto rv = pybind11::array_t<double>(
                    { self.nsites() }, { sizeof(double) }, self.pbegin(),
                    pybind11::cast(self.pbegin()));
                rv.attr("flags").attr("writeable") = false;
                return rv;
            },
            "Return positions as numpy array")
        .def_property_readonly(
            "nsites",
            [](const Sequence::VariantMatrix &self) { return self.nsites(); },
            "Number of positions")
        .def_property_readonly(
            "nsam",
            [](const Sequence::VariantMatrix &self) { return self.nsam(); },
            "Number of samples")
        .def_readonly_static("mask", &Sequence::VariantMatrix::mask,
                             "Reserved missing data state")
        .def("count_alleles",
             [](const Sequence::VariantMatrix &m) {
                 return Sequence::AlleleCountMatrix(m);
             })
        .def(
            "site",
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
        .def(
            "sample",
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
        //.def_buffer([](Sequence::VariantMatrix &m) -> py::buffer_info {
        //    return py::buffer_info(
        //        m.data(),            /* Pointer to buffer */
        //        sizeof(std::int8_t), /* Size of one scalar */
        //        py::format_descriptor<std::int8_t>::
        //            format(), /* Python struct-style format descriptor */
        //        2,            /* Number of dimensions */
        //        { m.nsites(), m.nsam() }, /* Buffer dimensions */
        //        { sizeof(std::int8_t)
        //              * m.nsam(), /* Strides (in bytes) for each index */
        //          sizeof(std::int8_t) });
        //})
        .def(
            "window",
            [](const Sequence::VariantMatrix &m, const double beg,
               const double end) {
                return Sequence::make_window(m, beg, end);
            },
            py::arg("beg"), py::arg("end"))
        .def(
            "slice",
            [](const Sequence::VariantMatrix &m, const double beg,
               const double end, const std::size_t i, const std::size_t j) {
                return Sequence::make_slice(m, beg, end, i, j);
            },
            py::arg("beg"), py::arg("end"), py::arg("i"), py::arg("j"))
        .def(py::pickle(
            [](const Sequence::VariantMatrix &m) {
                std::vector<std::int8_t> temp(
                    m.data(), m.data() + m.nsites() * m.nsam());
                std::vector<double> ptemp(m.pbegin(), m.pend());
                return py::make_tuple(std::move(temp), std::move(ptemp));
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
        .def(
            "__iter__",
            [](const Sequence::ConstColView &c) {
                return py::make_iterator(c.begin(), c.end());
            },
            py::keep_alive<0, 1>())
        .def(
            "as_list",
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
        .def(
            "__iter__",
            [](const Sequence::ColView &c) {
                return py::make_iterator(c.begin(), c.end());
            },
            py::keep_alive<0, 1>())
        .def(
            "as_list",
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
        .def(
            "__iter__",
            [](const Sequence::ConstRowView &r) {
                return py::make_iterator(r.begin(), r.end());
            },
            py::keep_alive<0, 1>())
        .def(
            "as_list",
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
        .def(
            "__iter__",
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

    py::class_<Sequence::StateCounts>(m, "StateCounts", py::buffer_protocol(),
                                      R"delim(
            Count the states at a site in a VariantMatrix.

            See :ref:`variantmatrix`
            )delim")
        .def(py::init<>())
        .def(py::init<std::int8_t>(), py::arg("refstate"))
        .def_readonly("counts", &Sequence::StateCounts::counts,
                      "The counts for each possible non-missing allelic state")
        .def_readonly("refstate", &Sequence::StateCounts::refstate,
                      "The reference state.")
        .def_readonly("n", &Sequence::StateCounts::n, "The sample size.")
        .def(
            "__iter__",
            [](const Sequence::StateCounts &sc) {
                return py::make_iterator(sc.counts.begin(), sc.counts.end());
            },
            py::keep_alive<0, 1>())
        .def("__len__",
             [](const Sequence::StateCounts &c) { return c.counts.size(); })
        .def("__getitem__",
             [](const Sequence::StateCounts &c, const std::size_t i) {
                 if (i >= c.counts.size())
                     {
                         throw std::invalid_argument("index out of range");
                     }
                 return c.counts[i];
             })
        .def("__call__",
             [](Sequence::StateCounts &c, Sequence::ConstRowView &r) { c(r); })
        .def("__call__", [](Sequence::StateCounts &c,
                            const Sequence::RowView &r) { c(r); })
        .def_buffer([](Sequence::StateCounts &c) -> py::buffer_info {
            return py::buffer_info(
                c.counts.data(),      /* Pointer to buffer */
                sizeof(std::int32_t), /* Size of one scalar */
                py::format_descriptor<std::int32_t>::
                    format(), /* Python struct-style format descriptor */
                1,            /* Number of dimensions */
                { c.counts.size() }, /* Buffer dimensions */
                {
                    sizeof(std::int32_t)
                    /* Strides (in bytes) for each index */
                });
        });

    m.def(
        "process_variable_sites",
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
        py::arg("m"), py::arg("refstates") = nullptr,
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

    m.def(
        "filter_haplotypes",
        [](Sequence::VariantMatrix &m, py::function f) {
            if (m.resizable())
                {
                    auto cpp_func = f.cast<
                        std::function<bool(const Sequence::ColView &)>>();
                    return Sequence::filter_haplotypes(m, cpp_func);
                }
            auto cpp_func = f.cast<
                std::function<bool(const Sequence::ConstColView &)>>();
            return Sequence::filter_haplotypes(m, cpp_func);
        },
        R"delim(
            Remove site data from a VariantMatrix

            :param m: A variant matrix
            :type m: :class:`libsequence.variant_matrix.VariantMatrix`
            :param f: A function
            :type f: callable

            See :ref:`variantmatrix` for details.

            .. note::

                If the the VariantMatrix was initially
                created from a numpy array, its internal state
                is changed to be based on C++ vectors, meaning
                that a copy happened internally.
            )delim",
        py::arg("m"), py::arg("f"));

    m.def(
        "filter_sites",
        //[](Sequence::VariantMatrix &m,
        //   const std::function<bool(const Sequence::ConstRowView &)> &f) {
        //    return Sequence::filter_sites(m, f);
        //},
        [](Sequence::VariantMatrix &m, py::function f) {
            if (m.resizable())
                {
                    auto cpp_func = f.cast<
                        std::function<bool(const Sequence::RowView &)>>();

                    return Sequence::filter_sites(m, cpp_func);
                }
            auto cpp_func = f.cast<
                std::function<bool(const Sequence::ConstRowView &)>>();

            return Sequence::filter_sites(m, cpp_func);
        },
        R"delim(
            Remove sample data from a VariantMatrix

            :param m: A variant matrix
            :type m: :class:`libsequence.variant_matrix.VariantMatrix`
            :param f: A function
            :type f: callable

            See :ref:`variantmatrix` for details.

            .. note::

                If the the VariantMatrix was initially
                created from a numpy array, its internal state
                is changed to be based on C++ vectors, meaning
                that a copy happened internally.
            )delim",
        py::arg("m"), py::arg("f"));

    // Various I/O for VariantMatrix in "ms" format
    m.def("ms_from_stdin", []() -> py::object {
        if (std::cin.eof())
            {
                return py::none();
            }
        auto vm = Sequence::from_msformat(std::cin);
        py::object o = py::cast(std::move(vm));
        return o;
    });

    py::class_<MockVM>(m, "MockVM")
        .def(py::init([](py::array_t<std::int8_t,
                                     py::array::c_style | py::array::forcecast>
                             data,
                         py::array_t<double> pos,
                         std::int8_t max_allele_value) {
            std::unique_ptr<Sequence::GenotypeCapsule> dp(
                new NumpyGenotypeCapsule(data));
            std::unique_ptr<Sequence::PositionCapsule> pp(
                new NumpyPositionCapsule(pos));
            return MockVM(std::move(dp), std::move(pp));
        }));
}
