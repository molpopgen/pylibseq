//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of pylibseq.
//
// pylibseq is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pylibseq is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pylibseq.  If not, see <http://www.gnu.org/licenses/>.
//

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Sequence/PolyTable.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/polySiteVector.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/stateCounter.hpp>

namespace py = pybind11;

// PYBIND11_MAKE_OPAQUE(Sequence::polySiteVector)

Sequence::polySiteVector
polySiteVector_from_list(py::list l)
{
    Sequence::polySiteVector temp;
    for (auto element : l)
        {
            py::tuple t = element.cast<py::tuple>();
            temp.emplace_back(t[0].cast<double>(), py::str(t[1]));
        }
    return temp;
}

PYBIND11_PLUGIN(polytable)
{
    py::module m("polytable", "Access to libsequence's polymorphism table "
                              "classes and related functions");

    // py::bind_vector<Sequence::polySiteVector>(m,"OpaquePolySiteVector");

    py::class_<Sequence::PolyTable>(m, "PolyTable",
                                    "Base class for polymorphism tables")
        .def("data", &Sequence::PolyTable::GetData)
        .def("pos", &Sequence::PolyTable::GetPositions)
        .def("assign",
             [](Sequence::PolyTable& p, const Sequence::polySiteVector& v) {
                 auto rv = p.assign(v.cbegin(), v.cend());
                 if (!rv)
                     throw std::runtime_error("assignment failure");
             },
             "Assign data from a tuple of (position,data), where data are "
             "encoded as an str.")
        .def("assign",
             [](Sequence::PolyTable& p, const std::vector<double>& pos,
                const std::vector<std::string>& d) {
                 auto rv = p.assign(pos, d);
                 if (!rv)
                     throw std::runtime_error("assignment failure");
             },
             "Assign data from a list of positions and a  list of strings "
             "representing variable positions.")
        .def("empty", &Sequence::PolyTable::empty,
             "Return whether or not container is empty.")
        .def("numsites", &Sequence::PolyTable::numsites,
             "Return number of variable sites (columns).")
        .def("size", &Sequence::PolyTable::size,
             "Return number of elements (rows).")
        .def("position", &Sequence::PolyTable::position,
             "Return value of the i-th position.", py::arg("i"))
        .def("__getitem__",
             [](const Sequence::PolyTable& pt, const std::size_t i) {
                 if (i >= pt.size())
                     {
                         throw py::index_error("index out of range");
                     }
                 return pt[i];
             },
             "Return i-th row.")
        .def("__len__", [](const Sequence::PolyTable& p) { return p.size(); },
             "Return number of elements (rows).")
        .def("__getstate__", [](const Sequence::PolyTable& p) {
            return py::make_tuple(p.GetPositions(), p.GetData());
        });

    py::class_<Sequence::PolySites, Sequence::PolyTable>(
        m, "PolySites", "A polymorphism table for Sequence data.  "
                        "A/G/C/T/N/0/1 is the alphabet for analysis.")
        .def(py::init<>())
        .def("__init__",
             [](Sequence::PolySites& p, const Sequence::polySiteVector& v) {
                 new (&p) Sequence::PolySites(v.cbegin(), v.cend());
             })
        .def("__setstate__", [](Sequence::PolySites& p, py::tuple t) {
            if (t.size() != 2)
                {
                    throw std::runtime_error("incorrect tuple length");
                }
            new (&p)
                Sequence::PolySites(t[0].cast<std::vector<double>>(),
                                    t[1].cast<std::vector<std::string>>());
        });

    py::class_<Sequence::SimData, Sequence::PolyTable>(
        m, "SimData", "A polymorphism table for "
                      "binary data.  0/1 = "
                      "ancestral/derived.")
        .def(py::init<std::vector<double>, std::vector<std::string>>())
        .def(py::init<>())
        .def("__init__",
             [](Sequence::SimData& d, const Sequence::polySiteVector& p) {
                 new (&d) Sequence::SimData(p.cbegin(), p.cend());
             })
        .def("__setstate__", [](Sequence::SimData& p, py::tuple t) {
            if (t.size() != 2)
                {
                    throw std::runtime_error("incorrect tuple length");
                }
            new (&p) Sequence::SimData(t[0].cast<std::vector<double>>(),
                                       t[1].cast<std::vector<std::string>>());
        });

// Expose functions. We use macros to avoid tedium
#define MAKE_POLYTABLE_MANIP_FUNCTION(NAME, FXN, TYPE)                        \
    m.def(NAME, &Sequence::FXN<TYPE>, py::arg("polytable"),                   \
          py::arg("skip_ancestral") = false, py::arg("ancestral") = 0,        \
          py::arg("gapchar") = "-");

    MAKE_POLYTABLE_MANIP_FUNCTION("removeMono", removeInvariantPos,
                                  Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeGaps", removeGaps, Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeAmbiguous", removeAmbiguous,
                                  Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMultiHits", removeMultiHits,
                                  Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMissing", removeMissing,
                                  Sequence::SimData);

    MAKE_POLYTABLE_MANIP_FUNCTION("removeMono", removeInvariantPos,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeGaps", removeGaps,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeAmbiguous", removeAmbiguous,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMultiHits", removeMultiHits,
                                  Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION("removeMissing", removeMissing,
                                  Sequence::PolySites);

    m.def("isValid", [](const Sequence::PolyTable& p) {
        return Sequence::polyTableValid(&p);
    });

    m.def("removeColumns",
          [](const Sequence::SimData& d,
             std::function<bool(const Sequence::stateCounter&)> f) {
              return Sequence::removeColumns(d, f);
          });

    m.def("removeColumns",
          [](const Sequence::PolySites& d,
             std::function<bool(const Sequence::stateCounter&)> f,
             const bool skip_ancestral, const unsigned anc,
             const char gapchar) { return Sequence::removeColumns(d, f); },
          py::arg("d"), py::arg("fxn"), py::arg("skip_ancestral") = false,
          py::arg("anc") = 0, py::arg("gapchar") = '-');
    return m.ptr();
}
