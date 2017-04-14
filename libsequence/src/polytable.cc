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

PYBIND11_PLUGIN(polytable)
{
    py::module m("polytable", "Access to libsequence's polymorphism table "
                              "classes and related functions");

    py::class_<Sequence::PolyTable>(m, "PolyTable",
                                    "Base class for polymorphism tables")
        .def("GetData", &Sequence::PolyTable::GetData)
        .def("GetPositions", &Sequence::PolyTable::GetPositions)
        .def("assign",
             [](Sequence::PolyTable& p, const Sequence::polySiteVector& v) {
                 p.assign(v.cbegin(), v.cend());
             })
        .def("assign",
             [](Sequence::PolyTable& p, std::vector<double> pos,
                std::vector<std::string> d) {
                 p.assign(std::move(pos), std::move(d));
             })
        .def("empty", &Sequence::PolyTable::empty)
        .def("numsites", &Sequence::PolyTable::numsites)
        .def("size", &Sequence::PolyTable::size)
        .def("__len__", [](const Sequence::PolyTable& p) { return p.size(); });

    py::class_<Sequence::PolySites, Sequence::PolyTable>(
        m, "PolySites", "A polymorphism table for Sequence data.  "
                        "A/G/C/T/N/0/1 is the alphabet for analysis.")
        .def("__init__",
             [](Sequence::PolySites& p, const Sequence::polySiteVector& v) {
                 new (&p) Sequence::PolySites(v.cbegin(), v.cend());
             });

    py::class_<Sequence::SimData, Sequence::PolyTable>(
        m, "SimData",
        "A polymorphism table for binary data.  0/1 = ancestral/derived.")
        .def(py::init<std::vector<double>, std::vector<std::string>>())
        .def("__init__",
             [](Sequence::SimData& d, const Sequence::polySiteVector& p) {
                 new (&d) Sequence::SimData(p.cbegin(), p.cend());
             });

    py::class_<Sequence::stateCounter>(m, "StateCounter")
        .def(py::init<char>(), py::arg("gapchar") = '-')
        .def_readonly("zero", &Sequence::stateCounter::zero)
        .def_readonly("one", &Sequence::stateCounter::one)
        .def_readonly("a", &Sequence::stateCounter::a)
        .def_readonly("g", &Sequence::stateCounter::g)
        .def_readonly("c", &Sequence::stateCounter::c)
        .def_readonly("t", &Sequence::stateCounter::t)
        .def_readonly("ndna", &Sequence::stateCounter::ndna)
        .def("nstates", &Sequence::stateCounter::nStates)
        .def("__call__",
             [](Sequence::stateCounter& c, const char ch) { c(ch); });

// Expose functions. We use macros to avoid tedium
#define MAKE_POLYTABLE_MANIP_FUNCTION(FXN, TYPE)                              \
    m.def(" FXN ", &Sequence::FXN<TYPE>);

    MAKE_POLYTABLE_MANIP_FUNCTION(removeInvariantPos, Sequence::SimData);
    MAKE_POLYTABLE_MANIP_FUNCTION(removeInvariantPos, Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION(removeGaps, Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION(removeAmbiguous, Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION(removeMultiHits, Sequence::PolySites);
    MAKE_POLYTABLE_MANIP_FUNCTION(removeMissing, Sequence::PolySites);

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
