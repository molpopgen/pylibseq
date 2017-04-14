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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <Sequence/PolyTable.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/polySiteVector.hpp>

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

    return m.ptr();
}
