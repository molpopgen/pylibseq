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

namespace py = pybind11;

PYBIND11_PLUGIN(polytable)
{
    py::module m("polytable","Access to libsequence's polymorphism table classes and related functions");

    return m.ptr();
}
