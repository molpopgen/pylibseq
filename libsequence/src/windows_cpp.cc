#include <pybind11/pybind11.h>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolyTableSlice.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(windows_cpp)
{
    py::module m("windows_cpp", "Details of sliding window containers.");

    using SimDataWindows = Sequence::PolyTableSlice<Sequence::SimData>;
    using PolySitesWindows = Sequence::PolyTableSlice<Sequence::PolySites>;

    py::class_<SimDataWindows>(m, "SimDataWindows")
        .def("__init__",
             [](SimDataWindows& sdw, const Sequence::polySiteVector& p,
                const double window_size, const double step_len,
                const double starting_pos, const double ending_pos) {
                 new (&sdw) SimDataWindows(p.cbegin(), p.cend(), window_size,
                                           step_len, starting_pos, ending_pos);
             })
        .def("__getitem__", [](const SimDataWindows& sdw,
                               const unsigned i) { return sdw[i]; })
        .def("__len__", [](const SimDataWindows& swd) { return swd.size(); });

    py::class_<PolySitesWindows>(m, "PolySitesWindows")
        .def("__init__",
             [](PolySitesWindows& sdw, const Sequence::polySiteVector& p,
                const double window_size, const double step_len,
                const double starting_pos, const double ending_pos) {
                 new (&sdw)
                     PolySitesWindows(p.cbegin(), p.cend(), window_size,
                                      step_len, starting_pos, ending_pos);
             })
        .def("__getitem__", [](const PolySitesWindows& sdw,
                               const unsigned i) { return sdw[i]; })
        .def("__len__",
             [](const PolySitesWindows& swd) { return swd.size(); });

    return m.ptr();
}
