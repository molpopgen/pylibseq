#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Sequence/FST.hpp>
#include <Sequence/PolyTable.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(fst)
{
    py::module m("fst", "Fst");

    py::object polytable
        = (py::object)py::module::import("libsequence.polytable")
              .attr("PolyTable");

    py::class_<Sequence::FST>(m, "Fst", "Fst")
        .def("__init__",
             [](Sequence::FST& fst, const Sequence::PolyTable& pt,
                const std::vector<unsigned>& config, const bool haveOutgroup,
				const unsigned outgroup)
			 {
				 new (&fst) Sequence::FST(&pt, config.size(), &config[0],nullptr,haveOutgroup,outgroup);
             },
             py::arg("polytable"), py::arg("config"),py::arg("haveOutgroup")=false,py::arg("outgroup")=0)
        .def("hsm", &Sequence::FST::HSM)
        .def("slatkin", &Sequence::FST::Slatkin)
        .def("hbk", &Sequence::FST::HBK)
        .def("piB", &Sequence::FST::piB)
        .def("piS", &Sequence::FST::piS)
        .def("piT", &Sequence::FST::piT)
        .def("piD", &Sequence::FST::piD)
        .def("shared", &Sequence::FST::shared)
        .def("fixed", &Sequence::FST::fixed)
        .def("priv", &Sequence::FST::Private);

    return m.ptr();
}
