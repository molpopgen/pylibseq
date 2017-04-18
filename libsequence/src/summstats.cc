#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/SummStats/nSL.hpp>
#include <Sequence/SummStats/Garud.hpp>
#include <Sequence/SummStats/lHaf.hpp>
#include <Sequence/Recombination.hpp>

namespace py = pybind11;

PYBIND11_PLUGIN(summstats)
{
    py::module m("summstats", "Summary statistics");

    py::class_<Sequence::PolySNP>(m, "PolySNP")
        .def("__init__",
             [](Sequence::PolySNP& newobj, const Sequence::PolyTable& pt,
                const bool haveOutgroup, const unsigned outgroup,
                const bool totMuts) {
                 new (&newobj)
                     Sequence::PolySNP(&pt, haveOutgroup, outgroup, totMuts);
             })
        .def("thetapi", &Sequence::PolySNP::ThetaPi)
        .def("thetaw", &Sequence::PolySNP::ThetaW)
        .def("thetah", &Sequence::PolySNP::ThetaH)
        .def("thetapi", &Sequence::PolySNP::ThetaPi)
        .def("numpoly", &Sequence::PolySNP::NumPoly)
        .def("numexternalmutations", &Sequence::PolySNP::NumExternalMutations)
        .def("tajimasd", &Sequence::PolySNP::TajimasD)
        .def("hprime", &Sequence::PolySNP::Hprime)
        .def("fulid", &Sequence::PolySNP::FuLiD)
        .def("fulif", &Sequence::PolySNP::FuLiF)
        .def("fulidstar", &Sequence::PolySNP::FuLiDStar)
        .def("fulifstar", &Sequence::PolySNP::FuLiFStar)
        .def("hapdiv", &Sequence::PolySNP::DandVH)
        .def("nhaps", &Sequence::PolySNP::DandVK)
        .def("wallsb", &Sequence::PolySNP::WallsB)
        .def("wallsq", &Sequence::PolySNP::WallsQ)
        .def("wallsbprime", &Sequence::PolySNP::WallsBprime)
        .def("rm", &Sequence::PolySNP::Minrec);

    py::class_<Sequence::PolySIM, Sequence::PolySNP>(m, "PolySIM")
        .def("__init__",
             [](Sequence::PolySIM& newobj, const Sequence::SimData& d) {
                 new (&newobj) Sequence::PolySIM(&d);
             });

    py::class_<Sequence::PairwiseLDstats>(m, "PairwiseLDstats")
        .def_readonly("i", &Sequence::PairwiseLDstats::i)
        .def_readonly("j", &Sequence::PairwiseLDstats::j)
        .def_readonly("rsq", &Sequence::PairwiseLDstats::rsq)
        .def_readonly("D", &Sequence::PairwiseLDstats::D)
        .def_readonly("Dprime", &Sequence::PairwiseLDstats::Dprime)
        .def_readonly("skipped", &Sequence::PairwiseLDstats::Dprime);

    py::class_<Sequence::GarudStats>(m,"GarudStats")
        .def_readonly("H1",&Sequence::GarudStats::H1)
        .def_readonly("H12",&Sequence::GarudStats::H12)
        .def_readonly("H2H1",&Sequence::GarudStats::H2H1);

    m.def("nSLiHS",
          [](const Sequence::SimData& d) { return Sequence::nSL_t(d); });
    m.def("lhaf", &Sequence::lHaf);
    m.def("std_nSLiHS",
          [](const Sequence::SimData& d, const double minfreq,
             const double binsize) {
              return Sequence::snSL(d, minfreq, binsize);
          },
          py::arg("d"), py::arg("minfreq") = 0.0, py::arg("binsize") = 0.05);
    m.def("ld", &Sequence::Recombination::Disequilibrium);
    m.def("garudStats",&Sequence::H1H12);

    return m.ptr();
}
