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

    py::class_<Sequence::PolySNP>(m, "PolySNP",
                                  "Class to calculate summary statistics.")
        .def("__init__",
             [](Sequence::PolySNP& newobj, const Sequence::PolyTable& pt,
                const bool haveOutgroup, const unsigned outgroup,
                const bool totMuts) {
                 new (&newobj)
                     Sequence::PolySNP(&pt, haveOutgroup, outgroup, totMuts);
             })
        .def("thetapi", &Sequence::PolySNP::ThetaPi,
             "Sum of site heterozygosity/mean pairwise differences.")
        .def("thetaw", &Sequence::PolySNP::ThetaW,
             R"d(Watterson's :math:`\hat\theta`)d")
        .def("thetah", &Sequence::PolySNP::ThetaH,
             R"d(Fay and Wu's :math:`\hat\theta_H`)d")
        .def(
            "thetal", &Sequence::PolySNP::ThetaL,
            R"d(:math:`\hat\theta` from site homozygosity. See eqn. 1 of PMC1456302.)d")
        .def("numpoly", &Sequence::PolySNP::NumPoly,
             "Total number of variable positions.")
        .def("numsingletons", &Sequence::PolySNP::NumExternalMutations,
             "Number of singletons.")
        .def("numexternalmutations", &Sequence::PolySNP::NumExternalMutations,
             "Number of derived singletons.")
        .def("tajimasd", &Sequence::PolySNP::TajimasD, "Tajima's D.")
        .def("hprime", &Sequence::PolySNP::Hprime,
             R"d(Normalized version of Fay and Wu's :math:`\hat\theta_H`)d")
        .def("fulid", &Sequence::PolySNP::FuLiD, "Fu and Li's D")
        .def("fulif", &Sequence::PolySNP::FuLiF, "Fi and Li's F")
        .def("fulidstar", &Sequence::PolySNP::FuLiDStar, "Fu and Li's D*")
        .def("fulifstar", &Sequence::PolySNP::FuLiFStar, "Fu and Li's F*")
        .def("hapdiv", &Sequence::PolySNP::DandVH, "Haplotype diversity")
        .def("nhaps", &Sequence::PolySNP::DandVK, "Number of haplotypes")
        .def("wallsb", &Sequence::PolySNP::WallsB, "Wall's B")
        .def("wallsq", &Sequence::PolySNP::WallsQ, "Wall's Q")
        .def("wallsbprime", &Sequence::PolySNP::WallsBprime, "Wall's B'")
        .def("rm", &Sequence::PolySNP::Minrec, "Hudson and Kaplan's lower "
                                               "bound on the number of "
                                               "recombination events.");

    py::class_<Sequence::PolySIM, Sequence::PolySNP>(m, "PolySIM",
                                                     R"delim(
            Class to calculate summary statistics for 0/1-encoded data
            in a :class:`libsequence.polytable.SimData` object.
            )delim")
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

    py::class_<Sequence::GarudStats>(m, "GarudStats")
        .def_readonly("H1", &Sequence::GarudStats::H1)
        .def_readonly("H12", &Sequence::GarudStats::H12)
        .def_readonly("H2H1", &Sequence::GarudStats::H2H1);

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
    m.def("garudStats", &Sequence::H1H12);

    return m.ptr();
}
