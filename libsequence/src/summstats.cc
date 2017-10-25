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
#include <Sequence/stateCounter.hpp>

namespace py = pybind11;

PYBIND11_MODULE(summstats, m)
{
    m.doc() = "Summary statistics";

    py::object polytable
        = (py::object)py::module::import("libsequence.polytable")
              .attr("PolyTable");

    py::class_<Sequence::PolySNP>(m, "PolySNP",
                                  "Class to calculate summary statistics.")
        .def("__init__",
             [](Sequence::PolySNP& newobj, const Sequence::PolyTable& pt,
                const bool haveOutgroup, const unsigned outgroup,
                const bool totMuts) {
                 new (&newobj)
                     Sequence::PolySNP(&pt, haveOutgroup, outgroup, totMuts);
             },
             py::arg("pt"), py::arg("haveOutgroup") = false,
             py::arg("outgroup") = 0, py::arg("totMuts") = true)
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
             py::arg("likeThorntonAndolfatto") = false,
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

			.. note::
				Public API shared with :class:`libsequence.summstats.PolySNP`.
            )delim")
        .def("__init__",
             [](Sequence::PolySIM& newobj, const Sequence::SimData& d) {
                 new (&newobj) Sequence::PolySIM(&d);
             });

    py::class_<Sequence::stateCounter>(
        m, "StateCounter", "Tally up the number of occurrences of value "
                           "polymorphism characters at a site.")
        .def(py::init<char>(), py::arg("gapchar") = '-')
        .def_readonly("zero", &Sequence::stateCounter::zero,
                      "Number of times '0' was seen at a site.")
        .def_readonly("one", &Sequence::stateCounter::one,
                      "Number of times '1' was seen at a site.")
        .def_readonly("a", &Sequence::stateCounter::a,
                      "Number of times 'A' or 'a' was seen at a site.")
        .def_readonly("g", &Sequence::stateCounter::g,
                      "Number of times 'G' or 'g' was seen at a site.")
        .def_readonly("c", &Sequence::stateCounter::c,
                      "Number of times 'C' or 'c' was seen at a site.")
        .def_readonly("t", &Sequence::stateCounter::t,
                      "Number of times 'T' or 't' was seen at a site.")
        .def_readonly(
            "ndna", &Sequence::stateCounter::ndna,
            "Number of times a non-DNA character was seen at a site.")
        .def("nstates", &Sequence::stateCounter::nStates,
             "The total number of character states observed at a site.")
        .def("__call__", [](Sequence::stateCounter& c, const std::string& s) {
            std::for_each(s.begin(), s.end(), [&c](const char ch) { c(ch); });
        });

    m.def(
        "nSLiHS",
        [](const Sequence::SimData& d, py::object core_snps, py::object gmap) {
            if (gmap.is_none())
                {
                    if (core_snps.is_none())
                        {
                            return Sequence::nSL_t(d);
                        }
                    else
                        {
                            return Sequence::nSL_t(
                                d, core_snps.cast<std::vector<std::size_t>>());
                        }
                }
            if (core_snps.is_none())
                {
                    return Sequence::nSL_t(
                        d, gmap.cast<std::unordered_map<double, double>>());
                }
            return Sequence::nSL_t(
                d, core_snps.cast<std::vector<std::size_t>>(),
                gmap.cast<std::unordered_map<double, double>>());
        },
        R"delim(
		"Raw"/unstandardized :math:`nS_L` and iHS from Ferrer-Admetlla et al. doi:10.1093/molbev/msu077.

		:param pt: A :class:`libsequence.polytable.PolyTable`
        :param core_snps: (None) Indexes of SNPs to analyze as "core" SNPs.
		:param gmap: (None) A dictionary relating each position in pt to its location on a genetic map.
		:return: A list of (nSL,iHS) tuples
		:rtype: list
    
		.. note:: Only :class:`libsequence.polytable.SimData` types currently supported

		)delim",
        py::arg("d"), py::arg("core_snps") = nullptr,
        py::arg("gmap") = nullptr);

    m.def("lhaf", &Sequence::lHaf,
          R"delim(
		:math:`l-HAF` from Ronen et al. DOI:10.1371/journal.pgen.1005527
    
		:param pt: A :class:`libsequence.polytable.PolyTable`
		:param l: The scaling factor for the statistic. See paper for details.
		:return: The :math:`l-HAF` statistic for each haplotype in pt
    
		:rtype: list

		.. note:: Only :class:`libsequence.polytable.SimData` types currently supported.
		)delim");

    m.def("ld",
          [](const Sequence::PolyTable& p, const bool have_outgroup,
             const unsigned outgroup, const unsigned mincount,
             const double maxd) {
              auto temp = Sequence::Recombination::Disequilibrium(
                  &p, have_outgroup, outgroup, mincount, maxd);
              // Before filling a py::list, let's get rid of skipped objects
              temp.erase(
                  std::remove_if(temp.begin(), temp.end(),
                                 [](const Sequence::PairwiseLDstats& s) {
                                     return s.skipped;
                                 }),
                  temp.end());
              py::list rv;
              for (auto&& ld : temp)
                  {
                      py::dict temp;
                      temp[py::str("i")] = py::float_(ld.i);
                      temp[py::str("j")] = py::float_(ld.j);
                      temp[py::str("rsq")] = py::float_(ld.rsq);
                      temp[py::str("D")] = py::float_(ld.D);
                      temp[py::str("Dprime")] = py::float_(ld.Dprime);
                      rv.append(std::move(temp));
                  }
              return rv;
          },
          R"delim(
		Return pairwise LD statistics.
		
		:param p: A :class:`libsequence.polytable.PolySites` or :class:`libsequence.polytable.SimData`.
		:param have_outgroup: (False) Boolean--is outgroup sequence present in p?
		:param outgroup: (0) The index of the outgroup sequence in p.
		:param mincount: Do not include site pairs with MAF < mincount.
		:param maxd: Do not include site pairs separated by >= maxd.

		:rtype: list of dict
		)delim",
          py::arg("p"), py::arg("have_outgroup") = false,
          py::arg("outgroup") = 0, py::arg("mincount") = 1,
          py::arg("maxd") = std::numeric_limits<double>::max());

    m.def("garudStats",
          [](const Sequence::SimData& d) {
              auto g = Sequence::H1H12(d);
              py::dict rv;
              rv[py::str("H1")] = py::float_(g.H1);
              rv[py::str("H12")] = py::float_(g.H12);
              rv[py::str("H2H1")] = py::float_(g.H2H1);
              return rv;
          },
          R"delim(
		Returns the H1, H12, and H2/H1 statistics from PMC4338236 as a dict.

		:param d: A :class:`libsequence.polytable.SimData`.
		)delim",
          py::arg("d"));
}
