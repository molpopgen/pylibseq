#include <Sequence/SimData.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/SimParams.hpp>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;
using namespace Sequence;

std::pair<double, double>
omega_max(const SimData& data)
{
    if (data.empty())
        {
            return make_pair(numeric_limits<double>::quiet_NaN(),
                             numeric_limits<double>::quiet_NaN());
        }
    vector<unsigned> mafs;
    for (auto i = data.sbegin(); i < data.send(); ++i)
        {
            unsigned c = count(i->second.begin(), i->second.end(), '1');
            unsigned maf = min(c, unsigned(data.size()) - c);
            mafs.push_back(maf);
        }
    PolySNP adata(&data);
    vector<PairwiseLDstats> ld = adata.Disequilibrium(2);

    if (ld.empty())
        {
            return make_pair(numeric_limits<double>::quiet_NaN(),
                             numeric_limits<double>::quiet_NaN());
        }
    double omega_max = numeric_limits<double>::min();
    double snp = numeric_limits<double>::quiet_NaN();
    unsigned S = data.numsites();
    for (unsigned pos = 1; pos < S - 1; ++pos)
        {
            if (mafs[pos] > 1)
                {
                    unsigned l = pos + 1;
                    double position = data.position(pos);
                    double srsqL = 0., srsqR = 0., sumrsqLR = 0.;
                    for (unsigned i = 0; i < ld.size(); ++i)
                        {
                            if (ld[i].i <= position)
                                {
                                    srsqL += ld[i].rsq;
                                }
                            else if (ld[i].j > position)
                                {
                                    srsqR += ld[i].rsq;
                                }
                            if (ld[i].i <= position && ld[i].j > position)
                                {
                                    sumrsqLR += ld[i].rsq;
                                }
                        }
                    double numerator
                        = (1. / (double(l * (l - 1)) / 2.
                                 + double((S - l) * (S - l - 1)) / 2.))
                          * (srsqL + srsqR);
                    double denominator
                        = (1. / (double(l * (S - l)))) * sumrsqLR;
                    double omega = numerator / denominator;
                    if (finite(omega) && omega > omega_max)
                        {
                            omega_max = omega;
                            snp = position;
                        }
                }
        }
    return make_pair(omega_max, snp);
}
