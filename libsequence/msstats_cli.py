# Python replacement for the C++ version of msstats.
from __future__ import print_function
import argparse
import libsequence
import libsequence.polytable as pt
import libsequence.summstats as sstats
import libsequence.citations as citations
import sqlite3
import collections
import pandas as pd


def make_parser():
    parser = argparse.ArgumentParser(description='Calculate summary statistics from data in "ms"-format',
                                     epilog=citations.LIBSEQUENCE)

    parser.add_argument(
        "-V", "--version", action='version',
        version='%(prog)s {}'.format(libsequence.__version__))
    parser.add_argument(
        "-v", "--verbose", action='store_true', help="Verbose output.  Extra info printed to standard error stream.")

    parser.add_argument("--garud", "-g", action='store_true',
                        help="Calculate H1, H12, etc.")
    parser.add_argument("--outfile", "-o", type=str,
                        help="sqlite3 output file name")
    return parser


ClassicStats = collections.namedtuple('ClassicStats', ['rep',
                                                       'thetapi', 'thetaw', 'thetah', 'tajd', 'S', 'singletons', 'dsingletons'])
GarudStatsT = collections.namedtuple(
    'GarudStatsT', ['rep', 'H1', 'H12', 'H2H1'])


def classic_stats(d,rep):
    ad = sstats.PolySIM(d)
    return ClassicStats(rep, ad.thetapi(),
                        ad.thetaw(),
                        ad.thetah(),
                        ad.tajimasd(),
                        ad.numpoly(),
                        ad.numsingletons(),
                        ad.numexternalmutations())


def msstats_main(arg_list=None):
    parser = make_parser()
    args = parser.parse_args(arg_list)

    d = pt.SimData()
    classic_stats_list = []
    garud_stats_list = []
    rep = 0
    while d.from_stdin() is True:
        classic_stats_list.append(classic_stats(d,rep))
        if args.garud is True:
            gstats = sstats.garudStats(d)
            garud_stats_list.append(GarudStatsT(
                rep, gstats['H1'], gstats['H12'], gstats['H2H1']))
        rep += 1
    df = pd.DataFrame(classic_stats_list).set_index('rep')
    if len(garud_stats_list) > 0:
        gdf = pd.DataFrame(garud_stats_list).set_index('rep')
        df = df.join(gdf)
    con = sqlite3.connect(args.outfile)
    df.to_sql('stats', con, index=False)
    con.close()
