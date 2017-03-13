#Python replacement for the C++ version of msstats.
from __future__ import print_function
import argparse
import libsequence
import libsequence.polytable as pt
import libsequence.summstats as sstats
import libsequence.console.citations as citations
import sqlite3
import pandas as pd

def make_parser():
    parser = argparse.ArgumentParser(description='Calculate summary statistics from data in "ms"-format',
            epilog=citations.LIBSEQUENCE)

    parser.add_argument(
    "-V", "--version", action='version',
    version='%(prog)s {}'.format(libsequence.__version__))
    parser.add_argument(
    "-v", "--verbose", action='store_true',help="Verbose output.  Extra info printed to standard error stream.")

    parser.add_argument("--garud","-g",action='store_true',help="Calculate H1, H12, etc.")
    parser.add_argument("--outfile","-o",type=str,help="sqlite3 output file name")
    return parser

def classic_stats(d):
    ad = sstats.PolySIM(d)
    return {'thetapi':ad.thetapi(),
            'thetaw':ad.thetaw(),
            'thetah':ad.thetah(),
            'tajd':ad.tajimasd(),
            'S':ad.numpoly(),
            'singletons':ad.numsingletons(),
            'dsingletons':ad.numexternalmutations()}

def msstats_main(arg_list=None):
    parser=make_parser()
    args=parser.parse_args(arg_list)

    d=pt.SimData()
    data=[]
    rep=0
    while pt.readSimData(d) is True:
        classic=classic_stats(d)
        if args.garud is True:
            gstats = sstats.garudStats(d)
            classic.update(gstats)
        classic['rep']=rep
        rep+=1
        data.append(classic)
    df=pd.DataFrame(data)
    con=sqlite3.connect(args.outfile)
    df.to_sql('stats',con,index=False)
    con.close()
