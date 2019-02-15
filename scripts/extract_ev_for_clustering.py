"""(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""
import argparse
import glob

from my_little_hhpred_reader import HHpredOutput


###USAGE
## takes hhsearch all vs all results directory
## creates file "all_vs_all.out"

def parse_results(HHRS, outname="all_vs_all.out", ev_cutoff=.001):
    RESULTS = []
    # parser = HHOutputParser(False)
    for hhr in HHRS:
        results = HHpredOutput(hhr)  # parser.parse_file(hhr)
        this = results.query  # _query_name
        for hit in results.hits:
            id = hit.target  # _id
            if id == this:
                continue
            # eva =  hit._evalue
            if hit.eval < ev_cutoff:
                RESULTS.append((this, id, hit.eval))
    print RESULTS
    with open(outname, "w") as out:
        out.write("\n".join(map(lambda x: "{}\t{}\t{}".format(*x), RESULTS)))

    return RESULTS


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Parser for HHsearch results, returns a tsv file with extracted e-values')
    parser.add_argument('-e', '--evalue', action="store", dest="EV_CUTOFF", type=float, default=0.0001,
                        help='Set the e-value cutoff, defaults to 1e-4')
    parser.add_argument('-o', '--output', action="store", dest="output", default='all_vs_all.out',
                        help='Output filename')
    parser.add_argument('working_directory', action="store",
                        help='Directory in which we should look for .hhr files to parse')
    args = parser.parse_args()

    HHRS = glob.glob(args.working_directory + "/*.hhr")

    parse_results(HHRS, args.output, args.EV_CUTOFF)

#    print HHRS
