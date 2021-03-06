"""(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""
import re
import sys


def parse_clusters(clfile, num_repr=3, obligatory=()):
    with open(clfile) as input:
        clstrs = re.split(">Cluster", input.read())  # [1:num_repr+1]
    #        clstrs = map(lambda x:x.split("\n")[1:],re.split(">Cluster",input.read())[1:NUM_REPR+1])
    #        print clstrs[0][0]
    reprs = []
    for clstr in clstrs[1:num_repr + 1]:
        seqs = re.findall("\t([0-9]+)aa, >([A-Z0-9a-z|_/-]+)\.\.\.", clstr)
        in_clstr = filter(lambda x: not any(y in x[1] for y in obligatory),
                          sorted(seqs, reverse=True, key=lambda x: int(x[0])))
        reprs.append(in_clstr[0][1])
    if obligatory:
        for clstr in clstrs[1:]:
            seqs = re.findall("\t([0-9]+)aa, >([A-Z0-9a-z|_/-]+)\.\.\.", clstr)
            for s in filter(lambda x: any(y in x[1] for y in obligatory), seqs):
                reprs.append(s[1])
    return reprs


def get_seqs(fafile, reprs):
    all_seqs = []
    with open(fafile) as input:
        all_seqs = re.split("\n>", input.read()[1:])
    out = []
    for seq in all_seqs:
        if any(re.match("{}\n".format(x), seq) for x in reprs):
            out.append(seq)
    return out


def main(clfile, fafile, outfile="", num_repr=3, obligatory=()):
    obligatory = filter(lambda x: x, obligatory)
    out = get_seqs(fafile, parse_clusters(clfile, num_repr, obligatory))
    if outfile:
        with open(outfile, "w", 0) as output:
            output.write(">" + ("\n>".join(out)))
    else:
        print ">" + ("\n>".join(out))


if __name__ == "__main__":
    clfile = sys.argv[1]
    fafile = sys.argv[2]
    main(clfile, fafile)
