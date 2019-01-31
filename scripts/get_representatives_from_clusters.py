"""(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""
import sys
import re

NUM_REPR = 3

def parse_clusters(clfile):
    with open(clfile) as input:
        clstrs = re.split(">Cluster",input.read())[1:NUM_REPR+1]
#        clstrs = map(lambda x:x.split("\n")[1:],re.split(">Cluster",input.read())[1:NUM_REPR+1])
#        print clstrs[0][0]
    reprs = []
    for clstr in clstrs:
        seqs = re.findall("\t([0-9]+)aa, >([A-Z0-9a-z_]+)\.\.\.",clstr)
        reprs.append(sorted(seqs, reverse =True, key=lambda x: int(x[0]))[0][1])
    return reprs

def get_seqs(fafile, reprs):
    all_seqs = []
    with open(fafile) as input:
        all_seqs = re.split("\n>",input.read()[1:])
    out = []
    for seq in all_seqs:
        if any(re.match("{}\n".format(x),seq) for x in reprs):
            out.append(seq)
    return out

def main(clfile,fafile,outfile=""):
    out = get_seqs(fafile, parse_clusters(clfile))
    if outfile:
        with open(outfile,"w",0) as output:
            output.write(">"+("\n>".join(out)))
    else:
        print ">" + ("\n>".join(out))


if __name__ == "__main__":
    clfile = sys.argv[1]
    fafile = sys.argv[2]
    main(clfile, fafile)
