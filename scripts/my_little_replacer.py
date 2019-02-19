"""(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""

import glob
import re
import sys


def from_save(filename):
    with open(filename) as save:
        clusters = eval(save.readline())
        labels = eval(save.readline())
        data = eval(save.readline())
    #        clusters,labels,data = map(eval,save.readlines())
    return data, labels, clusters


def read_fasta(fname):
    seqs = []
    with open(fname) as input:
        seqs = re.findall(">([A-Za-z0-9/|_-]+)\n([\.a-zA-Z\n-]+)", input.read())
    return fname.split("/")[-1].replace("_representatives.fa", ""), [(x[0], x[1].replace("\n", "")) for x in seqs]


def replace(seq, aligned):
    return "".join(["-" if s is None else seq[s - 1] for s in aligned])


def get(dict, key):
    _tmp = None
    for k, v in dict.items():
        if k == key:
            return v
        if key in k or k in key:
            _tmp = v
    if _tmp is not None: return _tmp
    raise IndexError("No such key:{}".format(key))


def main(sfile, rpdir, outfile=''):
    data, labels, clusters = from_save(sfile)
    representatives = {}
    for file in glob.glob(rpdir + "/*representatives.fa"):
        name, seqs = read_fasta(file)
        representatives[name] = seqs

    #print len(representatives.keys()),representatives.keys()
    #print representatives["PF03788"]

    if outfile:
        with open(outfile, "w", 0) as output:
            for l, d in zip(labels, data):
                for n, s in get(representatives, l):
                    output.write(">{}\n{}\n".format(n, replace(s, d)))
    else:
        for l, d in zip(labels, data):
            #            print l , get(representatives, l)
            #print l,len(get(representatives, l))
            for n, s in get(representatives, l):
                print ">{}\n{}".format(n, replace(s, d))


if __name__ == "__main__":
    if len(sys.argv) > 3:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[1], sys.argv[2])
