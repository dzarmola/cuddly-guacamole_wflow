"""(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""
import sys


def read_fasta(fname):
    cols = []
    names = []
    cnt = 1
    with open(fname) as input:
        id = ''
        seq = ''
        for line in input:
            if line[0] == ">":
                if id:
                    names.append(id)
                    for i, c in enumerate(seq):
                        if i < len(cols):
                            cols[i] += c
                        else:
                            cols.append(c)
                    seq = ''
                id = "{}_{}".format(line.strip(), cnt)
                cnt += 1
            else:
                seq += line.strip()
        if id:
            names.append(id)
            for i, c in enumerate(seq):
                if i < len(cols):
                    cols[i] += c
                else:
                    cols.append(c)
    return names, cols


def length(x):
    cnt = 0
    for i in x:
        cnt += 1
    return cnt


def remove_by_cutoff(cols, cutoff):
    lc = len(cols[0]) * 1.  # /100.
    new_cols = filter(lambda x: length(_ for _ in x if _ != "-") / lc >= cutoff, cols)
    # print sorted(map(lambda x: round(x,2),[length(_ for _ in x if _ != "-")/lc for x in cols]))
    #    print cutoff,lc
    #    print new_cols
    return new_cols


def remake_seqs(ncols):
    outseqs = ["" for x in ncols[0]]
    for col in ncols:
        for i, c in enumerate(col):
            outseqs[i] += c
    return outseqs


def write_out(names, strings, outname):
    with open(outname, "w", 0) as out:
        for n, s in zip(names, strings):
            # if n[-2]=="_": n = n[:-2]
            out.write("{}\n{}\n".format(n, s))
    coutname = outname.replace(".fa", "_cdhit.fa")
    with open(coutname, "w", 0) as out:
        for n, s in zip(names, strings):
            # if n[-2]=="_1": n = n[:-2]
            out.write("{}\n{}\n".format(n, s.replace("-", "").replace(".","")))


def main(fname, cutoff, family=''):
    outname = fname.replace(".fa", "_{}.fa".format(cutoff))
    if cutoff > 1:
        cutoff /= 100.
    names, cols = read_fasta(fname)
    if family:
        for i, n in enumerate(names):
            if family in n:
                break
            else:
                names[i] = ">{}_{}".format(family, n.strip(">"))
    #    print cols
    new_cols = remove_by_cutoff(cols, cutoff)
    #    outname = fname.replace(".fa","_{}.fa".format(cutoff))
    write_out(names, remake_seqs(new_cols), outname)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
