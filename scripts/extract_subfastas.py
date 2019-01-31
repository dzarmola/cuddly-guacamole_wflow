#!/usr/bin/python
"""(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl"""

from my_little_hhpred_reader import *
import json
import subprocess
import os,re
cwd = os.getcwd()

def coreToSeq(seq,seq_start,core_start):
    """ Given seq and index of the core in this sequence,
        return alignment index of the core"""
    core_start-=seq_start
    if core_start < 0:
        return (0,0)
    if core_start == 0:
        return (0,1)
    for i,s in enumerate(seq):
        if not core_start:
            return (i,1)
        if s not in ["-","."] and s==s.upper():
            core_start -= 1
    else:
        return (len(seq)-1,0)

def seqToCore(seq,seq_index,seq_start):
    """ Given seq_index of the core limit in the matching sequence,
    returns the index of the expected core in the given seq"""

    prefix = seq[:seq_index].replace("-","").replace(".","")
    seq_start += len(prefix)
    return seq_start

def overlap(p1,p2,s,e):
    i = set(range(p1,p2+1))
    j = set(range(s,e))
    return len(i.intersection(j))

def uuu(p, s,e):
    return p[overlap(p[0][0],p[0][1],s,e)<overlap(p[1][0],p[1][1],s,e)]

def parse_results(HHRS,cores):
    RESULTS = {'length':0,'hits':[],'name':""}
    for hhr in HHRS:
        results = HHpredOutput(hhr)
        this=results.query
        RESULTS['name'] = this
        for hit in results.hits:
            id = hit.target
            RESULTS['length'] = len(hit.match)
            if id==this:
                continue
            pfam_s = hit.q_cons_s
            pfam_e = hit.q_cons_e
            pdb_s = hit.t_cons_s
            pdb_e = hit.t_cons_e
            pfam_seq = hit.q_cons #hit.q_better_ref #hit._alignment._msa['query'].sequence
            pdb_seq = hit.t_cons #hit.t_better_ref #hit._alignment._msa[hit._id].sequence

            cs = cores[id.split("_")[0].lower()]
            cs = uuu(cs,pdb_s,pdb_e)
    
            i,er0 = coreToSeq(pdb_seq,pdb_s,cs[0])
            j,er1 = coreToSeq(pdb_seq,pdb_s,cs[1])

            s = seqToCore(pfam_seq,i,pfam_s)
            e = seqToCore(pfam_seq,j,pfam_s)
            match = hit.q_better_ref[i:j+1]
            pdb_core = hit.t_better_ref[i:j+1]
    
            def proper_length(seq):
                 return len(seq.replace("-","").replace(".",""))
            RESULTS['hits'].append((id,(s,e),(er0,er1),hit.eval,(s,match,e,pdb_s+i,pdb_core,pdb_s+proper_length(pdb_core)) ))

    return RESULTS

def overlap(s1,e1,s2,e2,cutoff=0.5):
    l1 = e1-s1
    l2 = e2-s2
    return  len( set(range(s1,e1+1)) & set(range(s2,e2+1))) > min(l1,l2)*cutoff

def combine(s1,e1,s2,e2):
    return [min(s1,e1,s2,e2), max(s1,e1,s2,e2)]

def fix_results(res):
    EXITABLE = 0
    seqs = []
    for k in res['hits']:
        if k[3] <= 0.001 and k[1][1]-k[1][0]>40:
            seqs.append([list(k[2]),list(k[4][:3]),k[3]])
    if not seqs:
        EXITABLE = 1
        for k in res['hits']:
            if k[3] < 0.01:
                seqs.append([list(k[2]),list(k[4][:3]),k[3]])
        if not seqs:
            exit("{} No cores <.001 found!".format(sys.argv[1]))


    seqs = [_[1] for _ in sorted(seqs, key=lambda x:x[2]) ] 
    s = seqs[0]
    kombajn = [s[:2]]
    r = [s[0],s[2]]
    for _ in seqs[1:]:
        _r = [_[0],_[2]]
        if overlap(*(r+_r)):
            kombajn.append(_[:2])
            r = combine(*(r+_r))

    return (r[0],r[1],EXITABLE)

def extract_seq_from_hhm(hhm,s,e):
    seq = ""
    match_state = re.compile("^([A-Z-])\s([0-9-]+)\s")
    with open(hhm) as input:
        for line in input:
            if match_state.match(line):
                aa,num = match_state.findall(line)[0]
                num = int(num)
                if s<=num and num <=e:
                    seq+=aa
    return seq


def main(hhr,json_file):

    with open(json_file) as infile:
        cores = json.load(infile)#{ "4n7w": [(10,87),(163,242)] }

    res = parse_results([hhr],cores)
    s,e,ok= fix_results(res)

    return s,e

if __name__ == "__main__":
    import sys
    hhr = sys.argv[1]
    with open(sys.argv[2]) as infile:
        positions = json.load(infile)

    fam = hhr.split("/")[1]

    for part in positions.keys(): #["core"]:#"N C core".split():
        cores = positions[part] # { "4n7w": [(10,87),(163,242)] }
        res = parse_results([hhr],cores)
        s,e,ok= fix_results(res)
        subprocess.call("{}/scripts/extract_fasta.sh {} {} {} {}".format(cwd,fam,s,e,part),shell=True)
