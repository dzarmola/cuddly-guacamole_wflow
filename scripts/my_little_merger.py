import time
import sys,glob
from my_little_hhpred_reader import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.patches import Rectangle as rect
from matplotlib.lines import Line2D as line
import argparse
EVAL_CUTOFF = 0.001

def load_data(files):

    matrix = {}

    hhrs = []
    lens = {}
    for file in files:
        h = HHpredOutput(file)
        hhrs.append(h)
        lens[h.query] = h.len
        for _h in h.hits:
            if _h.eval<EVAL_CUTOFF:
                matrix[h.query] = matrix.get(h.query,{}) 
                matrix[_h.target] = matrix.get(_h.target,{}) 
                if _h.eval < matrix[_h.target].get(h.query,[1000,0,0])[0]:
                    matrix[_h.target][h.query] = (_h.eval,_h.t_cons,_h.q_cons)
                    matrix[h.query][_h.target] = (_h.eval,_h.q_cons,_h.t_cons)

    return matrix,hhrs,lens

def sort_fams_by_dist(mat):
    fams = sorted(mat.keys(), key=lambda x: sum(_[0] for _ in mat[x].values())/len(mat[x].values()) if mat[x] else 1000 )
    for key in mat.values():
        for k in sorted(key,key=lambda x:x[0]):
            if k not in fams:
                fams.append(k)
    return fams

def read_in_clusters(file):
    cluster = {}
    with open(file) as input:
        for i,line in enumerate(input):
            for fam in line.split():
                cluster[fam] = i
    return cluster

def save_point(data,labels,clusters,lens):
#    timestr = time.strftime("%Y%m%d-%H%M%S")
#    with open("savepoint_{}.txt".format(timestr),"w",0) as output:
    global savepoint_name
    with open(savepoint_name,"w",0) as output:
        output.write(str(clusters)+"\n")
        output.write(str(labels)+"\n")
        output.write(str(data)+"\n")
        output.write(str(lens))

def from_save(filename):
    with open(filename) as save:
        clusters,labels,data,lens = map(eval,save.readlines())
    for k,v in lens.items():
        lens[k] = int(v)
    pretty_plot(data,labels,clusters,lens)

def main(files,cluster_file):
    matrix,hhrs,lens = load_data(files)
    clusters = read_in_clusters(cluster_file) if cluster_file else {}
    edges,starter = make_a_graph(hhrs)
    columns = go_through_graph(edges,starter)
    fqueue = sort_fams_by_dist(matrix)
    fix_columns(columns,fqueue)
    columns = sorter3(columns)
    data,labels = prepare_data(columns,clusters)
    save_point(data,labels,clusters,lens)
    pretty_plot(data,labels,clusters,lens)

def prepare_data(columns,clusters):
    labels = [x[0] for x in columns[0]]
    data = [[] for x in columns[0]]
    for c in columns:
        for i,x in enumerate(c):
            data[i].append(x[1])
    labels,data = zip(*sorted(zip(labels,data), key=lambda x:clusters.get(x[0],len(clusters)) ))
    return data,labels

def pretty_plot(data,labels,clusters,lens):
    global plot_name
    colors = cm.rainbow(np.linspace(0, 1, len(data)))

    fig = plt.figure()
    ax = plt.axes()

    width = 2
    height = 6
    ax.set_xlim(0, len(data[0])*width)
    ax.set_ylim(0, len(data)*height)

    label_colors = cm.gist_ncar(np.linspace(0, 1, len(set(clusters.values()))+1 ))

    for r,row in enumerate(data):
#        mr = min(_ for _ in row  if _ is not None)
        colors = cm.rainbow(np.linspace(0, 1, lens.get(labels[r],max(row))-2))
#        colors = np.append([0.,0.,0.,1.],np.append(colors,[0.,0.,0.,1.]))
        colors = np.vstack((np.array([0.,0.,0.,1.]),colors,np.array([0.,0.,0.,1.])))
        for i,pos in enumerate(row):
            if pos is not None:
                ax.add_patch(rect((i*width,r*height),width,height-2,color=colors[pos-1] )) #tu bylo pos -mr
                if labels[r] in clusters:
                    ax.add_line(line([i*width,(i+1)*width],[r*height,r*height] , lw=1,color = label_colors[clusters[labels[r]]]))
                    ax.add_line(line([i*width,(i+1)*width],[(r+1)*height-1,(r+1)*height-1] , lw =1, color = label_colors[clusters[labels[r]]]))
    plt.yticks( map(lambda x:x*height+height/2.,range(len(labels))), labels)
#    plt.rc('ytick', labelsize=6)
    ax.tick_params(axis = 'y', which = 'major', labelsize = 8)

    #plt.show()
    plt.tight_layout()

    plt.savefig(plot_name,dpi=300)
#    plt.savefig("merged_alignment_eval_{}.png".format(EVAL_CUTOFF),dpi=300)


def comparable(c1,c2):
    return any(c1[x][1] is not None and c2[x][1] is not None for x in xrange(len(c1)))

def half(it):
    lit = list(it)
    return sum(lit)*1./len(lit)>=.5

def first_lower(c1,c2):
    return all(c1[x][1] < c2[x][1]  for x in xrange(len(c1)) if c1[x][1] is not None and c2[x][1] is not None)


def sorter3(columns):
    columns = sorted(columns, key=lambda x: len([_ for _ in x if _[1] is not None]))
    csorted = []
    wstawiony=0

    while columns:
        c = columns.pop()
        csorted.append(c)
        no_change = 0
        while columns:
            col = columns.pop()
            found_comp = None
            for i,s in enumerate(csorted):
                if comparable(col,s):
                    found_comp = i
                    no_change = 0
                    if first_lower(col,s):
                        csorted.insert(i,col)
                        wstawiony+=1
                        break
            else:
                if found_comp is not None: 
                        csorted.insert(found_comp+1,col)
                        wstawiony+=1
                else:
                    columns.append(col)
                    no_change += 1
                if len(columns)+1 == no_change:
                    break
    return csorted

def fix_columns(columns, fams):
    for c,col in enumerate(columns):
        present = [_[0] for _ in col]
        for f in fams:
            if f not in present:
                col.append((f,None))
        columns[c] = sorted(col)

def make_a_graph(hhrs):
    edges = []
    popular = {}
    evals = []
    for q in hhrs:
        for t in q.hits:
          if t.eval < EVAL_CUTOFF:
            for qc,tc in t.match:
                if qc!=None and tc!=None:
                    edges.append(((q.query,qc),(t.target,tc)))
                    popular[(q.query,qc)] = popular.get((q.query,qc),0) + 1
                    popular[(t.target,tc)] = popular.get((t.target,tc),0) + 1
                    evals.append(len(t.match))
    evals,edges = map(list,zip(*sorted(zip(evals,edges))))
    fix_graph(edges,evals)
    return list(edges),sorted(popular.keys(),key=lambda x: popular[x], reverse=True)[0]

def fix_graph(edges,evals):
    delete = []
    for i,e in enumerate(edges):
        for j,f in enumerate(edges[i+1:]):
            if len(set(e+f)) == 3 and len(set([_[0] for _ in e+f])) == 2:
                if evals[i] > evals[j+i]:
                    delete.append(j+i)
                else:
                    delete.append(i)
    delete.sort()
    while delete:
        d = delete.pop()
        edges.pop(d)
        evals.pop(d)


def neigh(edges,vert):
    new = []
    delete = []
    for i,e in enumerate(edges):
        if vert in e:
            nv = e[1] if e[0]==vert else e[0]
            new.append(nv)
            delete.append(i)
    while delete:
        edges.pop(delete.pop())
    return new

def in_column(column,what):
    return what[0] in [_[0] for _ in column]

def bfs(edges,start,column):
    queue = [start]
    while queue:
        cur = queue.pop()
        column.append(cur)
        ns = neigh(edges,cur)
        for n in ns:
            if not in_column(column,n) and not in_column(queue,n):
                queue.append(n)

def go_through_graph(edges,starter):
    columns = []
    while edges:
        column = []
        bfs(edges,starter,column)
#        get_all(edges,starter,column)
        columns.append(column)
        if edges:
            starter = edges[0][0]
    return columns
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge hhpred alignments')
    parser.add_argument('hhrs', metavar='HHR', type=str, nargs='*',
                        help='.hhr files to be merged')
    parser.add_argument('--hhr_dir', type=str, default="", help='Directory from which all .hhr will be taken')
    parser.add_argument('--plot_from_save', nargs=1, type=str, default="", help='Savepoint from a previous run')
    parser.add_argument('--cluster_file', nargs=1, type=str, default="", help='Cluster defining file: each cluster in one line')
    parser.add_argument('--save_name', nargs=1, type=str, default="savepoint_{}.txt".format(time.strftime("%Y%m%d-%H%M%S")),
                         help='Custom savepoint name')
    parser.add_argument('--plot_name', nargs=1, type=str, default="merged_alignment_{}.png".format(time.strftime("%Y%m%d-%H%M%S")),
                         help='Custom plot name')
    args = parser.parse_args()
    
    savepoint_name = args.save_name
    if type(savepoint_name) is not str:
        savepoint_name = savepoint_name[0]

    plot_name = args.plot_name
    if type(plot_name) is not str:
        plot_name = plot_name[0]
    
    just_plot = args.plot_from_save
    if just_plot:
        just_plot = just_plot[0]

    cluster_file = args.cluster_file
    if cluster_file:
        cluster_file = cluster_file[0]

    if just_plot:
        from_save(just_plot)
    else:
        if not args.hhrs and not args.hhr_dir:
            exit("No .hhr files indicated")
        hhrs = args.hhrs
        d_hhrs = glob.glob("{}/*.hhr".format(args.hhr_dir))
        hhrs = list(set(hhrs + d_hhrs))
        main(hhrs,cluster_file)

