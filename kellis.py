from collections import defaultdict
from scipy.stats import binom_test
import os
import numpy as np

def parser(file_name):
    genes = defaultdict(str)
    Y_names = {}
    with open(file_name, 'r') as handle:
        for line in handle:
            if not (line[0] == "K" or line[0] == "Y"):
                continue

            try:
                name, seq = line.split()
            except:
                print line
                return genes
            if name[0] == 'K':
                name = "K"
            else:
                if name not in Y_names:
                    Y_names[name] = "Y%d" % len(Y_names)
                name = Y_names[name]

            genes[name] += seq

    return genes


def get_diss(genes):
    """get_diss
    
    Parameters
    ----------
    genes : list
      gene is a list with three strings, corresponding to the sequence
    
    Returns
    -------
    """
    if not len(genes) == 3:
        print genes.keys()
    assert len(genes) == 3
    gv = genes.values()
    if not len(gv[0]) == len(gv[1]) == len(gv[2]):
        print len(gv[0]) , len(gv[1]) , len(gv[2])
    n = min(len(gv[0]) , len(gv[1]) , len(gv[2]))

    nucleotides = 'ACGTacgt'

    sY = sK = sZ = 0
    for i in xrange(n):
        y, z, k = genes['Y0'][i], genes['Y1'][i], genes['K'][i]

        if y not in nucleotides:
            continue
        if z not in nucleotides:
            continue
        if k not in nucleotides:
            continue

        if y == z != k:
            sK += 1
        elif y == k != z:
            sZ += 1
        elif k == z != y:
            sY += 1

    return sK, sZ, sY
        
    
def do_binom_test(s):
    return binom_test(s[:2])

def do_binom_test2(s):
    return binom_test((s[0], s[2]))

def do_accel_test(s):
    return (s[1] + s[2]) / 2. / s[0]


def run_all_genes(directory="."):
    files = os.listdir(directory)
    alignments = [directory + "/" + f for f in files if f.endswith("t.txt")]
    
    alignments2 = []
    g1 = get_76_fast()
    for al in alignments:
        for gene in g1:
            if gene in al:
                alignments2.append(al)
                continue

    if False:
        svals, pvals, genes = [], [], []
        for i, a in enumerate(alignments):
            print i, a
            g = parser(a)
            s = get_diss(g)
            p = do_binom_test(s)
            pvals.append(p)

    genes = [parser(g) for g in alignments2]
    svals = [get_diss(g) for g in genes]
    avals = [do_accel_test(s) for s in svals]
    bvals = [do_binom_test2(s) for s in svals]
    pvals = [do_binom_test(s) for s in svals]
    return np.vstack((avals, bvals, pvals, np.transpose(svals)))

def get_76_fast():
    import pandas as pd
    x=pd.read_excel("../../S9_Trees/Duplicated_Pairs.xls")
    y = pd.read_excel("../../S9_Trees/Nucleotide_Divergence2.xls")
    to_keep = x['Unnamed: 11']=='76 fast' 
    x2 = x[to_keep]

    z1 = pd.merge(x2,y,how='inner',
                  left_on=('gene1'),
                  right_on=('name1'))
    z2 = pd.merge(x2,y,how='inner',
                  left_on=('gene1'),
                  right_on=('gene1'))
    z3 = pd.merge(x2,y,how='inner',
                  left_on=('gene1'),
                  right_on=('gene2'))
    z4 = pd.merge(x2,y,how='inner',
                  left_on=('gene1'),
                  right_on=('name2'))
    g1 = list(z1['name1'])
    g1.extend(list(z2['name1']))
    g1.extend(list(z3['name1']))
    g1.extend(list(z4['name1']))
    return np.unique(g1)


def shitty_sims(s, n=100000):
    from numpy.random import binomial as rbinom
    ss = s[0] + s[1]
    ss = np.array(ss, dtype=np.int)

    res = np.empty(n)
    for i in xrange(n):
        r = rbinom(ss, 0.5)
        res[i] = np.sum(np.min((r, ss-r), 0) < np.min(s, 0))

    return res



if __name__ == "__main__":
    p = run_all_genes('/data/kellis/S8_DupGenes/')
    np.savetxt("s.txt", p)


