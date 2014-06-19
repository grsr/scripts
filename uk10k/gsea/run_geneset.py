from pandas import *
import numpy as np
import sys
import os

if __name__ == '__main__':
    gene_pvals = read_csv(sys.argv[1], index_col=0)
    gene_sets = open(sys.argv[2], 'r')
    out_dir = sys.argv[3]
    os.mkdir(out_dir)
    for line in gene_sets:
        gene_set = line.strip().split('\t')
        name = gene_set[0]
        url = gene_set[1]
        genes = set(gene_set[2:])
        gene_pvals['present'] = gene_pvals.apply(lambda r: 1 if r.name in genes else 0, axis=1)
        gene_pvals[['p_adj', 'present']].to_csv(out_dir+'/'+name+'.txt', index=False, header=None, sep=" ")

