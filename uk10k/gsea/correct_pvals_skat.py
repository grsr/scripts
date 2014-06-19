import sys
from pandas import *

if __name__ == '__main__':
    # 1       861341  866445  SAMD11.etc.1    2.466915e-01    1       861263  866445  AL645608.1      ENSG00000268179 -
    df = read_table(sys.argv[1], names=['chr', 'start', 'end', 'id', 'p', 'g_chr', 'g_start', 'g_end', 'g_name', 'g_id', 'g_strand'])
    genes = {}

    for gene in df.g_name.unique():
        print 'processing',gene
        gene_vars = df[df.g_name==gene]
        pvals = gene_vars['p']       
        n_vars = len(pvals)
        if n_vars == 0:
            print gene,'has no variants?'
            continue
        p_raw = pvals.min()
        p_adj = 1 - ((1 - p_raw)**n_vars)
        gene_var1 = gene_vars.iloc[0]
        gene_length = gene_var1['g_end'] - gene_var1['g_start'] + 1
        genes[gene] = {'p_raw': p_raw, 'p_adj': p_adj, 'n_variants': n_vars, 'g_length': gene_length}

    res = DataFrame(genes).T
    res.to_csv(sys.argv[2])

