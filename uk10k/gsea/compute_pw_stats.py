from pandas import *
import numpy as np
import sys
from glob import glob
from os.path import splitext, basename, getsize

if __name__ == '__main__':
    res = {}
    all_nes_p = Series()
    all_nes_s = Series()
    #max_nes_p = Series()
    max_nes_p = None
    for f in glob(sys.argv[1]+'*.out'):
        pw = splitext(basename(f))[0]
        if getsize(f) == 0:
            print "no results for pathway", pw
            continue
        es = read_table(f, header=None, squeeze=True)
        es_s = es[0]
        es_p = es[1:]
        es_p_mean = es_p.mean()
        es_p_std = es_p.std()
        nes_s = (es_s - es_p_mean) / es_p_std
        nes_p = (es_p - es_p_mean) / es_p_std
        pval = len(nes_p[nes_p >= nes_s]) / float(len(nes_p))  
        res[pw] = {'pval': pval, 'nes_s': nes_s}
        all_nes_p = all_nes_p.append(nes_p)
        all_nes_s = all_nes_s.append(Series([nes_s]))
        #max_nes_p = max_nes_p.append(Series([nes_p.max()]))
        if max_nes_p is None:
            max_nes_p = nes_p.values
        else:
            for i, m in enumerate(nes_p):
                if max_nes_p[i] < m:
                    max_nes_p[i] = m

    for pw in res.keys():
        nes_star = res[pw]['nes_s']
        top = len(all_nes_p[all_nes_p >= nes_star]) / float(len(all_nes_p))
        bottom = len(all_nes_s[all_nes_s >= nes_star]) / float(len(all_nes_s))
        fdr = top / bottom
        fdr = 1 if fdr > 1 else fdr
        res[pw]['fdr'] = fdr
        #fwer = len(all_nes_p[all_nes_p >= nes_star]) / float(len(all_nes_p))
        fwer = len(max_nes_p[max_nes_p >= nes_star]) / float(len(max_nes_p))
        res[pw]['fwer'] = fwer

    df = DataFrame.from_dict(res, orient='index')
    df.index.name = 'pathway'
    
    df.sort('pval', inplace=True)

    df.to_csv(sys.argv[2])

