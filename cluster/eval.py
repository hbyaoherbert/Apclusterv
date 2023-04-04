import numpy as np
import pandas as pd
from scipy.special import comb
import argparse

def rand_index_score(clusters, classes):
    tp_plus_fp = comb(np.bincount(clusters), 2).sum()
    tp_plus_fn = comb(np.bincount(classes), 2).sum()
    A = np.c_[(clusters, classes)]
    tp = sum(comb(np.bincount(A[A[:, 0] == i, 1]), 2).sum()
             for i in set(clusters))
    fp = tp_plus_fp - tp
    fn = tp_plus_fn - tp
    tn = comb(len(A), 2) - tp - fp - fn

    print(tp,tn,fp,fn)
    return (tp + tn) / (tp + fp + fn + tn)


def c1tores(resfile,tax):
    taxdf = pd.read_csv(args.tax,sep='\t',header=0)
    ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))

    ctg2id = {}
    genus2id = {}
    ctg2genus = {}
    classes = []

    
    for ctgkey,genus in  zip(ctgkey,taxdf['genus']):
        if not genus in genus2id.keys():
            thisid = len(genus2id)
            genus2id[genus] = thisid

        thisid = len(ctg2id)
        ctg2id[ctgkey] = thisid
       
        ctg2genus[ctgkey] = genus2id[genus]
        classes.append(genus2id[genus])


    df = pd.read_csv(args.resfile,header=0,sep=',')
  
    clusters = np.ones(len(classes),dtype=int) * 100000

    for idx,row in df.iterrows():
        for ctg in row['Members'].split(','):
            if not ctg in ctg2genus.keys():
                continue
            


            clusters[ctg2id[ctg]] = idx
    
    return clusters,classes

def c1res(resfile,tax):
    cldf = pd.read_csv(resfile,header=0,sep=',')
    ctgset = set()
    for members in cldf['Members']:
        for ctg in members.split(','):
            ctgset.add(ctg)

    taxdf = pd.read_csv(args.tax,sep='\t',header=0)
    ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))

    ctg2id = {}
    genus2id = {}
    ctg2genus = {}
    classes = []

    
    for ctgkey,genus in  zip(ctgkey,taxdf['genus']):
        if not ctgkey in ctgset:
            continue
        if not genus in genus2id.keys():
            thisid = len(genus2id)
            genus2id[genus] = thisid

        thisid = len(ctg2id)
        ctg2id[ctgkey] = thisid
       
        ctg2genus[ctgkey] = genus2id[genus]
        classes.append(genus2id[genus])


    df = pd.read_csv(args.resfile,header=0,sep=',')
  
    clusters = np.ones(len(classes),dtype=int) * 100000

    for idx,row in df.iterrows():
        for ctg in row['Members'].split(','):
            if not ctg in ctg2genus.keys():
                continue
            


            clusters[ctg2id[ctg]] = idx
    c2size = {}
    for i in clusters:
        if i in c2size.keys():
            c2size[i] += 1
        else:
            c2size[i] = 1
    sensitivity = 0
    for i,num in c2size.items():
        if num>=2:
            sensitivity += num
    print("sensitivity:",sensitivity)
    return clusters,classes

parser = argparse.ArgumentParser()
parser.add_argument("resfile",type=str)
parser.add_argument("tax",type=str)
args = parser.parse_args()

clusters,classes = c1tores(args.resfile,args.tax)

print(rand_index_score(clusters,classes) )       

    

