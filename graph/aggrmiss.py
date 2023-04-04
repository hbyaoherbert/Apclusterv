import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("pcfile",type=str)
parser.add_argument("tabfile",type=str)
args = parser.parse_args()

pcdf = pd.read_csv(args.pcfile,sep=',',header=0)
pcdf = pcdf.fillna(-1)
prot2ctg = {}
for prot,contig,cluster in zip(pcdf['protein_id'],pcdf['contig_id'],pcdf['cluster']):
    if cluster == -1:
        continue
    ctg = contig.replace(' ','~')
    prot2ctg[prot] = ctg


def protmap(prot):
    if prot in prot2ctg.keys():
        return prot2ctg[prot]
    else:
        return -1

tabdf = pd.read_csv(args.tabfile,sep='\t',header=0)

keys = []
for idx,row in tabdf.iterrows():
    ctg1 = protmap(row['prot1'])
    ctg2 = protmap(row['prot2'])
    if ctg1 == -1 or ctg2 == -1:
        continue
    clkey = ctg1+'#'+ctg2
    keys.append(clkey)
tabdf['clkey'] = keys

outdf = tabdf.groupby('clkey')['e'].agg('min').reset_index();
outdf.to_csv('script/edgenovel.tab',header=None,index=False)

