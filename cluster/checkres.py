import argparse
import pandas as pd
import numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("resfile",type=str)
parser.add_argument("between",type=str)
#parser.add_argument("tax",type=str)

args = parser.parse_args()

betweendf = pd.read_csv(args.between,sep='\t',header=0)
node2between = dict(zip(betweendf['node'],betweendf['between']))

taxdf = pd.read_csv('tax.tab',sep='\t',header=0)
ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))
ctg2genus = dict(zip(ctgkey,taxdf['genus']))

cldf = pd.read_csv(args.resfile,sep=',',header=0)
generas = []
for idx,row in cldf.iterrows():
    genera = set()
    betweens = []
    for member in row['Members'].split(','):
        between = 0
        if member in node2between.keys():
            between = node2between[member]
        betweens.append(between)
        if member in ctg2genus.keys():
            genera.add(ctg2genus[member])
        
    if row['Size']>=20:
        print(row['Members'],genera,np.max(betweens),np.median(betweens))
    generas.append(str(list(genera)))

cldf['genus'] = generas

cldf.to_csv('viral_cluster_res',sep=',')
