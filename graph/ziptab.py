import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("pcfile",type=str)
parser.add_argument("clres",type=str)
parser.add_argument("tabfile",type=str)
args = parser.parse_args()

clres = pd.read_csv(args.clres,sep=',',header=0)
cluster = 0
annotated = set()
ctg2cl = {}
for members in list(clres['Members']):
    for member in members.split(' '):
        ctg2cl[member] = cluster
    cluster += 1


pcdf = pd.read_csv(args.pcfile,sep=',',header=0)
pcdf = pcdf.fillna(-1)
prot2cl = {}
for prot,contig,pc in zip(pcdf['protein_id'],pcdf['contig_id'],pcdf['cluster']):
    if cluster == -1:
        continue
    ctg = contig.replace(' ','~')
    if not ctg in ctg2cl.keys():
        continue
    prot2cl[prot] = ctg2cl[ctg]

def protmap(prot):
    global prot2cl
    if prot in prot2cl.keys():
        return prot2cl[prot]
    else:
        return -1
       
tabdf = pd.read_csv(args.tabfile,header=None,sep=' ')
tabdf.columns = ['prot1','prot2','e']
cluster1 = tabdf['prot1'].map(protmap)
cluster2 = tabdf['prot2'].map(protmap)
tabdf['cluster1'] = cluster1
tabdf['cluster2'] = cluster2
clarray = np.array(tabdf[['cluster1','cluster2']].values,dtype=int)
farray = tabdf.values
mask = np.all( (clarray[:,0:1]==-1) | (clarray[:,1:2]==-1),axis=1)
clarray = clarray[mask]
farray = farray[mask]


missmask = np.all((clarray[:,0:1]==-1) & (clarray[:,1:2]==-1),axis=1)
novel = farray[missmask]
reach = farray[np.logical_not(missmask)]

noveldf = pd.DataFrame(novel,columns = ['prot1','prot2','e','cluster1','cluster2'])
noveldf.to_csv('script/clnovel.tab',sep='\t',index=False) 
reachdf = pd.DataFrame(reach,columns = ['prot1','prot2','e','cluster1','cluster2'])
reachdf.to_csv('script/clreach.tab',sep='\t',index=False)

