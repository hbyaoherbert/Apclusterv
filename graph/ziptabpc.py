import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("pcfile",type=str)
parser.add_argument("tabfile",type=str)
args = parser.parse_args()


pcdf = pd.read_csv(args.pcfile,sep=',',header=0)
pcdf = pcdf.fillna(-1)
prot2pc = {}
for prot,contig,pc in zip(pcdf['protein_id'],pcdf['contig_id'],pcdf['cluster']):
    if pc == -1:
        continue
    pc = int(pc.split('_')[1])
    prot2pc[prot] = pc

def protmap(prot):
    global prot2pc
    if prot in prot2pc.keys():
        return prot2pc[prot]
    else:
        return -1
       
tabdf = pd.read_csv(args.tabfile,header=0,sep='\t')
tabdf.columns = ['prot1','prot2','hsp','e','score']
cluster1 = tabdf['prot1'].map(protmap)
cluster2 = tabdf['prot2'].map(protmap)
tabdf['pc1'] = cluster1
tabdf['pc2'] = cluster2
clarray = np.array(tabdf[['pc1','pc2']].values,dtype=int)
farray = tabdf.values
mask = np.all( (clarray[:,0:1]!=-1) & (clarray[:,1:2]==clarray[:,0:1]) ,axis=1)
clarray = clarray[mask]

farray = farray[mask]


outdf = pd.DataFrame(farray,columns = ['prot1','prot2','hsp','e','score','pc1','pc2'])

outdf.to_csv('graph/pc.tab',sep='\t',index=False)

