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

ctg2id = {}
id2ctg = {}

def mapctg(ctgid):
    if ctgid in id2ctg.keys():
        return id2ctg[ctgid]
    else:
        return -1
for prot,ctg,pc in zip(pcdf['protein_id'],pcdf['contig_id'],pcdf['cluster']):
    if pc == -1:
        continue
    if ctg in ctg2id.keys():
        ctgid = ctg2id[ctg]
        prot2ctg[prot] = ctgid
    else:
        nextid = len(ctg2id)
        ctg2id[ctg] = nextid
        id2ctg[nextid] = ctg
        prot2ctg[prot] = nextid

def protmap(prot):
    global prot2ctg
    if prot in prot2ctg.keys():
        return prot2ctg[prot]
    else:
        return -1

   
tabdf = pd.read_csv(args.tabfile,header=0,sep='\t')
ctg1 = tabdf['prot1'].map(protmap)
ctg2 = tabdf['prot2'].map(protmap)
tabdf['ctg1'] = ctg1
tabdf['ctg2'] = ctg2


#newdf = tabdf.groupby(['ctg1','ctg2']).agg({'e':['min'],'score':[np.sum]}).reset_index()
newdf = tabdf.groupby(['ctg1','ctg2']).agg({'hsp':[np.sum],'score':[np.sum]}).reset_index()
newdf.columns = ['ctg1','ctg2','hsp','score']

countdf = tabdf.groupby(['ctg1','ctg2']).size().reset_index(name="count")
newdf['count']  = countdf['count']
newdf.to_csv('tmp/aggrscore.csv',index=False,sep='\t')

ctgname1 = newdf['ctg1'].map(mapctg)
ctgname2 = newdf['ctg2'].map(mapctg)
newdf['ctgname1'] = ctgname1
newdf['ctgname2'] = ctgname2
newdf.to_csv('tmp/aggrdetailscore.csv',index=False,sep='\t')

id_ctg = list(id2ctg.items())
ctgdf = pd.DataFrame(id_ctg,columns=['id','ctg'])
ctgdf.to_csv('tmp/id2ctg',sep='\t',index=False)

