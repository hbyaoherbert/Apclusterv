import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("pcfile",type=str)
parser.add_argument("tabfile",type=str)
parser.add_argument("aggr",type=str)
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

def prot2key(row):
    return (row['prot1'],row['prot2'])



   
tabdf = pd.read_csv(args.tabfile,header=0,sep='\t')
ctg1 = tabdf['prot1'].map(protmap)
ctg2 = tabdf['prot2'].map(protmap)
tabdf['ctg1'] = ctg1
tabdf['ctg2'] = ctg2


#newdf = tabdf.groupby(['ctg1','ctg2']).agg({'e':['min'],'score':[np.sum]}).reset_index()
#newdf = tabdf.groupby(['ctg1','ctg2']).agg({'hsp':[np.sum],'score':[np.sum]}).reset_index()
#newdf.columns = ['ctg1','ctg2','hsp','score']

oceandf1 = tabdf.groupby(['prot1','ctg2']).size().reset_index()
oceandf1.columns = ['prot','ctg','count']
oceandf1 = oceandf1[oceandf1['count']>1]

oceandf2 = tabdf.groupby(['prot2','ctg1']).size().reset_index()
oceandf2.columns = ['prot','ctg','count']
oceandf2 = oceandf2[oceandf2['count']>1]


poi1 = tabdf[tabdf.prot1.isin(list(oceandf1['prot']))].reset_index()
poi1 = poi1[poi1.ctg2.isin(list(oceandf1['ctg']))].reset_index()
poi2 = tabdf[tabdf.prot2.isin(list(oceandf2['prot']))].reset_index()
poi2 = poi2[poi2.ctg1.isin(list(oceandf2['ctg']))].reset_index()

alignments1 =  set(poi1.apply(prot2key,axis=1))
alignments2 =  set(poi2.apply(prot2key,axis=1))
alignments = alignments1.union(alignments2)
print(len(alignments))


id_ctg = list(id2ctg.items())
#ctgdf = pd.DataFrame(id_ctg,columns=['id','ctg'])
#ctgdf.to_csv('cltest/id2ctg',sep='\t',index=False)


diamonddf = pd.read_csv('merged.self-diamond.tab',sep='\t',header=None)
diamonddf.columns = ['prot1','prot2','pid','alnlen','mismatch','gap','qstart','qend','tstart','tend','evalue','bitscore']
duplicate = diamonddf[diamonddf.apply(prot2key,axis=1).isin(alignments)].reset_index()

aln2vec = {}
aln2hsp = {}
data = []
for prot1,prot2,qstart,qend,bitscore in zip(duplicate['prot1'],duplicate['prot2'],duplicate['qstart'],duplicate['qend'],duplicate['bitscore']):
    ctg2 = prot2ctg[prot2]
    alnkey = (prot1,ctg2)
    if alnkey in aln2hsp.keys():
        aln2hsp[alnkey].append([qstart,qend,bitscore/(qend-qstart+1)])
    else:
        aln2hsp[alnkey] = [[qstart,qend,bitscore/(qend-qstart+1)]]
for alnkey in aln2hsp.keys():
    maxlen = 0
    for qstart,qend,score in aln2hsp[alnkey]:
        maxlen = max(maxlen,qend)
    scorevec = np.zeros(maxlen+1)
    for qstart,qend,score in aln2hsp[alnkey]:
        hsplen = qend-qstart+1
        vec = [score for i in range(hsplen)]
        scorevec[qstart:qend+1] = np.maximum(scorevec[qstart:qend+1],vec)
    aggrscore = np.sum(scorevec)
    prot,ctg = alnkey
    data.append([prot,ctg,aggrscore])

alndf = pd.DataFrame(data,columns = ['prot1','ctg2','score'])
alndf['ctg1'] = alndf['prot1'].map(prot2ctg)

alndf = alndf.groupby(['ctg1','ctg2']).agg({'score':[np.sum]}).reset_index()
alndf.columns = ['ctg1','ctg2','score']
#alndf.to_csv('oceannew.csv',sep='\t',index=False)



ctg1 = duplicate['prot1'].map(protmap)
ctg2 = duplicate['prot2'].map(protmap)
duplicate['ctg1'] = ctg1
duplicate['ctg2'] = ctg2
flatscore = duplicate.groupby(['ctg1','ctg2']).agg({'bitscore':[np.sum]}).reset_index()
flatscore.columns = ['ctg1','ctg2','score']
#flatscore.to_csv('oceanflat',sep='\t',index=False)

alndf['score'] = alndf['score'] - flatscore['score']
#alndf.to_csv('cltest/dupinc.csv',sep='\t',index=False)

aln2score = {}



aggrdf = pd.read_csv(args.aggr,sep='\t',header=0)
for ctg1,ctg2,score in zip(aggrdf['ctg1'],aggrdf['ctg2'],aggrdf['score']):
    ctgname1 = id2ctg[ctg1]
    ctgname2 = id2ctg[ctg2]
    aln2score[(ctgname1,ctgname2)] = score

for ctg1,ctg2,score in zip(alndf['ctg1'],alndf['ctg2'],alndf['score']):
    ctgname1 = id2ctg[ctg1]
    ctgname2 = id2ctg[ctg2]
    aln2score[(ctgname1,ctgname2)] += score


data = []
for ctg_pair,score in aln2score.items():
    ctg1 = ctg_pair[0]
    ctg2 = ctg_pair[1]
    data.append([ctg1,ctg2,score])
    '''
    if aln2score[(ctg1,ctg1)] > aln2score[(ctg2,ctg2)]:
        continue    
    
    data.append([ctg1,ctg2,aln2score[ctg_pair]/aln2score[(ctg1,ctg1)] ])
    '''

normdf = pd.DataFrame(data,columns=['ctg1','ctg2','score'])
normdf.to_csv("tmp/aggrnorm.csv",sep='\t',index=False)
