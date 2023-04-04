from sklearn.cluster import AffinityPropagation
from sklearn.cluster import *
import pandas as pd
import argparse
import numpy as np
import adj_ap
import networkx as nx
import json

parser = argparse.ArgumentParser()
parser.add_argument('c1abs')
parser.add_argument('cluster')
parser.add_argument('scorefile')
parser.add_argument('majorfile')


configfile = open('config','r')
var2value = {}
var2value['k'] = 1

for line in configfile:
    line = line.strip()
    info = line.split('=')
    if len(info) != 2:
        continue
    var = info[0].strip()
    value = float(info[1].strip())
    if var in var2value.keys():
        var2value[var] = value

args = parser.parse_args()


scoredf = pd.read_csv(args.scorefile,sep='\t',header=0)
pairscore = {}
for ctg1,ctg2,score in zip(scoredf['ctg1'],scoredf['ctg2'],scoredf['score']):
    ctg1 = ctg1.replace(' ','~')
    ctg2 = ctg2.replace(' ','~')

    pairscore[(ctg1,ctg2)] = score  

complements = []
c2size = {}
c2sim = {}
c2members = {}

def sumcluster(cluster):
    cldf = pd.read_csv(cluster,sep='\t',header=None)
    cldf.columns = ['i','rep','contig']
    rep2ctg = {}
    ctg2rep = {}
    for idx,row in cldf.iterrows():
        ctg2rep[row['contig']] = row['rep']
        if row['rep'] in rep2ctg.keys():
            rep2ctg[row['rep']].append(row['contig'])
        else:
            rep2ctg[row['rep']] = [row['contig']]
    print(len(rep2ctg))
    return rep2ctg,ctg2rep

def separate(rep2ctg,ctg2rep,c1abs):
    rep2edges = {}
    for rep in rep2ctg.keys():
        rep2edges[rep] = []

    edge1 = pd.read_csv(c1abs,sep='\t',header=None)
    edge1.columns = ['ctg1','ctg2','hsp','score']
    

    for ctg1,ctg2,hsp,score in zip(edge1['ctg1'],edge1['ctg2'],edge1['hsp'],edge1['score']):
        if not ctg1 in ctg2rep.keys():
            continue
        if not ctg2 in ctg2rep.keys():
            continue

        if ctg2rep[ctg1] != ctg2rep[ctg2]:
            continue
        if pairscore[(ctg1,ctg1)] < pairscore[(ctg2,ctg2)]:
            shortctg = ctg1
            longctg = ctg2
        else:
            shortctg = ctg2
            longctg = ctg1
            
        if not (shortctg,longctg) in pairscore.keys():
            continue
        similarity = pairscore[(shortctg,longctg)]/pairscore[(shortctg,shortctg)]
        rep2edges[ctg2rep[ctg1]].append([ctg1,ctg2,hsp,similarity])

    
    
    return rep2edges

def printresult(majorfile):

    with open('tmp/c2members','w') as f1:
        json.dump(c2members,f1)
    with open('tmp/c2size','w') as f2:
        json.dump(c2size,f2)
    with open('tmp/c2sim','w') as f3:
        json.dump(c2sim,f3)

    apexamplars = []
    c2membersf = {}
    for line in open(majorfile):
        line = line.strip()
        info = line.split('\t')
        ctgs = info[0].split(',')
        examplars = info[1].split(',')
        
        for ctg,examplar in zip(ctgs,examplars):
            apexamplars.append(examplar)
            if examplar in c2membersf.keys():
                c2membersf[examplar].append(ctg)
            else:
                c2membersf[examplar] = [ctg]

    for examplar,members in c2members.items():
        c2membersf[examplar] = members

    resdata = []

    for examplar,members in c2membersf.items():
        memberstr = ','.join(members)
        resdata.append(memberstr)
    resdf = pd.DataFrame(resdata,columns=['Members'])
    resdf.to_csv('testres.csv',sep=',',index=False)

    complementdf = pd.DataFrame(complements,columns = ['ctg'])
    complementdf.to_csv('tmp/complements',index=False)

def runcluster():
    rep2ctg,ctg2rep = sumcluster(args.cluster)
    rep2edges = separate(rep2ctg,ctg2rep,args.c1abs)

    for rep,edges in rep2edges.items():

        
        ctg2id,id2ctg = construct(edges)
        if len(ctg2id) < 5:
            continue
        #outfile = open('sres/{}.res'.format(rep),'w')
        aptest(edges,ctg2id,id2ctg)
        #outfile.close()
    printresult(args.majorfile)
    #attach(args.scorefile,args.majorfile,0.15)

def construct(edges):
    ctg2id = {}
    id2ctg = {}

    for ctg1,ctg2,hsp,score in edges:
        if not ctg1 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg1] = thisid
            id2ctg[thisid] = ctg1
            
        if not ctg2 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg2] = thisid
            id2ctg[thisid] = ctg2 
    return ctg2id,id2ctg
'''
def construct(crossgenus,ingenus,rep2edges):


    cldf = pd.read_csv(cluster,sep='\t',header=None)
    cldf.columns = ['i','rep','contig']
    ctgset = set(cldf['contig'])
    edges = []
    edge1 = pd.read_csv(crossgenus,sep='\t',header=None)
    edge1.columns = ['ctg1','ctg2','e','score','count','genus1','genus2']
    
    for ctg1,ctg2,e,score in zip(edge1['ctg1'],edge1['ctg2'],edge1['e'],edge1['score']):
        if not ctg1 in ctgset:
            continue
        if not ctg2 in ctgset:
            continue

        if not ctg1 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg1] = thisid
            id2ctg[thisid] = ctg1
            
        if not ctg2 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg2] = thisid
            id2ctg[thisid] = ctg2 
        
        edges.append([ctg1,ctg2,e,score])
    
    edge2 = pd.read_csv(ingenus,sep='\t',header=None)
    edge2.columns = ['ctg1','ctg2','e','score','count','genus1','genus2']
    for ctg1,ctg2,e,score in zip(edge2['ctg1'],edge2['ctg2'],edge2['e'],edge2['score']):
        if not ctg1 in ctgset:
            continue
        if not ctg2 in ctgset:
            continue

        if not ctg1 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg1] = thisid
            id2ctg[thisid] = ctg1
            
        if not ctg2 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg2] = thisid
            id2ctg[thisid] = ctg2 
            
        edges.append([ctg1,ctg2,e,score])
    #print(edges)    
    return edges
'''
def aptest(edges,ctg2id,id2ctg):
    
    N = len(id2ctg)
    mtx = np.zeros((N,N))

    simg =  nx.Graph()
    for ctg1,ctg2,hsp,score in edges:
        ctgid1 = ctg2id[ctg1]
        ctgid2 = ctg2id[ctg2]
        
        sim = score
   
        #simg.add_edge(ctg1,ctg2,weight=sim)
        mtx[ctgid1,ctgid2] = sim
        mtx[ctgid2,ctgid1] = sim
    #print(mtx,mtx.shape,np.max(mtx),np.mean(mtx))
    #mtx = mtx/np.max(mtx)
    for i in range(N-1):
        for j in range(i+1,N):
            if mtx[i,j] == 0:
                continue
            simg.add_edge(i,j,weight=mtx[i,j])
        

    #print(simg.edges(data=True))
    clcoef = nx.clustering(simg,weight='weight')
    between = nx.betweenness_centrality(simg)
    
    preference = np.zeros(N)
    #preference = np.median(mtx)
    #print(clcoef)
    pstep = np.zeros(N)


    for node,coef in clcoef.items():
        #print(node,coef,between[node])
        preference[node] = coef - var2value['k']*between[node]
        pstep[node] = 0.01
    preference = np.clip(preference,0,1)
    #preference = np.mean(mtx,axis=0)
    #pstep = (np.max(mtx) - np.min(mtx) )/100
    #print(ctg2id)
    #cluster_centers_indices, labels = affinity_propagation(mtx, max_iter=1000,preference=preference, random_state=39)
    print(preference)
    cluster_centers_indices, labels = adj_ap.adj_affinity_propagation(mtx,damping=0.2,pstep=pstep, max_iter=2000,preference=preference, random_state=39)
    centerstr  = ''

    print(cluster_centers_indices)
    for center in cluster_centers_indices:
        centerstr+=id2ctg[center]+','
        complements.append(id2ctg[center])
        c2size[id2ctg[center]] = 0
        c2sim[id2ctg[center]] = -1
        c2members[id2ctg[center]] = []
    print(labels)


    idx = 0
    for label in labels:
        examplar = id2ctg[cluster_centers_indices[label]]
        ctg = id2ctg[idx]
        c2members[examplar].append(ctg)
        if cluster_centers_indices[label] != idx:
            c2sim[examplar] = max(c2sim[examplar],mtx[cluster_centers_indices[label],idx])
        c2size[examplar] += 1
        idx += 1
        
    centerstr = centerstr.rstrip(',')
    #outfile.write("{}\n{}\n{}\n".format(centerstr,N,np.mean(preference)))
    '''
    for i in range(len(labels)):
        betweeness = 0
        if id2ctg[i] in between.keys():
            betweeness = between[id2ctg[i]]
        print(labels[i],id2ctg[i],betweeness)
    '''


def minsim(simlist,pairscore):
    res = 10
    if len(simlist) == 1:
        return 1
    for query,num in simlist:
        for target,num in simlist:
            if query == target:
                continue
            if pairscore[(query,query)] < pairscore[(target,target)]:
                shortctg = query
                longctg = target
            else:
                shortctg = target
                longctg = query
            
            if not (shortctg,longctg) in pairscore.keys():
                continue
            similarity = pairscore[(shortctg,longctg)]/pairscore[(shortctg,shortctg)]

            res = min(similarity,res)
    return res
def attach(scorefile,majorfile,complementfile,ratio):


    
    complementdf = pd.read_csv(complementfile,sep=',',header=0)
    complements = list(complementdf['ctg'])

    with open('tmp/c2members') as f1:
        c2members = json.load(f1)
    with open('tmp/c2size') as f2:
        c2size = json.load(f2)
    with open('tmp/c2sim') as f3:
        c2sim = json.load(f3)
    

    apexamplars = []
    c2membersf = {}
    for line in open(majorfile):
        line = line.strip()
        info = line.split('\t')
        ctgs = info[0].split(',')
        examplars = info[1].split(',')
        
        for ctg,examplar in zip(ctgs,examplars):
            apexamplars.append(examplar)
            if examplar in c2membersf.keys():
                c2membersf[examplar].append(ctg)
            else:
                c2membersf[examplar] = [ctg]

    for examplar,members in c2members.items():
        c2membersf[examplar] = members

    resdata = []

    for examplar,members in c2membersf.items():
        memberstr = ','.join(members)
        resdata.append(memberstr)
    resdf = pd.DataFrame(resdata,columns=['Members'])
  

    major = list(set(apexamplars))
    
    scoredf = pd.read_csv(scorefile,sep='\t',header=0)
    pairscore = {}
    for ctg1,ctg2,score in zip(scoredf['ctg1'],scoredf['ctg2'],scoredf['score']):
        ctg1 = ctg1.replace(' ','~')
        ctg2 = ctg2.replace(' ','~')

        pairscore[(ctg1,ctg2)] = score    
    branches = []
    m2c = {}
    for mcenter in major:
        for ccenter in complements:
            if pairscore[(mcenter,mcenter)] < pairscore[(ccenter,ccenter)]:
                shortctg = mcenter
                longctg = ccenter
            else:
                shortctg = ccenter
                longctg = mcenter
            
            if not (shortctg,longctg) in pairscore.keys():
                continue
            similarity = pairscore[(shortctg,longctg)]/pairscore[(shortctg,shortctg)]
            if similarity >= ratio:
                #print(mcenter,ccenter)
                branches.append([mcenter,ccenter,similarity])

    valid_branches = []
    for mcenter,ccenter,similarity in branches:
        if mcenter in m2c.keys():
            m2c[mcenter].append((ccenter,similarity))
        else:
            m2c[mcenter] = [(ccenter,similarity)]
    for mcenter in m2c.keys():
        m2c[mcenter] = sorted(m2c[mcenter],key=lambda d:d[1])
    for mcenter in m2c.keys():
        if  len(m2c[mcenter]) == 1:
            valid_branches.append((mcenter,m2c[mcenter][0][0],m2c[mcenter][0][1]))
            continue

        N = len(m2c[mcenter])

        for i in range(N):
            
            minscore = minsim(m2c[mcenter][i+1:N],pairscore)
            if minscore !=10 and minscore >= ratio:
                break

        for j in range(i,N):        
            valid_branches.append((mcenter,m2c[mcenter][j][0],m2c[mcenter][j][1]))
          
                
            
    if True:
        taxdf = pd.read_csv('tax.tab',sep='\t',header=0)
        ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))
        ctg2genus = dict(zip(ctgkey,taxdf['genus']))
        num = 0
        wrong= 0
        for mcenter,ccenter,similarity in branches:
            if not mcenter in ctg2genus.keys():
                continue
            if not ccenter in ctg2genus.keys():
                continue

            if ctg2genus[mcenter] == ctg2genus[ccenter]:
                if similarity < c2sim[ccenter] or c2sim[ccenter] <=ratio:
                    continue
                num += c2size[ccenter]
                
                print(c2size[ccenter],similarity,c2sim[ccenter])
            else:
                #print(mcenter,ccenter,ctg2genus[mcenter],ctg2genus[ccenter],similarity,c2sim[ccenter] )
                if similarity < c2sim[ccenter] or c2sim[ccenter] <= ratio:
                    continue
                print(mcenter,ccenter,ctg2genus[mcenter],ctg2genus[ccenter],similarity,c2sim[ccenter] )
                wrong += c2size[ccenter]
        print(num,wrong)
        
     
            
runcluster()
