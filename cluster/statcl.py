import pandas as pd 
import argparse
import numpy as np
import networkx as nx

parser = argparse.ArgumentParser()
parser.add_argument('cluster',type=str)
parser.add_argument('ingenus',type=str)
parser.add_argument('crossgenus',type=str)
parser.add_argument('clstat',type=str)

args = parser.parse_args()



def construct(crossgenus,ingenus,cluster):
    cldf = pd.read_csv(cluster,sep='\t',header=None)
    cldf.columns = ['i','rep','contig']
    ctgset = set(cldf['contig'])
    edges = []
    edge1 = pd.read_csv(crossgenus,sep='\t',header=None)
    edge1.columns = ['ctg1','ctg2','e','count','genus1','genus2']
    
    for ctg1,ctg2,e in zip(edge1['ctg1'],edge1['ctg2'],edge1['e']):
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
        
        edges.append([ctg1,ctg2,e])
    
    edge2 = pd.read_csv(ingenus,sep='\t',header=None)
    edge2.columns = ['ctg1','ctg2','e','count','genus1','genus2']
    for ctg1,ctg2,e in zip(edge2['ctg1'],edge2['ctg2'],edge2['e']):
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
            
        edges.append([ctg1,ctg2,e])
    #print(edges)    
    return edges
#stratification of nodes according to betweeness centrality, which is a global topology measurement that describes how likely a node is to be a bridge to different clustering components
def parsecl(cluster,ingenus,crossgenus,clstat):
    cldf = pd.read_csv(cluster,sep='\t',header=None)
    

    cldf.columns = ['i','rep','contig']
    ctg2id = {}
    id2ctg = {}

    rep2members = {}
    member2rep = {}

    rep2edges = {}
    rep2in = {}
    rep2cross = {}
    for ctgid,repid,ctg in zip(cldf['i'],cldf['rep'],cldf['contig']):
        ctg2id[ctg] = ctgid 
        id2ctg[ctgid] = ctg 
        member2rep[ctg] = repid
        if repid in rep2members.keys():
            rep2members[repid].append(ctg)

        else:
            rep2cross[repid] = []
            rep2edges[repid] = []
            rep2members[repid] = [ctg]

    edges = []
    edge1 = pd.read_csv(crossgenus,sep='\t',header=None)
    edge1.columns = ['ctg1','ctg2','e','count','genus1','genus2']
    
    G = nx.Graph()
    for ctg1,ctg2,e in zip(edge1['ctg1'],edge1['ctg2'],edge1['e']):
        if not ctg1 in ctg2id.keys():
            continue
        if not ctg2 in ctg2id.keys():
            continue
        rep1 = member2rep[ctg1]
        rep2 = member2rep[ctg2]
        
        G.add_edge(ctg1,ctg2,weight=e)
        rep2edges[rep1].append([ctg1,ctg2,e])
       
        rep2cross[rep1].append([ctg1,ctg2,e])

    edge2 = pd.read_csv(ingenus,sep='\t',header=None)
    edge2.columns = ['ctg1','ctg2','e','count','genus1','genus2']
    for ctg1,ctg2,e in zip(edge2['ctg1'],edge2['ctg2'],edge2['e']):
        if not ctg1 in ctg2id.keys():
            continue
        if not ctg2 in ctg2id.keys():
            continue
        e = max(1e-300,e)
        rep1 = member2rep[ctg1]
        rep2 = member2rep[ctg2]
        if rep1 != rep2:
            continue
        
        G.add_edge(ctg1,ctg2,weight=e)
        rep2edges[rep1].append([ctg1,ctg2,e])
        

    
    node2between = nx.betweenness_centrality(G)
   
    rep2e = {}
    outfile = open(clstat,'w')
    for node,between in node2between.items():
        outfile.write(node+'\t'+str(between)+'\n')
    outfile.close()
    for rep,edges in rep2edges.items():
        between = []
        evalue = []
        identity = []
        aln = []
        for ctg1,ctg2,e in edges:
            evalue.append(e)
        if len(evalue) == 0:
            continue
        rep2e[rep] = np.median(evalue)
    for rep,members in rep2members.items():
        between = []
        for member in members:
            if member in node2between.keys():
                between.append(node2between[member])
        
        if len(between) ==0:
            continue
        print(rep,rep2e[rep],np.median(between),np.max(between),len(rep2edges[rep]),len(rep2cross[rep]))


parsecl(args.cluster,args.ingenus,args.crossgenus,args.clstat)
