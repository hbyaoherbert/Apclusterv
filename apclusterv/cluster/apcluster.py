import argparse
import pandas as pd
import numpy as np


#parser = argparse.ArgumentParser()
#parser.add_argument("resfile",type=str)
#parser.add_argument("aggr",type=str)
#parser.add_argument("map",type=str)

#args = parser.parse_args()

#usage: filter_edge(args)
#usage: construct_graph('tmp/'+args.aggr+".abs",'tmp/'+args.aggr+".rep")

ctg2idx = {}
idx2ctg = {}
parent = []
replen = []

def readmcl(filename):
	cl2members = {}
	mclres = pd.read_csv(filename,sep=',',header=0)
	for idx,row in mclres.iterrows():
		cl2members[idx] = row['Members'].split(',')
	return cl2members


def filter_edge(args):

	mapdf = pd.read_csv(args.map,sep=',',header=0)
	prot2id = dict(zip(mapdf['protein_id'],mapdf.index))
	ctg2id = dict(zip(mapdf['contig'],mapdf['ctgid']))
	prot2ctg = dict(zip(mapdf.index,mapdf['ctgid']))


	aggr = pd.read_csv(args.aggr,sep='\t',header=0)

	edge2score = {}

	cl2edges = {}
	ctg2cl = {}
	pabs = []

	rep2members = readmcl(args.resfile)
	for rep,members in rep2members.items():
		for member in members:
			ctg2cl[ctg2id[member]] = rep
	for score,ctg1,ctg2 in zip(aggr['alnscore'],aggr['ctg1'],aggr['ctg2']):
		
		if not ctg1 in ctg2cl.keys():
			pabs.append([ctg1,ctg2,score])
			continue
		if not ctg2 in ctg2cl.keys():
			pabs.append([ctg1,ctg2,score])
			continue

		cl1 = ctg2cl[ctg1]
		cl2 = ctg2cl[ctg2]
		if cl1 == cl2:
			pabs.append([ctg1,ctg2,score])


	pabsdf = pd.DataFrame(pabs,columns = ['node1','node2','score'])

	pabsdf.to_csv(args.aggr+".abs",index=False,sep='\t')




def find_rep(ctgid):
    rep = ctgid
    while parent[rep] != rep:
        rep = parent[rep]
    
    while parent[ctgid] != rep:
        nextnode = parent[ctgid]
        parent[ctgid] = rep
        ctgid = nextnode
    
    return rep

def find_rep_merge(ctgid,rep_after_merge):
    

    rep = ctgid
    while parent[rep] != rep:
        rep = parent[rep]
    if rep == rep_after_merge:
        return rep

    if replen[rep] + replen[rep_after_merge] > 200:
        return rep
    while parent[ctgid] != rep:
        nextnode = parent[ctgid]
        parent[ctgid] = rep_after_merge
        ctgid = nextnode

    replen[rep_after_merge] += replen[rep]
    parent[rep] = rep_after_merge

    return rep

def construct(complement_e):
    edges = []
    edge1 = pd.read_csv(complement_e,sep='\t',header=0)
    


    for ctg1,ctg2,score in zip(edge1['node1'],edge1['node2'],edge1['score']):

        if not ctg1 in ctg2idx.keys():
            thisid = len(ctg2idx)
            ctg2idx[ctg1] = thisid
            idx2ctg[thisid] = ctg1
            parent.append(thisid)
            replen.append(1)
        if not ctg2 in ctg2idx.keys():
            thisid = len(ctg2idx)
            ctg2idx[ctg2] = thisid
            idx2ctg[thisid] = ctg2 
            parent.append(thisid)
            replen.append(1)
        edges.append([ctg1,ctg2,score])
    edges = sorted(edges,key=lambda d:d[2],reverse=True)
    
   
    return edges
'''
def readlen(filename):
    ctg2len = {}
    for line in open(filename,'r'):
        line = line.strip()
        info = line.split('\t')
        ctg2len[info[0]] = int(info[1])
    return ctg2len
    

'''
def cluster(edges):
    for ctg1,ctg2,score in edges:
        
        if score < 0.3:
            continue

        ctgid1 = ctg2idx[ctg1]
        ctgid2 = ctg2idx[ctg2]

        longctg = max(ctgid1,ctgid2)
        shortctg = min(ctgid1,ctgid2)

        #if aln/shortlen < 0.75 or pid <90:
            #continue

        
        replong =  find_rep(longctg)
        repshort = find_rep_merge(shortctg,replong)
        #print(shortctg,longctg,repshort,replong,shortlen,longlen,aln,pid)

def printcluster(outname):
    outfile = open(outname,'w')
    for i in range(len(parent)):
        rep = parent[i]
        while rep!=parent[rep]:
            rep = parent[rep]
        
        outfile.write('\t'.join([str(i),str(rep),str(idx2ctg[i])])+'\n')
            
            
def construct_graph(edgedata,clusterout):
#parser = argparse.ArgumentParser()
#parser.add_argument('complement_e')
#parser.add_argument('complement_rep')
#args = parser.parse_args()

	edges = construct(edgedata)
	cluster(edges)
	printcluster(clusterout) 
	
