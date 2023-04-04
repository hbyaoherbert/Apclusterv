import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('complement_e')
parser.add_argument('complement_rep')

args = parser.parse_args()
ctg2id = {}
id2ctg = {}
parent = []
replen = []

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
    edge1 = pd.read_csv(complement_e,sep='\t',header=None)
    edge1.columns = ['ctg1','ctg2','hsp','score']


    for ctg1,ctg2,score in zip(edge1['ctg1'],edge1['ctg2'],edge1['score']):

        if not ctg1 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg1] = thisid
            id2ctg[thisid] = ctg1
            parent.append(thisid)
            replen.append(1)
        if not ctg2 in ctg2id.keys():
            thisid = len(ctg2id)
            ctg2id[ctg2] = thisid
            id2ctg[thisid] = ctg2 
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
        
        if score < 50:
            continue

        ctgid1 = ctg2id[ctg1]
        ctgid2 = ctg2id[ctg2]

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
        
        outfile.write('\t'.join([str(i),str(rep),id2ctg[i]])+'\n')
            

edges = construct(args.complement_e)
cluster(edges)
printcluster(args.complement_rep) 
