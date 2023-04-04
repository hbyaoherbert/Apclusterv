import pandas as pd
import argparse


def find_rep(parent,clusterid):
    rep = clusterid
    while parent[rep] != rep:
        rep = parent[rep]
    
    while parent[clusterid] != rep:
        nextnode = parent[clusterid]
        parent[clusterid] = rep
        clusterid = nextnode
    
    return rep

def find_rep_merge(parent,clusterid,rep_after_merge):
    

    rep = clusterid
    while parent[rep] != rep:
        rep = parent[rep]
    if rep == rep_after_merge:
        return rep
    while parent[clusterid] != rep:
        nextnode = parent[clusterid]
        parent[clusterid] = rep_after_merge
        clusterid = nextnode
    parent[rep] = rep_after_merge
    return rep

def printcluster(parent,cl2members,outdir):
    rep2members = {}
    member2rep = {}
    rep2edges = {}
    cl2rep = {}
    for i in parent.keys():
        rep = parent[i]
        while rep!=parent[rep]:
            rep = parent[rep]
        cl2rep[i] = rep
        if rep in rep2members.keys():
            rep2members[rep] = rep2members[rep].union(set(cl2members[i]))
        else:
            rep2members[rep] = set(cl2members[i])
    return rep2members,cl2rep
        #print('\t'.join([str(i),str(rep),str(len(cl2members[i])),str(len(rep2members[rep]))   ]))
    '''
    for rep,members in rep2members.items():
        print(rep+'.'+'#'.join(rep2cl[rep]))
        #print(rep)
        rep2edges[rep] = []
        outfile = open(outdir+'/'+rep+'.nodes','w')
        
        for member in members:
            member2rep[member] = rep
            outfile.write(member+'\n')
        outfile.close()
    '''
    return rep2members,cl2rep
    '''
    edges = pd.read_csv(edgelist,sep=' ',header=None)
    edges.columns = ['n1','n2','w']
   
    for n1,n2,w in zip(edges['n1'],edges['n2'],edges['w']):
        if not n1 in member2rep.keys():
            continue
        if not n2 in member2rep.keys():
            continue
        rep1 = member2rep[n1]
        rep2 = member2rep[n2]
        if rep1 == rep2:
            rep2edges[rep1].append([n1,n2,w])
    for rep,edges in rep2edges.items():
        outfile = open(outdir+'/'+rep+'.edges','w')
        for n1,n2,w in edges:
            outfile.write(n1+'\t'+n2+'\t'+str(w)+'\n')
        outfile.close()
   '''

def collectoverlap(clusterfile):
    
    
    parent = {}
    cl2members = {}
    cldf = pd.read_csv(clusterfile,sep=',',header=0,index_col=0)
    g2cl = {}
    

    for idx,row in cldf.iterrows():
        
        cl2members[str(idx)] = row['Members'].split(' ')
        for genome in row['Members'].split(' '):
            if genome in g2cl.keys():
                g2cl[genome].append(idx)
            else:
                g2cl[genome] = [idx]
    unionset = set()

    g2union = {}
    union2cl = {}
    data = []
    for g in g2cl.keys():
        clusterset = set()
        if len(g2cl[g]) == 1:
            continue
        unionidx = []
        for clusteridx in g2cl[g]:
            parent[str(clusteridx)] = str(clusteridx)
            unionidx.append(str(clusteridx))

        unionkey = '#'.join(unionidx)
        if unionkey in unionset:
            continue
        unionset.add(unionkey)
        
    for unionkey in unionset:
        clusters = unionkey.split('#')
        cl1 = clusters[0]
        cl2 = clusters[1]
        replong =  find_rep(parent,cl2)
        repshort = find_rep_merge(parent,cl1,replong)

        for i in range(2,len(clusters)):
            cl1 = clusters[i-1]
            cl2 = clusters[i]
            replong =  find_rep(parent,cl2)
            repshort = find_rep_merge(parent,cl1,replong)
         

    

    rep2members,cl2rep = printcluster(parent,cl2members,"tem/overlap")
    
    for cluster,members in cl2members.items():
        if cluster in cl2rep.keys():
            continue
        rep2members[cluster] = members
        cl2rep[cluster] = cluster
    return rep2members,cl2rep

    

