from sklearn.cluster import AffinityPropagation
from sklearn.cluster import *
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
import numpy.ma as ma
import pandas as pd
import argparse
import adj_ap

parser = argparse.ArgumentParser()
parser.add_argument('pcfile',type=str)
parser.add_argument('clusterlist',type=str)
parser.add_argument('tax',type=str)
args = parser.parse_args()

def readnode(nodefile):
    node2id = {}
    nodes = []
    for line in open(nodefile,'r'):
        line = line.strip()
        nodes.append(line)
        thisid = len(node2id)
        node2id[line] = thisid
    return nodes


def genvec(pcfile,clusterlist,tax):
    pcdf = pd.read_csv(pcfile,sep=',',header=0)
    pcdf = pcdf.dropna()
    ctgkey = pcdf['contig_id'].map(lambda s:s.replace(' ','~'))
    pcdf['ctgkey'] = ctgkey

    for line in open(clusterlist,'r'):
        line = line.strip()
        nodefile = 'cloverlap/'+line+'.nodes'
        nodes = readnode(nodefile)
        edgefile = 'cloverlap/'+line+'.edges'

        clusterdf = pcdf[pcdf['ctgkey'].isin(nodes)]
        
    
        vecdf = pd.crosstab(clusterdf['ctgkey'],clusterdf['cluster'])
        #vecdf.to_csv(nodefile+'.vec',sep=',')
        #print(len(vecdf.columns))
        id2node = {}
        node2id = {}
        nodes = []
        for ctgkey in vecdf.index:
            thisid = len(id2node)
            id2node[thisid]= ctgkey
            node2id[ctgkey] = thisid
            nodes.append(ctgkey)
        vecnp = np.clip(vecdf.values,0,1)
        mtx= -euclidean_distances(vecnp, squared=True)
        edges = pd.read_csv(edgefile,sep='\t',header=None)
        edges.columns = ['n1','n2','w']
        mask = ma.make_mask(mtx)
         
        for n1,n2,w in zip(edges['n1'],edges['n2'],edges['w']):
            id1 = node2id[n1]
            id2 = node2id[n2]
            mask[id1,id2] = False
            mask[id2,id1] = False
        if len(mtx)<50:
            continue
        

        mtx[mask] = -len(vecdf.columns)
         
        print(np.median(mtx),np.min(mtx),np.max(mtx))
        nmax = np.max(mtx)
        nmin = np.min(mtx)
        
        
        #ap2 = AffinityPropagation(affinity="euclidean",random_state = 39,max_iter=5000,preference = nmin)
        #ap2.fit(vecdf.values)
        preference = np.median(mtx) *len(mtx)/10
        pstep = (nmax-nmin)/10
        cluster_centers_indices, labels = adj_ap.adj_affinity_propagation(mtx,damping=0.5,pstep=pstep, max_iter=2000,preference=preference, random_state=39)
        print(nodefile,len(nodes),len(cluster_centers_indices))
        taxdf = pd.read_csv(tax,sep='\t',header=0)
        ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))
        ctg2g = dict(zip(ctgkey,taxdf['genus']))
        mtxdf = pd.DataFrame(mtx,columns = vecdf.index)
        mtxdf.index = vecdf.index
        mtxdf.to_csv(nodefile+'.euc',index=True,sep=',')

        
        outfile = open(nodefile+'.eucap','w')
        for i in range(len(mtx)):
            ctg = id2node[i]
            if ctg in ctg2g.keys():
                genus = ctg2g[ctg]
            else:
                genus = 'not assigned'
            outfile.write(ctg+'\t'+genus+'\t'+str(labels[i])+'\t'+id2node[labels[i]]+'\n')
        outfile.close()
        


def ap(nodefile,edgefile,tax):
    nodes = []
    taxdf = pd.read_csv(tax,sep='\t',header=0)
    ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))
    ctg2g = dict(zip(ctgkey,taxdf['genus']))


    node2id = {}
    for line in open(nodefile,'r'):
        line = line.strip()
        nodes.append(line)
        thisid = len(node2id)
        node2id[line] = thisid

    N = len(nodes)
    mtx = np.ones((N,N)) * -1    
    edges = pd.read_csv(edgefile,sep='\t',header=0)
    edges.columns = ['n1','n2','w']
    nmax = np.max(edges['w'])
    nmin = np.min(edges['w'])
    for n1,n2,w in zip(edges['n1'],edges['n2'],edges['w']):
        nid1 = node2id[n1]
        nid2 = node2id[n2]
        mtx[nid1,nid2] = w/nmax
        mtx[nid2,nid1] = w/nmax

    preference = nmin/nmax
    #mtx = mtx/np.max(edges['w'])
    cluster_centers_indices, labels = affinity_propagation(mtx, max_iter=1000,preference=preference, random_state=39)
    print(nodefile,N,len(cluster_centers_indices))
    outfile = open(nodefile+'.ap','w')
    for i in range(len(mtx)):
        ctg = nodes[i]
        if ctg in ctg2g.keys():
            genus = ctg2g[ctg]
        else:
            genus = 'not assigned'
        outfile.write(ctg+'\t'+genus+'\t'+str(labels[i])+'\n')
    outfile.close()


def runapadj(simg,clcoef,between):
    mtx = np.array(nx.adjacency_matrix(simg).todense())
    

    preference = np.zeros(N)
    #preference = np.median(mtx)
    #print(clcoef)
    pstep = np.zeros(N)
    for node,coef in clcoef.items():
        print(node,coef,between[node])
        preference[node] = coef - 5*between[node]
        pstep[node] = between[node]
    #preference = np.mean(mtx,axis=0)
    pstep = (np.max(mtx) - np.min(mtx) )/100
    print(ctg2id)
    #cluster_centers_indices, labels = affinity_propagation(mtx, max_iter=1000,preference=preference, random_state=39)
    cluster_centers_indices, labels = adj_ap.adj_affinity_propagation(mtx,damping=0.2,pstep=pstep, max_iter=2000,preference=preference, random_state=39)

#for line in open(args.clusterlist,'r'):
    #line = line.strip()
    #ap('cloverlap/'+line+'.nodes','cloverlap/'+line+'.edges',args.tax)

if __name__ == "__main__":
    genvec(args.pcfile,args.clusterlist,args.tax)
