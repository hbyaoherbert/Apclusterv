import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("aggr",type=str)
parser.add_argument("id2ctg",type=str)
parser.add_argument("info",type=str)
args = parser.parse_args()




def find_rep(ctgid,parent):
    rep = ctgid
    while parent[rep] != rep:
        rep = parent[rep]
    
    while parent[ctgid] != rep:
        nextnode = parent[ctgid]
        parent[ctgid] = rep
        ctgid = nextnode
    
    return rep

def find_rep_merge(ctgid,rep_after_merge,parent):
    

    rep = ctgid
    while parent[rep] != rep:
        rep = parent[rep]
    if rep == rep_after_merge:
        return rep
    while parent[ctgid] != rep:
        nextnode = parent[ctgid]
        parent[ctgid] = rep_after_merge
        ctgid = nextnode
    parent[rep] = rep_after_merge
    return rep

def printcluster(parent,id2ctg,outname):
    rep2ctg = {}
    for i in parent.keys():
        rep = parent[i]
        while rep!=parent[rep]:
            rep = parent[rep]
        if rep in rep2ctg.keys():
            rep2ctg[rep].append(i)
        else:
            rep2ctg[rep] = [i]

    data = []        
    for rep,ctgs in rep2ctg.items():
          
        repctg = id2ctg[rep]
        members = [id2ctg[ctgid] for ctgid in ctgs]
        if len(members) > 10:
            print(members)
        nmember = len(members)
        members = ",".join(members)
        data.append([repctg,members,nmember])

        df = pd.DataFrame(data,columns = ['cluster','members','num'])
        df.to_csv(outname,sep='\t',index=False)
        

def parseaggr(aggr,id2ctg,info):
    df = pd.read_csv(info,sep=',',header=0)
    parent = {}
    iddf = pd.read_csv(id2ctg,sep='\t',header=0)
    ctgname = iddf['ctg'].map(lambda s:s.replace(' ','~'))

    ctg2id = dict(zip(ctgname,iddf['id']))
    id2ctg = dict(zip(iddf['id'],ctgname))

    for ctg,ctgid in ctg2id.items():
        parent[ctgid] = ctgid
    ctgkey = df['contig_id'].map(lambda s:s.replace(' ','~'))
    ctg2nprot = dict(zip(ctgkey,df['proteins']))
    ratio = 0.15
    aggrdf = pd.read_csv(aggr,sep='\t',header=0)
    edges = []
    
    for ctgkey1,ctgkey2,count in zip(aggrdf['ctgname1'],aggrdf['ctgname2'],aggrdf['count']):
        ctgkey1 = ctgkey1.replace(' ','~')
        ctgkey2 = ctgkey2.replace(' ','~')
        reflen = min(ctg2nprot[ctgkey1],ctg2nprot[ctgkey2])
        if reflen<=10:
            continue
        if count*1.0/reflen >=ratio:
            ctgid1 = ctg2id[ctgkey1]
            ctgid2 = ctg2id[ctgkey2]
            edges.append([ctgid1,ctgid2,ratio])
            id1 = min(ctgid1,ctgid2)
            id2 = max(ctgid1,ctgid2)
            rep_after_merge = find_rep(id2,parent)
            find_rep_merge(id1,rep_after_merge,parent)
    '''
    for ctgid1,ctgid2,count,reflen in edges:
        cluster1 = ctg2cl[ctgid1]
        cluster2 = ctg2cl[ctgid2]
        if cluster1 == cluster2：
            cl2edges[cluster1].append(count)
        else:
            cl2outedges[cluster1].append(count)
            cl2outedges[cluster2].append(count)    
    '''
    printcluster(parent,id2ctg,"cltest/simplecl")

parseaggr(args.aggr,args.id2ctg,args.info)
