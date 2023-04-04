import argparse
import pandas as pd
import numpy as np
from scipy import stats
import networkx as nx
import adj_ap
import cloverlap
parser = argparse.ArgumentParser()
parser.add_argument("resfile",type=str)
parser.add_argument("aggr",type=str)
parser.add_argument("out",type=str)


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

aggr = pd.read_csv(args.aggr,sep='\t',header=0)




def runapadj(mtx,node2id,clcoef,between,ngenera):
    #mtx = np.array(nx.adjacency_matrix(simg).todense())
    
    N = len(mtx)
    preference = np.zeros(N)
    #preference = np.median(mtx)
    #print(clcoef)
    pstep = np.zeros(N)
    id2node = {}
    for node,coef in clcoef.items():
        nodeid = node2id[node]
        id2node[nodeid] = node
        #print(node,coef,between[node])
        preference[nodeid] = coef - var2value['k']*between[node]
        pstep[nodeid] = 0.01
    preference = np.clip(preference,0,1)
    preference = preference * np.median(mtx)
    #preference = np.mean(mtx,axis=0)
    #pstep = (np.max(mtx) - np.min(mtx) )/100
    #print(ctg2id)
    #cluster_centers_indices, labels = affinity_propagation(mtx, max_iter=1000,preference=preference, random_state=39)
    cluster_centers_indices, labels = adj_ap.adj_affinity_propagation(mtx,damping=0.2,pstep=pstep, max_iter=2000,preference=preference, random_state=39)
    print(len(cluster_centers_indices),ngenera)
    
    centers = [ id2node[centerid] for centerid in cluster_centers_indices]

    clusters = [id2node[cluster_centers_indices[i]] for i in labels]

    return centers,clusters




edge2score = {}

cl2edges = {}
ctg2cl = {}

for hsp,score,ctg1,ctg2 in zip(aggr['hsp'],aggr['score'],aggr['ctgname1'],aggr['ctgname2']):
    ctg1 = ctg1.replace(' ','~')
    ctg2 = ctg2.replace(' ','~')
    edge2score[ctg1+'#'+ctg2] = (score,hsp)

rep2members,cl2rep = cloverlap.collectoverlap(args.resfile)
for rep,members in rep2members.items():
    for member in members:
        ctg2cl[member] = rep
    cl2edges[rep] = []

'''
for idx,row in cldf.iterrows():
    cl2edges[idx] = []
    for member in row['Members'].split(' '):
        ctg2cl[member] = idx
'''


edgeset = set()
for line in open('c1.ntw','r'):
    line = line.strip()
    info = line.split(' ')
    node1 = info[0]
    node2 = info[1]

    if not node1 in ctg2cl.keys():
        continue
    if not node2 in ctg2cl.keys():
        continue

    cl1 = ctg2cl[node1]
    cl2 = ctg2cl[node2]
    if cl1 != cl2:
        continue


    key1 = node1+'#'+node2
    key2 = node2+'#'+node1


    if key1 in edgeset or key2 in edgeset:
        continue

    edgeset.add(key1)
    
    
    if key1 in edge2score.keys():
        edgekey = key1
    elif key2 in edge2score.keys():
        edgekey = key2
    else:
        continue
    score,hsp = edge2score[edgekey]

    cl2edges[cl2rep[cl1]].append([node1,node2,score,hsp])
    

generas = []

pure = []
chimeric = []
betweendist = []
pure4viz = []
chimeric4viz = []

apres = open(args.out,'w')

for cluster in rep2members.keys():
    genera = set()
    #print(np.min(cl2edges[idx]),np.max(cl2edges[idx]),np.median(cl2edges[idx]))
    
    simg = nx.Graph()
    #simg2 = nx.Graph() #unsymmetric similarity score
    
    scores = []
    if len(cl2edges[cluster]) == 0:
        continue
    for node1,node2,score,hsp in cl2edges[cluster]:
        
        simg.add_edge(node1,node2,weight=score/hsp)
        #simg2.add_edge(node1,node2,weight=score)
        scores.append(score/hsp)
    
    '''
    mtx2 = np.array(nx.adjacency_matrix(simg2).todense())
    N = len(mtx2)
    selfscore = [edge2score[node+'#'+node][0] for node in simg2.nodes]
    selfscore = np.array(selfscore).reshape(N,1)
   
    #print(selfscore)
    mtx_unsymmetric = mtx2/selfscore
    idx = 0
    node2id = {}
    for node in simg2.nodes:
        node2id[node]  = idx
        idx += 1
    
    '''
    idx = 0
    node2id = {}
    for node in simg.nodes:
        node2id[node]  = idx
        idx += 1
    mtx = np.array(nx.adjacency_matrix(simg).todense())
    
    clustering = nx.clustering(simg,weight="weight")
    clcoef = list(clustering.values())

    bc = nx.betweenness_centrality(simg)


    bcnpy = np.array(list(bc.values()))

    betweens = bcnpy[~np.isnan(bcnpy)]
    betweenlist = betweens.tolist()
    betweendist = betweendist + betweenlist
    
    #betweendist = np.concatenate( (betweens,betweendist),axis=None)
    
    if len(clcoef) == 0 or len(betweens)==0:
        continue


    betweenh = betweens[betweens > 0.001]
    
    mom = np.median(scores)/np.mean(scores) 
    mstat = abs(1-mom)
    
    mstat = np.mean(scores)
    
    #mstat = len(betweenh)
    if len(rep2members[cluster])>=10:
      
        #if len(betweenh) == 0:
            #continue
        centers,subclusters = runapadj(mtx,node2id,clustering,bc,1)
        
        nodestr = ','.join(simg.nodes)
        centerstr = ','.join(centers)
        clusterstr = ','.join(subclusters)
        apres.write("{}\t{}\n".format(nodestr,clusterstr))
        
        #print(row['Members'],len(genera),mstat,np.max(betweens),len(betweenh))
    else:
        nodestr = ','.join(rep2members[cluster])
        clusterstr = ','.join(rep2members[cluster])
        apres.write("{}\t{}\n".format(nodestr,clusterstr))


apres.close()
'''
betweendist = np.array(betweendist)
np.save("cltest/"+args.out+".betweendist.npy",betweendist)
#cldf['genus'] = generas
puredf = pd.DataFrame(pure4viz,columns=['nGenera','momstat','between'])
puredf.to_csv("cltest/"+args.out+".pure4viz.csv",sep='\t',index=False)
chimeric4viz = pd.DataFrame(chimeric4viz,columns=['nGenera','momstat','between'])
chimeric4viz.to_csv("cltest/"+args.out+".chimeric4viz.csv",sep='\t',index=False)
'''

#cldf.to_csv('viral_cluster_res',sep=',')
