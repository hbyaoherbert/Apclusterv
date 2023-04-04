import networkx as nx
import numpy as np
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('edgefile',type=str)
parser.add_argument('aggr',type=str)
args = parser.parse_args()
aggr = pd.read_csv(args.aggr,sep='\t',header=0)

edge2e = {}
for e,score,ctg1,ctg2 in zip(aggr['mine'],aggr['score'],aggr['ctgname1'],aggr['ctgname2']):
    ctg1 = ctg1.replace(' ','~')
    ctg2 = ctg2.replace(' ','~')
    edge2e[ctg1+'#'+ctg2] = e

G = nx.Graph()
edgeset = set()
for line in open('c1.ntw','r'):
    line = line.strip()
    info = line.split(' ')
    node1 = info[0]
    node2 = info[1]
    key1 = node1+'#'+node2
    key2 = node2+'#'+node1


    if key1 in edgeset or key2 in edgeset:
        continue
    edgeset.add(key1)

    if key1 in edge2e.keys():
        
        G.add_edge(node1,node2,weight=edge2e[key1])
    elif key2 in edge2e.keys():
    
        G.add_edge(node1,node2,weight=edge2e[key2])

    continue
    w = float(info[2])
    w = np.exp(-w)
    if node1+'#'+node2 in edgeset:
        continue
    if node2+'#'+node1 in edgeset:
        continue
    edgeset.add(node1+'#'+node2)
    G.add_edge(node1,node2,weight=w)
#n = len(G.nodes())
#m = len(G.edges())
#print( m*2 / (n*(n-1)))
#between = nx.edge_betweenness_centrality(G)


between = nx.betweenness_centrality(G,weight="weight")
#cluster = nx.clustering(G)
print('node\tbetween')
for node,num in between.items():
    print(node+'\t'+str(num))

