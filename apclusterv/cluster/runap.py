from sklearn.cluster import AffinityPropagation
from sklearn.cluster import *
import pandas as pd
import numpy as np
from apclusterv.cluster import adj_ap
import networkx as nx
import json
import argparse

#parser = argparse.ArgumentParser()
#parser.add_argument('c1abs')
#parser.add_argument('cluster')
#usage: runcluster(args)

var2value = {}
var2value['k'] = 10




complements = []
c2size = {}
c2sim = {}
c2members = {}


def sumcluster(cluster):
	cldf = pd.read_csv(cluster,sep=',',header=0)

	rep2ctg = {}
	ctg2rep = {}
	for idx,row in cldf.iterrows():
		ctg2rep[row['contig']] = row['rep']
		if row['rep'] in rep2ctg.keys():
			rep2ctg[row['rep']].append(row['contig'])
		else:
			rep2ctg[row['rep']] = [row['contig']]
	#print(len(rep2ctg))
	return rep2ctg,ctg2rep

def separate(rep2ctg,ctg2rep,c1abs):
	rep2edges = {}
	for rep in rep2ctg.keys():
		rep2edges[rep] = []

	edge1 = pd.read_csv(c1abs,sep=',',header=0)
	
	

	for ctg1,ctg2,score in zip(edge1['ctg1'],edge1['ctg2'],edge1['ratioscore']):
		if not ctg1 in ctg2rep.keys():
			continue
		if not ctg2 in ctg2rep.keys():
			continue

		if ctg2rep[ctg1] != ctg2rep[ctg2]:
			continue
		
		rep2edges[ctg2rep[ctg1]].append([ctg1,ctg2,score])

	
	
	return rep2edges

def printresult():

	with open('tmp/c2members','w') as f1:
		json.dump(c2members,f1)
	with open('tmp/c2size','w') as f2:
		json.dump(c2size,f2)
	with open('tmp/c2sim','w') as f3:
		json.dump(c2sim,f3)

	
def runcluster(args):
	var2value['k'] = args.apk
	rep2ctg,ctg2rep = sumcluster(args.repfile)
	rep2edges = separate(rep2ctg,ctg2rep,args.edgefile)

	for rep,edges in rep2edges.items():

		
		ctg2id,id2ctg = construct(edges)
		if len(ctg2id) < 5:
			continue
		#outfile = open('sres/{}.res'.format(rep),'w')
		aptest(edges,ctg2id,id2ctg)
		#outfile.close()
	printresult()
	

def construct(edges):
	ctg2id = {}
	id2ctg = {}

	for ctg1,ctg2,score in edges:
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
	for ctg1,ctg2,score in edges:
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
		#print('k',var2value['k'])
		#preference[node] = coef - var2value['k']*between[node]
		preference[node] = coef
		pstep[node] = 0.01
	preference = np.clip(preference,0,1)
	#preference = np.mean(mtx,axis=0)
	#pstep = (np.max(mtx) - np.min(mtx) )/100
	#print(ctg2id)
	#cluster_centers_indices, labels = affinity_propagation(mtx, max_iter=1000,preference=preference, random_state=39)
	#print(preference)
	cluster_centers_indices, labels = adj_ap.adj_affinity_propagation(mtx,damping=0.2,pstep=pstep, max_iter=2000,preference=preference, random_state=39)
	centerstr  = ''

	#print(cluster_centers_indices)
	for center in cluster_centers_indices:
		centerstr+=str(id2ctg[center])+','
		complements.append(id2ctg[center])
		c2size[id2ctg[center]] = 0
		c2sim[id2ctg[center]] = -1
		c2members[id2ctg[center]] = []

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
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument('edgefile')
	parser.add_argument('repfile')
	args = parser.parse_args()
	runcluster(args)
