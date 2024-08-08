import argparse
import pandas as pd
import json

#parser = argparse.ArgumentParser()
#parser.add_argument("alnscore",type=str)
#parser.add_argument("pcscore",type=str)
#parser.add_argument("map",type=str)

#args = parser.parse_args()
def bfs(graph, nodes):
	
	seen = set()
	result = []   #List of Lists to hold the final result 
	for node in nodes:
		if node not in seen:
			components = set()
			leaves = [node]
			while leaves:
				leaf = leaves.pop()
				seen.add(leaf)
				components.add(leaf)
				for connected_node in graph[leaf]: 
					if connected_node not in seen: leaves.append(connected_node)
			print(len(components))
			result.append(components) 
	return result
def connect(args):
	alndf = pd.read_csv(args.alnscore,sep='\t',header=0)
	edge2aln = {}
	for ctg1,ctg2,score in zip(alndf['ctg1'],alndf['ctg2'],alndf['ratioscore']):
		ctga = min(ctg1,ctg2)
		ctgb = max(ctg1,ctg2)
		edge2aln[(ctga,ctgb)] = score 
	

	mapdf = pd.read_csv(args.map,sep=',',header=0)
	prot2id = dict(zip(mapdf['protein_id'],mapdf.index))
	ctg2id = dict(zip(mapdf['contig'],mapdf['ctgid']))
	prot2ctg = dict(zip(mapdf.index,mapdf['ctgid']))



	df = pd.read_csv(args.pcscore,header=None,sep='\t')
	ctgset = set()
	graph = {}
	outscore = []
	for idx,row in df.iterrows():
		ctg1 = row[0]
		ctg2 = row[1]
	
		if ctg1==ctg2:
			continue
		ctg1 = ctg2id[ctg1]
		ctg2 = ctg2id[ctg2]

		ctga = min(ctg1,ctg2)
		ctgb = max(ctg1,ctg2)
		key = (ctga,ctgb)
		if key in edge2aln.keys():
			alnscore = edge2aln[key]
			outscore.append([ctga,ctgb,alnscore])


		ctgset.add(ctg1)
		ctgset.add(ctg2)
		if not ctg1 in graph.keys():
			graph[ctg1] = []
		if not ctg2 in graph.keys():
			graph[ctg2] = []
		graph[ctg1].append(ctg2)
		graph[ctg2].append(ctg1)

	scoredf = pd.DataFrame(outscore,columns=['ctg1','ctg2','ratioscore'])
	scoredf.to_csv(args.pcscore+'.aln',index=False,sep=',')
	results = bfs(graph,ctgset)
	print(len(graph))
	print(len(results))
	outdata = []
	idx = 0
	for components in results:
		for ctg in components:
			outdata.append([idx,ctg])
		idx += 1
	outdf = pd.DataFrame(outdata,columns=['rep','contig'])
	outdf.to_csv(args.pcscore+'.rep',sep=',')
