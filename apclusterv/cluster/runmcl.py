import argparse
import networkx as nx
import pandas as pd
import subprocess

#import markov_clustering as mc
#parser = argparse.ArgumentParser()
#parser.add_argument("edges",type=str)
#parser.add_argument("out",type=str)
#parser.add_argument("inflation",type=str)
#args = parser.parse_args()


def run_mcl(args):
	edgedf = pd.read_csv(args.edges,sep='\t',header=None)
	edgedf.columns = ['node1','node2','weight']
	edgeset = set()
	distg = nx.Graph()

	outfile = open('tmp/input2mcl','w')

	for n1,n2,w in zip(edgedf['node1'],edgedf['node2'],edgedf['weight']):
		key1 = n1+'#'+n2
		key2 = n2+'#'+n1
		if key1 in edgeset or key2 in edgeset:
			continue
		edgeset.add(key1)
		outfile.write(n1+'\t'+n2+'\t'+str(w)+'\n')

		distg.add_edge(n1,n2,weight=w)
		
	nodes = list(distg.nodes())
#mtx = nx.to_scipy_sparse_array(distg)
	outfile.close()

	inflation = args.inflation
	subprocess.run(['mcl', args.edges,'--abc','-I', str(inflation),'-o','mclout'])
	'''
	result = mc.run_mcl(mtx,inflation=inflation)
	clusters = mc.get_clusters(result)
	'''
	mclout = open('mclout','r')
	
	data = []
	for line in mclout:
		cl = line.strip().split('\t')
		if len(cl) < 2:
			continue
		members = []
		for member in cl:
			
			members.append(member)
		memberstr = ','.join(members)
		data.append([len(members),memberstr])
	mclout.close()
	df = pd.DataFrame(data,columns=['Size','Members'])
	df.to_csv('tmp/mcl'+str(inflation))
