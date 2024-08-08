import pandas as pd
from scipy.stats import hypergeom
import argparse
import numpy as np
from multiprocessing import Pool
from timeit import default_timer as timer
import subprocess
import os
import pkg_resources
'''
parser = argparse.ArgumentParser()
parser.add_argument('map',type=str)
parser.add_argument('aggr',type=str)
parser.add_argument('prot',type=str)
args = parser.parse_args()
'''

scorep = {}
def crossctg(pc,ctgset):
	ctglist = list(ctgset)
	N = len(ctglist)
	num = 0
	for i in range(N-1):
		for j in range(i+1,N):
			ctg1 = ctglist[i]
			ctg2 = ctglist[j]
			ctga = min(ctg1,ctg2)
			ctgb = max(ctg1,ctg2)
			ctgkey = (ctga,ctgb)
			if not ctgkey in scorep.keys():
				continue
			num +=1
			scorep[ctgkey] += 1
	#print(pc,num)

def aggrpc(df,ctg2id):
	df = df.dropna()
	df = df.sort_values(by='cluster')
	current_pc = ''
	ctgset = set()
	for prot,ctg,pc in zip(df['protein_id'],df['contig_id'],df['cluster']):
		ctg = ctg.replace(' ','~')
		ctg = ctg2id[ctg]
		if pc == current_pc:
			ctgset.add(ctg)
		else:
			if current_pc !='':
				crossctg(current_pc,ctgset)
			ctgset = set()
			current_pc = pc
	crossctg(current_pc,ctgset)


'''
parser = argparse.ArgumentParser()
parser.add_argument('map',type=str)
parser.add_argument('aggr',type=str)
parser.add_argument('prot',type=str)
args = parser.parse_args()
'''
def calglobal(args):

	mapdf = pd.read_csv(args.map,sep=',',header=0)
	prot2id = dict(zip(mapdf['protein_id'],mapdf.index))
	ctg2id = dict(zip(mapdf['contig'],mapdf['ctgid']))
	id2ctg = dict(zip(mapdf['ctgid'],mapdf['contig']))

	ctg2npc = {}
	df = pd.read_csv(args.prot,sep=',',header=0)
	df = df.sort_values(by='contig_id')
	current_contig = ''
	pcset = set()

	allpc = set()
	for contig,pc in zip(df['contig_id'],df['cluster']):
		allpc.add(pc)
		contig = contig.replace(' ','~')
		contig = ctg2id[contig]
		if contig == current_contig:
			pcset.add(pc)
		else:
			if current_contig != '':
				ctg2npc[current_contig] = len(pcset)
			current_contig = contig
			pcset = set()
			pcset.add(pc)



	ctg2npc[current_contig] = len(pcset)

	aggrdf = pd.read_csv(args.aggr,sep='\t',header=0)
	totalnpc = len(allpc)


	n_ctg = len(ctg2npc)

	T = n_ctg*(n_ctg-1)/2

	for ctg1,ctg2,count in zip(aggrdf['ctg1'],aggrdf['ctg2'],aggrdf['alnscore']):


		ctga = min(ctg1,ctg2)
		ctgb = max(ctg1,ctg2)
		ctgkey = (ctga,ctgb)
		scorep[ctgkey] = 0
	'''
	
	n = min(ctg2npc[ctg1],ctg2npc[ctg2])
	N = max(ctg2npc[ctg1],ctg2npc[ctg2])
	prob = hypergeom.sf(count-1,totalnpc,n,N)
	score = np.nan_to_num(-np.log10(prob) - np.log10(T))
	'''
	#print(ctg1,ctg2,scorep)

#print(len(scorep))



	df = pd.read_csv(args.prot,sep=',',header=0)
	print("start aggr")
	aggrpc(df,ctg2id)
	scores = []
	ctgkeys = []
	ns = []
	Ns = []


	for ctgkey,score in scorep.items():
		if score >= 2:
		
			
			scores.append(score-1)
	
			ctg1,ctg2 = ctgkey
			key = "{}#{}".format(ctg1,ctg2)
			ctgkeys.append(key)
			n = min(ctg2npc[ctg1],ctg2npc[ctg2])
			N = max(ctg2npc[ctg1],ctg2npc[ctg2])
			ns.append(n)
			Ns.append(N)

	outdf = pd.DataFrame()
	outdf["ctgkey"] = ctgkeys
	outdf["k"] = scores
	outdf["n"] = ns
	outdf["N"] = Ns
	outdf.to_csv("tmp/rtest.csv",sep=',',index=False)

	dist = pkg_resources.get_distribution('apclusterv')
 
# 获取包的资源文件路径
	
	rfile = dist.location+'/apclusterv/scripts/phypertest.R'
	subprocess.run(['Rscript', rfile,'tmp/rtest.csv'])
	#Rscript /path/to/your/script.R tmp/rtest.csv
	resdf = pd.read_csv('tmp/phyper.csv',header=0,index_col=0)
	resdf['similarity'] = resdf['res'].map(lambda prob:np.nan_to_num(-np.log10(prob) - np.log10(T)))
	outfile = open(args.out,'w')
	cutoff = 10
	if args.mode == "complete":
		cutoff = 5
	for ctgkey,sim in zip(resdf['ctgkey'],resdf['similarity']):
		if sim < cutoff:
			continue
		info = ctgkey.split('#')
		ctg1 = int(info[0])
		ctg2 = int(info[1])
		ctgname1 = id2ctg[ctg1]
		ctgname2 = id2ctg[ctg2]
		outfile.write(str(ctgname1)+'\t'+str(ctgname2)+'\t'+str(sim)+'\n')
		

	outfile.close()
	
