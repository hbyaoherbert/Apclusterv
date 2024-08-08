import pandas as pd 
import argparse
from apclusterv.utils.protein_mcl import *
from apclusterv.exec.maincluster import *
import pkg_resources
#mcl: mcl result
#mapfile: prot2contigmap
#protmap:  protein cluster map
#protein_id,contig_id,cluster

#def make_protein_clusters_mcl(blast_fp, out_p, inflation=2):



def parsemcl(mcl,mapfile):
	pcidx = 0
	mapdf = pd.read_csv(mapfile,sep=',',header=0)
	prot2ctg = dict(zip(mapdf['protein_id'],mapdf['contig']))
	data = []
	for line in open(mcl,'r'):
		line = line.strip()
		info = line.split('\t')
		
		if len(info)<2:
			continue
		pc = "PC_"+str(pcidx)
		for protein in info:
			#print(protein,prot2ctg[protein],pc)
			data.append([protein,prot2ctg[protein],pc])
		pcidx += 1
		
	df = pd.DataFrame(data,columns=['protein_id','contig_id','cluster'])
	df.to_csv('tmp/protein_cluster_v2.csv',index=False)
	
def parsecsv(csvname):
	mapdf = pd.read_csv(csvname,sep=',',header=0)
	mapdf['contig'] = mapdf['contig_id'].map(lambda s:s.replace(' ','~'))
	prot2id = dict(zip(mapdf['protein_id'],mapdf.index))
	ctgset = set(mapdf['contig'])
	ctg2id = {}
	for idx,ctg in enumerate(ctgset):
		ctg2id[ctg] = idx
    

	mapdf['ctgid'] = mapdf['contig'].map(lambda x:ctg2id[x])
	prot2ctg = dict(zip(mapdf.index,mapdf['ctgid']))

	mapdf[['protein_id','contig','ctgid']].to_csv('tmp/prot_contig_id.csv',sep=',',index=False)		
	
def main():
	if not os.path.exists("tmp"):
		os.makedirs("tmp")
	parser = argparse.ArgumentParser()
	parser.add_argument("-contig",type=str,help='contig dna file for clustering')
	parser.add_argument("-protein",type=str,help='protein sequences translated from the contigs')
	parser.add_argument("-csv",type=str,help="mapping file for proteins and contigs")
	parser.add_argument("-t",type=int,default=1,help="number of threads for alignment, 1 by default")
	parser.add_argument("-evalue",type=float,default=0.0001,help="evalue threshold for protein-protein alignment, 0.0001 by default")
	parser.add_argument("-alnnum",type=int,default=25,help="max target sequence number for protein-protein alignment, 25 by default")
	parser.add_argument("-i",type=float,default=3,help="inflation value for mcl,3 by default" )
	parser.add_argument("-r",type=float,default=4,help="relative simialrity cut-off parameter in the cluster integration step. 4 by default")

	
	args = parser.parse_args()
	if args.protein is not None:
		protfile = args.protein
		if args.csv is None:
			print("you need specify -csv argument if -protein is used")
			exit()
		parsecsv(args.csv)
	else:
		if args.contig is None:
			print("you need specify -contig or -protein")
			exit()
		protfile = args.contig+".faa"
	
	
	print("Running diamond alignment")
	db_fp = make_diamond_db(protfile, 'tmp', args.t)
	run_diamond(protfile, db_fp, args.t, args.evalue, args.alnnum, 'tmp/protein.diamond.tab')
	
	print("Creating protein family")
	
	mclres = make_protein_clusters_mcl('tmp/protein.diamond.tab', 'tmp')

	parsemcl(mclres,"tmp/prot_contig_id.csv")

	clusterargs = Arg()
	setattr(clusterargs,'map',"tmp/prot_contig_id.csv")
	setattr(clusterargs,'tab','tmp/protein.diamond.tab')
	setattr(clusterargs,'prot',"tmp/protein_cluster_v2.csv")
	setattr(clusterargs,'inflation',args.i)
	setattr(clusterargs,'apk',args.r)
	maincluster(clusterargs)


if __name__ == "__main__":
	print("please run with apclusterv")
	main()
	
