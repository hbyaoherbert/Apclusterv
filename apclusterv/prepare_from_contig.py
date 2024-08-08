from Bio import SeqIO
import subprocess
import argparse
import pandas as pd
import os

def parsefaa(transfile):
	prot = open(transfile,'r')
	data = []
	for record in SeqIO.parse(prot,'fasta'):
		pos = record.id.rfind('_')
		contig = record.id[0:pos]
		data.append([record.id,contig])
	mapdf = pd.DataFrame(data,columns = ['protein_id','contig_id'])

	mapdf['contig'] = mapdf['contig_id'].map(lambda s:s.replace(' ','~'))

	prot2id = dict(zip(mapdf['protein_id'],mapdf.index))
	ctgset = set(mapdf['contig'])
	ctg2id = {}
	for idx,ctg in enumerate(ctgset):
		ctg2id[ctg] = idx
    

	mapdf['ctgid'] = mapdf['contig'].map(lambda x:ctg2id[x])
	prot2ctg = dict(zip(mapdf.index,mapdf['ctgid']))

	mapdf[['protein_id','contig','ctgid']].to_csv('tmp/prot_contig_id.csv',sep=',',index=False)

def parsecsv(csvlist):
	dflist = []
	for csvfile in csvlist:
		df = pd.read_csv(csvfile,sep=',',header=0)
		dflist.append(df)
		
	mapdf = pd.concat(dflist)
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
	parser.add_argument("contig",help="contig fasta file")
	args = parser.parse_args()
	print("Predict genes (protein sequences) from contig dna sequence")
	#subprocess.run(["prodigal", "-i", args.contig, "-p", "meta" ,"-a", args.contig+".faa", "-o", "genelog"])
	#parsefaa(args.contig+".faa")
	
def fromcsvlist(csvlist):
	parsecsv(csvlist)	
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	#parser.add_argument("contig",help="contig fasta file")
	parser.add_argument("csv1",help="contig fasta file")
	parser.add_argument("csv2",help="db fasta file")
	args = parser.parse_args()
	fromcsvlist([args.csv1,args.csv2])
	
	
	
	
	
	
