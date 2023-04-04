import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument("ctgscore",type=str)
parser.add_argument("ctgid",type=str)
parser.add_argument("c1edge",type=str)
parser.add_argument("out",type=str)
args = parser.parse_args()


#ctg2genus = dict(zip(taxdf['contig'],taxdf['genus']))
ctgid = pd.read_csv(args.ctgid,sep='\t',header=0)
id2ctg = dict(zip(ctgid['id'],ctgid['ctg']))


scoredf = pd.read_csv(args.ctgscore,sep='\t',header=0)


edges = pd.read_csv(args.c1edge,sep=' ',header=None)
edges.columns = ['node1','node2','weight']
edgeset = set()
for node1,node2 in zip(edges['node1'],edges['node2']):
    edgeset.add(node1+'#'+node2)

outfile = open(args.out,'w')
for idx,row in scoredf.iterrows():
    ctgid1 = id2ctg[row['ctg1']].replace(' ','~')
    ctgid2 = id2ctg[row['ctg2']].replace(' ','~')
    if not ctgid1+'#'+ctgid2 in edgeset:
         outfile.write('\t'.join([ctgid1,ctgid2,str(row['hsp']),str(row['score'])])+'\n')

outfile.close()
