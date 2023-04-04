import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("clres",type=str)
parser.add_argument("pcfile",type=str)
args = parser.parse_args()

cldf = pd.read_csv(args.clres,sep=',',header=0)

member2cl = {}
for idx,row in cldf.iterrows():
    for member in row['Members'].split(' '):
        member2cl[member] = idx

pc = pd.read_csv(args.pcfile,sep=',',header=0)
pc = pc.dropna()

def mapctg(ctg):
    ctg = ctg.replace(' ','~')
    if ctg in member2cl.keys():
        return member2cl[ctg]
    else:
        return -1
pc['clres'] = pc['contig_id'].map(mapctg)
pc[['protein_id','contig_id','cluster','clres']].to_csv(args.pcfile+'.cl',sep=',',index=False)
    


