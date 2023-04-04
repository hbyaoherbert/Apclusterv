import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("refdb",type=str)
parser.add_argument("vdb",type=str)
parser.add_argument("outtax",type=str)
args = parser.parse_args()

r2genus = {}
for line in open(args.refdb,'r'):
    line = line.strip()
    info = line.split('\t')
    rid = info[0]
    taxstr = info[1]
    genus = taxstr.split('#')[1]
    if genus == 'not assigned':
        continue
    r2genus[rid] = genus

vtax = pd.read_csv(args.vdb,sep=',',header=0)
vtax = vtax.fillna(-1)
for contig,genus in zip(vtax['Organism/Name'],vtax['genus']):
    if genus == -1:
        continue
    r2genus[contig] = genus


r_genus = list(r2genus.items())
taxdf = pd.DataFrame(r_genus,columns = ['contig','genus'])
taxdf.to_csv(args.outtax,sep='\t',index=False)

