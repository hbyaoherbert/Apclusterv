import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("vtax",type=str)
parser.add_argument("rtax",type=str)
parser.add_argument('edgefile',type=str)
args = parser.parse_args()

r2genus = {}


for line in open(args.rtax,'r'):
    line = line.strip()
    info = line.split('\t')
    rid = info[0]
    taxstr = info[1]
    genus = taxstr.split('#')[1]
    if genus == 'not assigned':
        continue
    r2genus[rid] = genus
 
vtax = pd.read_csv(args.vtax,sep=',',header=0)
v2genus = {}
vtax = vtax.fillna(-1)
for contig,genus in zip(vtax['Organism/Name'],vtax['genus']):
    if genus == -1:
        continue
    v2genus[contig] = genus

def gmap(genome):
    if genome in r2genus.keys():
        return r2genus[genome]
    if genome in v2genus.keys():
        return v2genus[genome]
    return "not assigned"

for line in open(args.edgefile,'r'):
    line = line.strip()
    info = line.split(' ')
    node1 = info[0]
    node2 = info[1]
    genus1 = gmap(node1)
    genus2 = gmap(node2)

    if genus1 != genus2:
        print('\t'.join([node1,node2,info[2],genus1,genus2]))



