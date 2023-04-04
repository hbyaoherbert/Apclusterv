import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('resfile',type=str)
parser.add_argument('vtax',type=str)
parser.add_argument('ginfo',type=str)

args = parser.parse_args()

vtax = pd.read_csv(args.vtax,sep=',',header=0)
vtax = vtax.fillna(-1)
v2family = {}
for genome,family in zip(vtax['Organism/Name'],vtax['family']):
    if family == -1:
        continue
    v2family[genome.replace(' ','~')] = family
gpd = pd.read_csv(args.ginfo,sep='\t',header=0)
g2host = {} 
for gid,family,hostinfo in zip(gpd['GPD_id'],gpd['Predicted_phage_taxon'],gpd['Host_range_taxon']):
    hostrange = set()
    for host in hostinfo.split(','):
        taxparts = host.split('/')
        genus = taxparts[-2]
        hostrange.add(genus)

    g2host[gid] = hostrange
    v2family[gid] = family

resdf = pd.read_csv(args.resfile,sep=',',header=0)
num = 0
total = 0


fp = 0
assign = 0
for members in resdf['Members']:
    clrange = set()
    f2num = {}
    gpdnum = 0
    if len(members.split(',') ) < 2:
        continue


    for member in members.split(','):
        if not member in v2family.keys():
            continue
        family = v2family[member]
        if family in f2num.keys():
            f2num[family] += 1
        else:
            f2num[family] = 1
        gpdnum += 1
        if not member in g2host.keys():
            continue
        
        clrange = clrange.union(g2host[member])
        
        
    #print(f2num)
    #if len(clrange) == 0:
        #continue
    if gpdnum <2:
        continue
    assign += gpdnum

    if len(f2num)>1:
        f2numlist = list(f2num.items()) 
        f2numlist = sorted(f2numlist,key=lambda d:d[1],reverse=True)
        tp = f2numlist[0][1]
        fp += gpdnum - tp

    
    total += 1
    if len(clrange) > 1:
        num +=1
    #print(len(clrange))
     
print(assign,fp,(assign-fp)/assign)
print(num,total)




