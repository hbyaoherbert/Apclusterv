import pandas as pd 
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("pcfile",type=str)
parser.add_argument("clres",type=str)
parser.add_argument("tabfile",type=str)
parser.add_argument("ctgid",type=str)
args = parser.parse_args()

clres = pd.read_csv(args.clres,sep=',',header=0)
cluster = 0

ctg2cl = {}
for members in list(clres['Members']):
    for member in members.split(' '):
        ctg2cl[member] = cluster
    cluster += 1

ctgid = pd.read_csv(args.ctgid,sep='\t',header=0)
id2ctg = dict(zip(ctgid['id'],ctgid['ctg']))
tabdf = pd.read_csv(args.tabfile,sep='\t',header=0)
def separate():
    cl2data = {}
    for cl in set(ctg2cl.values()):
        cl2data[cl] = []
    for idx,row in tabdf.iterrows():
        ctgid1 =  id2ctg[row['ctg1']].replace(' ','~')
        ctgid2 =  id2ctg[row['ctg2']].replace(' ','~')
        cl1 = -1
        cl2 = -1
        if ctgid1 in ctg2cl.keys():
            cl1 = ctg2cl[ctgid1]
        if ctgid2 in ctg2cl.keys():
            cl2 = ctg2cl[ctgid2]
        if cl1 != -1 and cl1 == cl2:
            cl2data[cl1].append(row)
    summary = open('cldata/summary','w')
    summary.write('\t'.join(['cl','mine','maxe','mincount','maxcount']))
    summary.write('\n')
    for cl,data in cl2data.items():
        outfile = open('cldata/'+str(cl)+'.data','w')
        outfile.write('\t'.join(tabdf.columns))
        outfile.write('\n')
        mine = 1
        maxe = 0
        mincount = 1000
        maxcount = 0
        for row in data:
            mine = min(row['mine'],mine)
            maxe = max(row['mine'],maxe)
            mincount = min(row['count'],mincount)
            maxcount = max(row['count'],maxcount)
            ctgid1 =  id2ctg[row['ctg1']].replace(' ','~')
            ctgid2 =  id2ctg[row['ctg2']].replace(' ','~')
            outfile.write('\t'.join([ctgid1,ctgid2,str(row['mine']),str(row['count']) ]))
            outfile.write('\n')

        outfile.close()
        summary.write('\t'.join([str(cl),str(mine),str(maxe),str(mincount),str(maxcount)]))
        summary.write('\n')
    summary.close()

def rmcluster():
    errordf = pd.read_csv('cltest/aggrtest.ingenus',sep='\t',header=0)
    for idx,row in errordf.iterrows():
        ctgid1 =  id2ctg[row['ctg1']].replace(' ','~')
        ctgid2 =  id2ctg[row['ctg2']].replace(' ','~')
        cl1 = -1
        cl2 = -1
        if ctgid1 in ctg2cl.keys():
            cl1 = ctg2cl[ctgid1]
        if ctgid2 in ctg2cl.keys():
            cl2 = ctg2cl[ctgid2]
        if cl1 != -1 or cl2 !=-1:
            continue
        print('\t'.join([ctgid1,ctgid2,str(row['e']),str(row['count']),row['genus1'],row['genus2']]))

rmcluster()

