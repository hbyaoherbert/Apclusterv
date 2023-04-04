from sklearn.cluster import AffinityPropagation
from sklearn.cluster import *
import pandas as pd
import argparse
import numpy as np
import adj_ap
import networkx as nx
import json

parser = argparse.ArgumentParser()
parser.add_argument('scorefile')
parser.add_argument('majorfile')
parser.add_argument('complementfile')
parser.add_argument('ratio')
args = parser.parse_args()


configfile = open('config','r')
var2value = {}
var2value['k'] = 1

for line in configfile:
    line = line.strip()
    info = line.split('=')
    if len(info) != 2:
        continue
    var = info[0].strip()
    value = float(info[1].strip())
    if var in var2value.keys():
        var2value[var] = value

configfile.close()
def attach(scorefile,majorfile,complementfile,ratio):

    
    complementdf = pd.read_csv(complementfile,sep=',',header=0)
    complements = list(complementdf['ctg'])

    with open('tmp/c2members') as f1:
        c2members = json.load(f1)
    with open('tmp/c2size') as f2:
        c2size = json.load(f2)
    with open('tmp/c2sim') as f3:
        c2sim = json.load(f3)
    

    apexamplars = []
    c2membersf = {}
    assign = set()
    for line in open(majorfile):
        line = line.strip()
        info = line.split('\t')
        ctgs = info[0].split(',')
        examplars = info[1].split(',')
        
        for ctg,examplar in zip(ctgs,examplars):
            apexamplars.append(examplar)
            if examplar in c2membersf.keys():
                c2membersf[examplar].append(ctg)
            else:
                c2membersf[examplar] = [ctg]
            assign.add(ctg)

  

    major = list(set(apexamplars))
    
    scoredf = pd.read_csv(scorefile,sep='\t',header=0)
    pairscore = {}
    for ctg1,ctg2,score in zip(scoredf['ctg1'],scoredf['ctg2'],scoredf['score']):
        ctg1 = ctg1.replace(' ','~')
        ctg2 = ctg2.replace(' ','~')

        pairscore[(ctg1,ctg2)] = score    
    branches = []
    m2c = {}
    for mcenter in major:
        for ccenter in complements:
            if pairscore[(mcenter,mcenter)] < pairscore[(ccenter,ccenter)]:
                shortctg = mcenter
                longctg = ccenter
            else:
                shortctg = ccenter
                longctg = mcenter
            
            if not (shortctg,longctg) in pairscore.keys():
                continue
            similarity = pairscore[(shortctg,longctg)]/pairscore[(shortctg,shortctg)]
            if similarity >= ratio:
                #print(mcenter,ccenter)
                branches.append([mcenter,ccenter,similarity])
    
    merged = set()
    repset = set()
    for mcenter,ccenter,similarity in branches:
        repset.add(mcenter)
        if similarity < c2sim[ccenter] or c2sim[ccenter] <=ratio:
                    continue
        memberset = set(c2membersf[mcenter]).union(set(c2members[ccenter]))    

        assign = assign.union(set(c2members[ccenter]))
        c2membersf[mcenter] = list(memberset)
        merged.add(ccenter)

    for examplar,members in c2members.items():
        if examplar in merged:
            continue
        if c2sim[examplar] < ratio:
            continue
        memberset = set(members)
        memberset = memberset.difference(assign)
        memberlist = list(memberset)
        if len(memberset) == 0:
            continue
        repset.add(examplar)
        c2membersf[examplar] = memberlist

    repdata = pd.DataFrame(list(repset),columns=['representative'])
    repdata.to_csv('representative'+str(var2value['k'])+'.csv')
    resdata = []

    for examplar,members in c2membersf.items():

        
        memberstr = ','.join(members)
        resdata.append(memberstr)
    resdf = pd.DataFrame(resdata,columns=['Members'])

    resdf.to_csv('finalres'+str(var2value['k'])+'.csv',index=False)
    '''
    valid_branches = []
    for mcenter,ccenter,similarity in branches:
        if mcenter in m2c.keys():
            m2c[mcenter].append((ccenter,similarity))
        else:
            m2c[mcenter] = [(ccenter,similarity)]
    for mcenter in m2c.keys():
        m2c[mcenter] = sorted(m2c[mcenter],key=lambda d:d[1])
    for mcenter in m2c.keys():
        if  len(m2c[mcenter]) == 1:
            valid_branches.append((mcenter,m2c[mcenter][0][0],m2c[mcenter][0][1]))
            continue

        N = len(m2c[mcenter])

        for i in range(N):
            
            minscore = minsim(m2c[mcenter][i+1:N],pairscore)
            if minscore !=10 and minscore >= ratio:
                break

        for j in range(i,N):        
            valid_branches.append((mcenter,m2c[mcenter][j][0],m2c[mcenter][j][1]))
    '''      
                
            
    if False:
        taxdf = pd.read_csv('tax.tab',sep='\t',header=0)
        ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))
        ctg2genus = dict(zip(ctgkey,taxdf['genus']))
        num = 0
        wrong= 0

        for mcenter,ccenter,similarity in branches:
            if not mcenter in ctg2genus.keys():
                continue
            if not ccenter in ctg2genus.keys():
                continue   
            if ctg2genus[mcenter] == ctg2genus[ccenter]:
                if similarity < c2sim[ccenter] or c2sim[ccenter] <=ratio:
                    continue
                num += c2size[ccenter]
                
                print(c2size[ccenter],similarity,c2sim[ccenter])
            else:
                #print(mcenter,ccenter,ctg2genus[mcenter],ctg2genus[ccenter],similarity,c2sim[ccenter] )
                if similarity < c2sim[ccenter] or c2sim[ccenter] <=ratio:
                    continue
                print(mcenter,ccenter,ctg2genus[mcenter],ctg2genus[ccenter],similarity,c2sim[ccenter] )
                wrong += c2size[ccenter]
        print(num,wrong)


        
attach(args.scorefile,args.majorfile,args.complementfile,float(args.ratio))