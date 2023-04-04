import argparse
import pandas as pd
import numpy as np
from scipy import stats
from scipy.special import comb
parser = argparse.ArgumentParser()
parser.add_argument("normscore",type=str)
parser.add_argument("id2ctg",type=str)
parser.add_argument("clusterfile",type=str)


args = parser.parse_args()

def rand_index_score(clusters, classes):
    tp_plus_fp = comb(np.bincount(clusters), 2).sum()
    tp_plus_fn = comb(np.bincount(classes), 2).sum()
    A = np.c_[(clusters, classes)]
    tp = sum(comb(np.bincount(A[A[:, 0] == i, 1]), 2).sum()
             for i in set(clusters))
    fp = tp_plus_fp - tp
    fn = tp_plus_fn - tp
    tn = comb(len(A), 2) - tp - fp - fn

    print(tp,tn,fp,fn)
    return (tp + tn) / (tp + fp + fn + tn)

def summarize(clusterfile,examplar2rep):
    apexamplars = []
    
    classes = []
    ctgset = set()

    for line in open(clusterfile,'r'):
        line = line.strip()
        info = line.split('\t')
        ctgs = info[0].split(',')
        examplars = info[1].split(',')

        for i in range(len(ctgs)):
            
            ctg = ctgs[i]
            if not ctg in ctg2genus.keys():
                continue
            ctgset.add(ctg)
            examplar = examplars[i]
            #examplar2rep[examplar] = examplar
            apexamplars.append(examplar)
            classes.append(genus2id[ctg2genus[ctg]])

    

    clusters = np.ones(len(classes),dtype=int)*100000


    exset = set(apexamplars)
    finalrep = {}
    for examplar in exset:
        rep = examplar
        while examplar2rep[rep] != rep:
            print(rep,examplar2rep[rep] )
            rep = examplar2rep[rep] 

        finalrep[examplar] = ctg2id[rep]


    for idx in range(len(classes)):
        examplar = apexamplars[idx]

        repid = finalrep[examplar]
        clusters[idx] = repid

    clusters = clusters.tolist()
    
    c2size = {}
    for i in clusters:
        if i in c2size.keys():
            c2size[i] += 1
        else:
            c2size[i] = 1
    sensitivity = 0
    for i,num in c2size.items():
        if num>=2:
            sensitivity += num
    print("sensitivity:",sensitivity)
    for ctg,genus in ctg2genus.items():
        if ctg in ctgset:
            continue
        classes.append(genus2id[genus])
        clusters.append(100000)
    
    print(classes[0:20],clusters[0:20])

    print(rand_index_score(clusters,classes))

               
    


id_ctg = pd.read_csv(args.id2ctg,sep='\t',header=0)

ctgs = id_ctg['ctg'].map(lambda s:s.replace(' ','~'))
id2ctg = dict(zip(id_ctg['id'],ctgs))
ctg2id = dict(zip(ctgs,id_ctg['id']))



taxdf = pd.read_csv('tax.tab',sep='\t',header=0)
ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))
ctg2genus = dict(zip(ctgkey,taxdf['genus']))
genus2id = {}
for genus in  list(taxdf['genus']):
    if not genus in genus2id.keys():
        thisid = len(genus2id)
        genus2id[genus] = thisid




aln2score = {}

normdf = pd.read_csv(args.normscore,sep='\t',header=0)

for ctg1,ctg2,score in zip(normdf['ctg1'],normdf['ctg2'],normdf['score']):
    #ctgname1 = id2ctg[ctg1]
    #ctgname2 = id2ctg[ctg2]
    ctg1 = ctg1.replace(' ','~')
    ctg2 = ctg2.replace(' ','~')
    aln2score[(ctg1,ctg2)] = score


'''
for ctg1,ctg2,score in zip(dupincdf['ctg1'],dupincdf['ctg2'],dupincdf['score']):
    ctgname1 = id2ctg[ctg1]
    ctgname2 = id2ctg[ctg2]
    aln2score[(ctgname1,ctgname2)] += score


data = []
for ctg_pair,score in aln2score.items():
    data.append([ctg_pair[0],ctg_pair[1],score])

normdf = pd.DataFrame(data,columns=['ctg1','ctg2','score'])
normdf.to_csv("cltest/aggrnorm.csv",sep='\t',index=False)



'''

gset = set()

g2cluster = {}
idx = 0

examplar2rep = {}
for line in open(args.clusterfile,'r'):
    line = line.strip()
    info = line.split('\t')
    ctgs = info[0].split(',')
    examplars = info[1].split(',')
    for genome in examplars:
        gset.add(genome)
        g2cluster[genome] = idx
        examplar2rep[genome] = genome
    idx += 1 

'''
for line in open(args.apres,'r'):
    line = line.strip()
    for genome in line.split(','):
        gset.add(genome)
        g2cluster[genome] = idx
    idx += 1
'''

g2id = {}
id2g = {}
for genome in gset:
    thisid = len(g2id)
    g2id[genome] = thisid
    id2g[thisid] = genome

N = len(gset)
ratio = 0.15



for ratio in np.arange(0.1,0.3,0.02):
    ingenus = []
    crossgenus = []
    iexamplar2rep = {}
    for examplar,rep in examplar2rep.items():
        iexamplar2rep[examplar] = rep
    for i in range(N):
        for j in range(N):
            x = id2g[i]
            y = id2g[j]

            if not (x,y) in aln2score.keys():
                continue

            if not x in ctg2genus.keys():
                continue
            if not y in ctg2genus.keys():
                continue

            if aln2score[(x,x)] < aln2score[(y,y)]:
                shortctg = x
                longctg = y
            else:
                longctg = x
                shortctg = y
            if not (shortctg,longctg) in aln2score.keys():
                continue

            score = aln2score[(shortctg,longctg)] / aln2score[(shortctg,shortctg)]
        #scoret = aln2score[(y,x)] / aln2score[(x,x)]
            genus1 = ctg2genus[x]
            genus2 = ctg2genus[y]
        
            if genus1 == genus2:
               
                
                if score >= ratio and g2cluster[x] != g2cluster[y]:
                    iexamplar2rep[shortctg] = longctg
                    #print(shortctg,longctg)
                    ingenus.append(score)
            elif score >=ratio:
                
                if g2cluster[x] != g2cluster[y]:
                    #print(x,y,score,scoret)
                    crossgenus.append(score)
                    #print(shortctg,longctg)
                    iexamplar2rep[shortctg] = longctg
    #print(len(ingenus),len(crossgenus))
    #print(np.mean(ingenus),np.mean(crossgenus))
    summarize(args.clusterfile,iexamplar2rep)
    





