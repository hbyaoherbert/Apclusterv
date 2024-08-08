import numpy as np
import pandas as pd
from scipy.special import comb
import argparse
from sklearn import metrics
def rand_index_score(clusters, classes):
    tp_plus_fp = comb(np.bincount(clusters), 2).sum()
    tp_plus_fn = comb(np.bincount(classes), 2).sum()
    A = np.c_[(clusters, classes)]
    tp = sum(comb(np.bincount(A[A[:, 0] == i, 1]), 2).sum()
             for i in set(clusters))
    fp = tp_plus_fp - tp
    fn = tp_plus_fn - tp
    tn = comb(len(A), 2) - tp - fp - fn
    N = tp+tn+fp+fn

    print(tp/N,tn/N,fp/N,fn/N)
    return (tp + tn) / (tp + fp + fn + tn)


  
def c1res(resfile,tax):
    cldf = pd.read_csv(resfile,header=0,sep=',')
    ctgset = set()
    for members in cldf['Members']:
        for ctg in members.split(','):
            if ctg.find("NODE_")>=0:
                info = ctg.split('_')
                ctg = info[1]
            ctgset.add(ctg)

    taxdf = pd.read_csv(args.tax,sep='\t',header=0)
    ctgkey = taxdf['contig'].map(lambda s:s.replace(' ','~'))

    ctg2id = {}
    genus2id = {}
    ctg2genus = {}
    classes = []

    
    for ctgkey,genus in  zip(ctgkey,taxdf['genus']):
        if ctgkey.find("NODE_")>=0:
            info = ctgkey.split('_')
            ctgkey = info[1]
        if not ctgkey in ctgset:
            continue
        
        if not genus in genus2id.keys():
            thisid = len(genus2id)
            genus2id[genus] = thisid

        thisid = len(ctg2id)
        ctg2id[ctgkey] = thisid
       
        ctg2genus[ctgkey] = genus2id[genus]
        classes.append(genus2id[genus])


    df = pd.read_csv(args.resfile,header=0,sep=',')
  

    #clusters = np.ones(len(classes),dtype=int) * 100000
    N = len(classes)
    Nres = len(df)
 
    clusters = [i+N for i in range(N)]
    
    outdata = []
    
    nctg = 0
    for idx,row in df.iterrows():
        genus2num = {}
        
        for ctg in row['Members'].split(','):
            if ctg.find("NODE_")>=0:
                info = ctg.split('_')
                ctg = info[1]
            if not ctg in ctg2genus.keys():
                continue
            
            nctg +=1

            genus = ctg2genus[ctg]
            genus2num.setdefault(genus,0)
            genus2num[genus]+=1
            
            clusters[ctg2id[ctg]] = idx
        data = []
        genus2numlist = genus2num.items()
        genus2numlist = sorted(genus2numlist,key=lambda d:d[1],reverse=True)
        for genus,num in genus2numlist:
            data.append(str(genus))
            data.append(str(num))
        outstr = ",".join(data)
        outdata.append([idx,outstr])
    outdf = pd.DataFrame(outdata,columns=['cl','tax'])
    outdf.to_csv(resfile+'.tax',index=False)
    c2size = {}
    for i in clusters:
        if i in c2size.keys():
            c2size[i] += 1
        else:
            c2size[i] = 1
    sensitivity = 0
    

    return clusters,classes

parser = argparse.ArgumentParser()
parser.add_argument("resfile",type=str)
parser.add_argument("tax",type=str)
args = parser.parse_args()

clusters,classes = c1res(args.resfile,args.tax)

print('RI: ',rand_index_score(clusters,classes) )       
print('ARI: ',metrics.adjusted_rand_score(classes, clusters))
    

