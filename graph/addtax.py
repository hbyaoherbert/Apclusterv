import pandas as pd
import argparse
import ast
parser = argparse.ArgumentParser()
parser.add_argument('cluster',type=str)
parser.add_argument('tax',type=str)
args = parser.parse_args()

cldf = pd.read_csv(args.cluster,sep=',',header=0,index_col=0)
rlist = ['genus','family','order','class','phylum']
dbrlist = ['species','genus','family','order','class','phylum']

taxdf = pd.read_csv(args.tax,sep='\t',header=0)
ctgkey = taxdf['contig'].map(lambda d:d.replace(' ','~'))
ctg2g = dict(zip(ctgkey,taxdf['genus']))



df = cldf.copy()
taxcol = []
memcol = []
for idx,row in cldf.iterrows():
    for level in ['genus']:
        #taxlist = ast.literal_eval(row[level])
        mems = []
        taxset = set()
        tax2num = {}
        total = 0
        for member in row['Members'].split(' '):
            
            if not member in ctg2g.keys():
                continue
            tax = ctg2g[member]
            if tax in tax2num.keys():
                tax2num[tax] += 1
            else:
                tax2num[tax] = 1
            total += 1
            mems.append(member+'#'+ctg2g[member])
             
            
            taxset.add(ctg2g[member])
        taxcol.append(str(list(taxset)))
        if len(tax2num) == 0:
            continue
        tax2num = list(tax2num.items())
        tax2num = sorted(tax2num,key=lambda d:d[1],reverse=True)
        if tax2num[0][1]/total < 0.7:
            print(idx,mems,taxset)
    
        

df['genus'] = taxcol

df.to_csv('c1.clusters.tax',sep=',',index=True)

