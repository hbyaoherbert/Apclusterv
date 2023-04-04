import pandas as pd 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("alnfile",type=str)
parser.add_argument("edf",type=str)
args = parser.parse_args()


df = pd.read_csv(args.alnfile,sep='\t',header=None)
edf = df[[0,1,3,10,11]]
pset = set()
data = []
for prot1,prot2,pid,e,score in zip(edf[0],edf[1],edf[3],edf[10],edf[11]):
   
    #p1 = min(prot1,prot2)
    #p2 = max(prot1,prot2)
    #key = p1+'#'+p2
    #if key in pset:
        #continue
    #pset.add(key)
    data.append([prot1,prot2,pid,e,score])
df = pd.DataFrame(data,columns=['prot1','prot2','hsp','e','score'])

df.to_csv(args.edf,sep='\t',index=False)

