import pandas as pd

df1 = pd.read_csv('cltest/aggrscore.csv',sep='\t',header=0)
df2 = pd.read_csv('cltest/aggrtest.csv',sep='\t',header=0)

cset = set()
for ctg1,ctg2 in zip(df1['ctg1'],df1['ctg2']):
    cset.add((ctg1,ctg2))

for ctg1,ctg2 in zip(df2['ctg1'],df2['ctg2']):
    if not (ctg1,ctg2) in cset:
        if not (ctg2,ctg1) in cset:
            print(ctg1,ctg2)





