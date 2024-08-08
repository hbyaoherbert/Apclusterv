import pandas as pd

import numpy as np

#parser = argparse.ArgumentParser()
#parser.add_argument('map',type=str)
#parser.add_argument('tab',type=str)
#args = parser.parse_args()

def aggrctgscore(args):
	mapdf = pd.read_csv(args.map,sep=',',header=0)
	prot2id = dict(zip(mapdf['protein_id'],mapdf.index))
	ctg2id = dict(zip(mapdf['contig'],mapdf['ctgid']))
	prot2ctg = dict(zip(mapdf.index,mapdf['ctgid']))

	df = pd.read_csv(args.tab,sep='\t',header=None)
	df = df[[0,1,3,4,5,11]]

	df.columns = ['prot1','prot2','match','mis','gap','alnscore']
	id1series = df['prot1'].map(lambda x:prot2id[x])
	id2series = df['prot2'].map(lambda y:prot2id[y])

	df['prot1'] = id1series
	df['prot2'] = id2series

	df['ctg1'] = df['prot1'].map(lambda x:prot2ctg[x])
	df['ctg2'] = df['prot2'].map(lambda y:prot2ctg[y])
	#print(df.head(100))
	df.to_csv('tmp/dvp.ctg.aln.tab',index=False,sep='\t')


	newdf = df.groupby(['ctg1','ctg2']).agg({'match':[np.sum],'mis':[np.sum],'gap':[np.sum],'alnscore':[np.sum]}).reset_index()
	newdf.columns = ['ctg1','ctg2','match','mis','gap','alnscore']

	key2alnscore = {}
	key2ratioscore = {}
	for row in newdf.itertuples():
		ctg1 = row[1]
		ctg2 = row[2]
		match = row[3]
		mis = row[4]
		gap = row[5]
		alnscore = row[6]		

		ratioscore = match/(match+mis+gap)
		ctga = min(ctg1,ctg2)
		ctgb = max(ctg1,ctg2)
		key = "{}#{}".format(ctga,ctgb)
		if key in key2alnscore.keys():
			key2alnscore[key] = max(key2alnscore[key],alnscore)
			key2ratioscore[key] = max(key2ratioscore[key],ratioscore)

		else:
			key2alnscore[key] = alnscore
			key2ratioscore[key] = ratioscore

	data = []
	for key,alnscore in key2alnscore.items():
		ratioscore = key2ratioscore[key]
		info = key.split('#')
		ctg1 = info[0]
		ctg2 = info[1]
		data.append([ctg1,ctg2,alnscore,ratioscore])

	outdf = pd.DataFrame(data,columns=['ctg1','ctg2','alnscore','ratioscore'])

	outdf.to_csv('tmp/ctgscore.csv',sep='\t',index=False)


