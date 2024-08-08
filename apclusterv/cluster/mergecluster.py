import pandas as pd
import numpy as np
import json
import argparse



def filter(examdf,col,supported):
	ctg2anno = {}
	ctg2annonum = {}
	prot2best = {}
	
	#Notice that an ap component is a cluster of proteins. Now we will aggregate protein-protein similarities to the genome-genome level.
	
	outdata = []
	outedgedata = []
	for qprot,tprot,score,ctg1,ctg2 in zip(examdf['qprot'],examdf['tprot'],examdf[col],examdf['ctg1'],examdf['ctg2']):
		if not tprot in supported:
			continue
		
		if (qprot,ctg2) in prot2best.keys():
			continue
		outedgedata.append([qprot,tprot,score,ctg1,ctg2])
		prot2best[(qprot,ctg2) ] = score
		ctg2anno.setdefault(ctg1,{})
		ctg2anno[ctg1].setdefault(ctg2,0)
		
		ctg2annonum.setdefault(ctg1,{})
		ctg2annonum[ctg1].setdefault(ctg2,0)
		
		ctg2anno[ctg1][ctg2] += score
		ctg2annonum[ctg1][ctg2] += 1
		
	for ctg1,anno in ctg2anno.items():
		annonum = ctg2annonum[ctg1]
		for ctg2,score in anno.items():
			
			outdata.append([ctg1,ctg2,score,annonum[ctg2]])

	outedgedf = pd.DataFrame(outedgedata,columns = ['qprot','tprot','score','ctg1','ctg2'])
	outedgedf.to_csv('examedge.csv',index=False)
	outdf = pd.DataFrame(outdata,columns = ['ctg1','ctg2','score','number']) 
	outdf.to_csv('tmp/exam',index=False,sep=',')
		
def pick(alnfile,tset,pairset,t2cut,member2ap):
	df = pd.read_csv(alnfile,sep='\t')
	
	df = df[df['ctg2'].isin(tset)]
	
	#print(len(df))
	whitelist = {}
   
	examlist = []    
	supported =  set()
	for row in df.itertuples():
		if row[7] == row[8]:
			continue
		pairkey = (row[7],row[8])
		if not int(row[8]) in tset:
			continue
		cutoff = t2cut[int(row[8])]
		ratio = int(row[3])*1.0/( int(row[3]) + int(row[4]) + int(row[5]))
		#if row[2] == 101982:
			#print(ratio,cutoff)
			
		#we will include genomes that haven't been assiged to any clusters to the waiting list, as long as they can align to a component center with similarity over its cut-off.
		if not pairkey in pairset:
			if ratio >= cutoff:
			#if True:
				supported.add(row[2])
				examlist.append([row[1],row[2],row[3],row[6],row[7],row[8]])
			continue
		
		
		if ratio >= cutoff:
		#if True:    
			supported.add(row[2])
			#if int(row[6]) == 15046:
				#print(row[7],row[1],row[2],row[3],row[6])
				
			#whitelist:for each ap center, all alignments over the component-wise cut-off will support the center as a confident representative
			if pairkey in whitelist.keys():    
				whitelist[pairkey].append([row[1],row[2],row[3]])
			else:
				whitelist[pairkey] = [[row[1],row[2],row[3]]]
			
	outdata = []
	for pairkey,alns in whitelist.items():
		ctg1,ctg2 = pairkey
		for prot1,prot2,match in alns:
			outdata.append([ctg1,ctg2,prot1,prot2,match,])
			
	outdf = pd.DataFrame(outdata,columns=['query','target','qprot','tprot','match'])
	outdf = outdf.sort_values(by=['qprot','match'],ascending=[True,False])
	outdf.to_csv('whitelist.csv',index=False)        
			
	examdf = pd.DataFrame(examlist,columns = ['qprot','tprot','match','score','ctg1','ctg2'])
	
	filter(examdf,'score',supported)
	
		
def prepare(alnfile,apfile,scorefile):
	
	
	
		
	ap = open(apfile,'r')
	c2members = json.load(ap)           
		   
	scoredf = pd.read_csv(scorefile,sep='\t')
	ap2members = {}
	ap2scores = {}
	member2ap = {}
	for c,members in c2members.items():
		ap2members[int(c)] = members
		ap2scores[int(c)] = []
		for member in members:
			member2ap[member] = int(c)
		
	for ctg1,ctg2,score in zip(scoredf['ctg1'],scoredf['ctg2'],scoredf['ratioscore']):
		if ctg1 == ctg2:
			continue
		if ctg1 in ap2members.keys():
			ap2scores[ctg1].append(score)
		if ctg2 in ap2members.keys():
			ap2scores[ctg2].append(score)
	for ap in ap2scores.keys():
		ap2scores[ap] = np.mean(ap2scores[ap])
		
	#for each ap compoment, the mean values of edges to the representive node is regarded as the cut-off for this component.
	pairset = set()
	for ap,members in ap2members.items():
		for member in members:
			if member == ap:
				continue
			pairkey = (member,ap)
			pairset.add(pairkey)
		
	pick(alnfile,ap2members.keys(),pairset,ap2scores,member2ap)   
			

#parser = argparse.ArgumentParser()
#parser.add_argument('alnfile',type=str)
#parser.add_argument('apfile',type=str)
#parser.add_argument('scorefile',type=str)
#args = parser.parse_args()
#prepare(args.alnfile,args.apfile,args.apscore)




#parser = argparse.ArgumentParser()
#parser.add_argument('map',type=str)

#parser.add_argument('exam',type=str)
#parser.add_argument('c2members',type=str)
#parser.add_argument('mcl',type=str)
#args = parser.parse_args()
#id 16192
#a previously unassigned genome (query) is now assigned to the closet genome (target), by aggregating the bitscores of alignments between a protein on the query genome and a protein in the APC set on the target genome. To make such assignment more confident, we require that at least 50% of all the proteins on the query genome can align to the target genome [].  
def assign(ctg1,anno):
	anno = sorted(anno,key=lambda x:x[1],reverse=True)
	answer = anno[0]
	score = answer[1]
	number = answer[2]
	
	advantage = -1
	second_target = -1 
	if len(anno)>1:
		second = anno[1]
		advantage = score - second[1]
		second_target = second[0]
	return score,number,answer[0],second_target,advantage
	
'''
(vcontact2env) yaohaobin@ubuntu:~/share/viral/dvp$ python exec/eval.py newscore/publish_apres.csv fulltax.tab
sensitivity: 6551
108505.0 21346968.0 65075.0 45480.0
0.9948736503541589
0.6599209898811542
(vcontact2env) yaohaobin@ubuntu:~/share/viral/dvp$ python exec/eval.py mcl3 fulltax.tab
sensitivity: 4236
53459.0 8937893.0 26228.0 24298.0
0.9944120015775484
0.6762677859891169
if only a small portion of genomes in a step-2 cluster are supported by a large number of newly-admitted genomes, we will extract the supported minorites from the cluster and let them form a new cluster together with the newly-admitted ones. More specifically, we will split the step-2 genome cluster when supported genomes is less than k=1/4, while the new cluster will be larger than the remaining old cluster. 

'''
def exportres(member2new,member2mcl,ctg2advantage,admitdf,id2ctg,ctg2id,args):
	cl2members = {}
	for member,cl in member2mcl.items():
		cl2members.setdefault(cl,[])
		cl2members[cl].append(member)
	new2members = {}   
	
	
	for member,new in member2new.items():
		new2members.setdefault(new,[])
		new2members[new].append(member)
	
	data = []

	admitdata = []
	for cl,members in cl2members.items():
		pruneset = judgecluster(ctg2advantage,cl,admitdf,args.apk)
		for ctgid in pruneset:
			ctgname = id2ctg[ctgid]
			admits_ind = new2members[ctgname]
			data.append([len(admits_ind),','.join(admits_ind)])

		nosupport = []
		supported = []
		
		outmembers = []

		for member in members:
			
			if member in new2members.keys():
				ctgid = ctg2id[member]
				if ctgid in pruneset:
					continue
				supported.append(member)
				admitnum = len(new2members[member])
				outmembers.append(member+"#"+str(admitnum))
			else:
				outmembers.append(member+"#0")
				nosupport.append(member)
		memberstr = ",".join(outmembers)
		admitdata.append([cl,memberstr])

		newlist = []
		newlist += supported
		for new in supported:
			newlist += new2members[new]        
				
		#print(len(supported),len(nosupport))
		if len(supported)>=1 and len(supported)*3 < len(nosupport):
		#if len(newlist)>len(nosupport) and len(supported)*3 < len(nosupport):
		#if False:
			#print("newcluster")
			oldstr = ",".join(nosupport)
			data.append([len(nosupport),oldstr])
			
				 
			#print(len(newlist),len(nosupport))
			newstr = ",".join(newlist)
			data.append([len(newlist),newstr])
		else:
			for new in supported:
				 members += new2members[new]
			memberstr = ",".join(members)
			data.append([len(members),memberstr])
	mergedf = pd.DataFrame(admitdata,columns = ['cl','Members'])
	mergedf.to_csv('mergestat.csv',index=False)	 
	resultdf = pd.DataFrame(data,columns=['Size','Members'])
	resultdf.to_csv('tmp/cluster_result.'+str(args.inflation)+'.'+str(args.apk)+'.csv')

def judgecluster(ctg2advantage,querycl,admitdf,apk):
	admit2advantage = {}
	for ctg1,ctg2,cl in zip(admitdf['ctg1'],admitdf['ctg2'],admitdf['cl']):
		if not ctg1 in ctg2advantage.keys():
			continue
		if not ctg2 in ctg2advantage.keys():
			continue
		if querycl != cl:
			continue
		thisadvantage = ctg2advantage[ctg1][2]
		if thisadvantage <0:
			continue
		admit2advantage.setdefault(ctg2,[])
		admit2advantage[ctg2].append(thisadvantage)
		center_best = ctg2advantage[ctg2][0]
		
		#print(ctg1,ctg2,thisadvantage,center_best)
	pruneset = set()
	for admitter in admit2advantage.keys():
		#print(admitter,np.mean(admit2advantage[admitter]),ctg2advantage[admitter])
		center_best = ctg2advantage[admitter][0]
		if center_best in admit2advantage.keys():
			partner_ratio = np.mean(admit2advantage[center_best])
			this_ratio = np.mean(admit2advantage[admitter])
			if partner_ratio + this_ratio and this_ratio<partner_ratio > float(apk):
				#print('prune',admitter,center_best,this_ratio,partner_ratio,querycl)
				pruneset.add(admitter)
	return pruneset
		

def finalres(args):
	apfile = open(args.apfile,'r')
	c2members = json.load(apfile)
	ap2members = {}
	member2c = {}
	for c,members in c2members.items():
		ap2members[int(c)] = members
		for member in members:
			member2c[member] = int(c)
		
	mapdf = pd.read_csv(args.map,sep=',',header=0)

	ctg2id = dict(zip(mapdf['contig'],mapdf['ctgid']))
	id2ctg = dict(zip(mapdf['ctgid'],mapdf['contig']))

	ctg2nprot = {}
	for prot,ctg in zip(mapdf['protein_id'],mapdf['ctgid']):
		ctg2nprot.setdefault(ctg,0)
		ctg2nprot[ctg] += 1




	examdf = pd.read_csv('tmp/exam',header=0,sep=',')
	examdf.columns = ['ctg1','ctg2','score','number']

	ctg2anno = {}
	for ctg1,ctg2,score,number in zip(examdf['ctg1'],examdf['ctg2'],examdf['score'],examdf['number']):

	#print(ctg1,ctg2,score,tax1,tax2)

		if ctg1 in ctg2anno.keys():
			ctg2anno[ctg1].append([ctg2,score,number])
		else:
			ctg2anno[ctg1] = [[ctg2,score,number]]

	mcldf = pd.read_csv(args.mcl,sep=',',header=0)
	member2mcl = {}

	member2new = {}
	mclnum = 0
	admitdata = []

	incluster = []
	for idx,row in mcldf.iterrows():
		mclnum += 1
		members = row['Members'].split(',')
		for member in members:
			member2mcl[member.replace(' ','~')] = idx

	ncut = 20  
	p = 0
	n = 0
	np = 0
	nomcl = 0

	pscore = 0
	npscore = 0
	
	admit = 0
	for ctg1 in ctg2anno.keys():
		
		ctgname1 = id2ctg[ctg1]
		
			
		
		
		
		#ctax = "unsign"
		#if member2c[ctg1] in id2tax.keys():
			#ctax = id2tax[member2c[ctg1]]
		#print(ctg1,tax1)
		score,number,ctg2,second_target,advantage = assign(ctg1,ctg2anno[ctg1])
		incluster.append([ctg1,ctg2,second_target,advantage])



		if ctg1 in member2c.keys():
			continue
		ctg2name = id2ctg[ctg2]
		ctg1name = id2ctg[ctg1]
		if ctgname1 in member2mcl.keys():
			continue
		nprot = ctg2nprot[ctg1]
		if nprot < 2:
			 continue
		if ctg2name in member2mcl.keys():
			if number/nprot >0.5:
				nprot2 = ctg2nprot[ctg2]
				member2new[ctg1name] = ctg2name
				admitdata.append([ctg1,ctg2,member2mcl[ctg2name],number,nprot,nprot2,advantage])
				
				admit +=1

		 

	
	#print(admit)    
	
	#for member,mcl in member2mcl.items():
		#if not member in member2new.keys():
			#member2new[member] = mcl
	
	incluster_df = pd.DataFrame(incluster,columns=['ctg1','ctg2','second','advantage'])
	incluster_df.to_csv('inclusternum.csv',index=False)
	ctg2advantage = {}
	for ctg1,ctg2,second,advantage in zip(incluster_df['ctg1'],incluster_df['ctg2'],incluster_df['second'],incluster_df['advantage']):
		ctg2advantage[ctg1] = [ctg2,second,advantage]

	admitdf = pd.DataFrame(admitdata,columns=['ctg1','ctg2','cl','number','nprot1','nprot2','advantage'])
	admitdf.to_csv('admitpair.csv',index=False)
	exportres(member2new,member2mcl,ctg2advantage,admitdf,id2ctg,ctg2id,args)

