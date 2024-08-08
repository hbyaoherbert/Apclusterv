import pandas as pd
import argparse
import json

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
	
	return score,number,answer[0]
	
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
def exportres(member2new,member2mcl):
	cl2members = {}
	for member,cl in member2mcl.items():
		cl2members.setdefault(cl,[])
		cl2members[cl].append(member)
	new2members = {}   
	
	
	for member,new in member2new.items():
		new2members.setdefault(new,[])
		new2members[new].append(member)
	
	data = []
	for cl,members in cl2members.items():
		nosupport = []
		supported = []
		
		
		for member in members:
			if member in new2members.keys():
				supported.append(member)
			else:
				nosupport.append(member)
				
		newlist = []
		newlist += supported
		for new in supported:
			newlist += new2members[new]        
				
		print(len(supported),len(nosupport))
		if len(newlist)>len(nosupport) and len(supported)*3 < len(nosupport):
			print("newcluster")
			oldstr = ",".join(nosupport)
			data.append([len(nosupport),oldstr])
			
				 
			print(len(newlist),len(nosupport))
			newstr = ",".join(newlist)
			data.append([len(newlist),newstr])
		else:
			for new in supported:
				 members += new2members[new]
			memberstr = ",".join(members)
			data.append([len(members),memberstr])
	 
	resultdf = pd.DataFrame(data,columns=['Size','Members'])
	resultdf.to_csv('publish_apres.csv')

def finalres(args)
	apfile = open(args.c2members,'r')
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




	examdf = pd.read_csv(args.exam,header=None,sep=' ',skiprows=1)
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
		if ctgname1 in member2mcl.keys():
			continue
		nprot = ctg2nprot[ctg1]
		if nprot < 2:
			 continue
			
		if ctg1 in member2c.keys():
			continue
		
		
		#ctax = "unsign"
		#if member2c[ctg1] in id2tax.keys():
			#ctax = id2tax[member2c[ctg1]]
		#print(ctg1,tax1)
		score,number,ctg2 = assign(ctg1,ctg2anno[ctg1])
		ctg2name = id2ctg[ctg2]
		ctg1name = id2ctg[ctg1]
		
		if ctg2name in member2mcl.keys():
			if number/nprot >0.5:
	   
				member2new[ctg1name] = ctg2name
				
				
				admit +=1

		 

	
	print(admit)    
	
	#for member,mcl in member2mcl.items():
		#if not member in member2new.keys():
			#member2new[member] = mcl
	exportres(member2new,member2mcl)
	
