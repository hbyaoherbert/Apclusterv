from apclusterv.cluster.calscore import *
from apclusterv.cluster.runmcl import *
from apclusterv.newscore.ratiosim import *
from apclusterv.cluster.ctgconnect import *
from apclusterv.cluster.runap import *
from apclusterv.cluster.mergecluster import *

class Arg:
	def __init__(self):
		pass

def maincluster(args):
	#parser = argparse.ArgumentParser()
	#parser.add_argument('map',type=str)
	#parser.add_argument('tab',type=str)
	#parser.add_argument('prot',type=str)
	
	
	#args = parser.parse_args()
	
	
	print("calculate contig-contig alignment and ratio score")
	aggrctgscore(args)
	
	print("generate input for mcl clustering")
	
	gsimargs = Arg()
	setattr(gsimargs, 'map', args.map)
	setattr(gsimargs,'aggr','tmp/ctgscore.csv')
	setattr(gsimargs,'prot',args.prot)
	setattr(gsimargs,'out','tmp/shared_protein.csv')
	
	calglobal(gsimargs)
	
	print("run mcl clustering")
	mclargs = Arg()
	inflation = args.inflation
	setattr(mclargs, 'edges', 'tmp/shared_protein.csv')
	setattr(mclargs,'inflation',inflation)
	run_mcl(mclargs)
	
	
	print("generate input for affinity propagation")
	
	
	apargs = Arg()
	
	setattr(apargs, 'map', args.map)
	setattr(apargs,'alnscore','tmp/ctgscore.csv')
	setattr(apargs,'pcscore','tmp/shared_protein.csv')


	#parser.add_argument("resfile",type=str)
	#parser.add_argument("aggr",type=str)
	#parser.add_argument("map",type=str)
	#filter_edge(apargs)
	#construct_graph(apargs.aggr+".abs",apargs.aggr+".rep")
	
	
	connect(apargs)
	
	print("run affinity propagation")
	runapargs = Arg()
	setattr(runapargs,'edgefile','tmp/shared_protein.csv.aln')
	setattr(runapargs,'repfile','tmp/shared_protein.csv.rep')
	setattr(runapargs,'apk',args.apk)
	runcluster(runapargs)
	
	
	print("merge cluster")
	
	exargs = Arg()
	setattr(exargs,'alnfile','tmp/dvp.ctg.aln.tab')
	setattr(exargs,'apfile','tmp/c2members')
	setattr(exargs,'apscore','tmp/ctgscore.csv')
	setattr(exargs,'map',args.map)
	setattr(exargs,'mcl','tmp/mcl'+str(args.inflation))
	setattr(exargs,'inflation',args.inflation)
	setattr(exargs,'apk',args.apk)

	prepare(exargs.alnfile,exargs.apfile,exargs.apscore)

	print("export result")
	finalres(exargs)
	#prepare(args.alnfile,args.apfile,args.apscore)
