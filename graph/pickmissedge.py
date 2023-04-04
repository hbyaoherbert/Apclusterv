import pandas as pd 
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("pcfile",type=str)
parser.add_argument("missgenus",type=str)
parser.add_argument("tabfile",type=str)
args = parser.parse_args()

missgenus = pd.read_csv(args.missgenus,sep=',',skiprows=1,header=None)
missgenus.columns = ['protein_id','contig_id','pc','genus']
pcset = set(missgenus['pc'])

pcdf = pd.read_csv(args.pcfile,sep=',',header=0)
pcdf = pcdf.fillna(-1)
pc_ctg = [[
for protein,contig,pc in zip(pcdf['protein_id'],pcdf['contig_id'],pcdf['cluster']):
    if pc in pcset:
        pc_ctg.append([protein,contig])
    
