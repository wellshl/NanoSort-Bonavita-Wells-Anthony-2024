import argparse
import config
import json
from pipeline import classSam
import numpy as np

parser=argparse.ArgumentParser(description='HTS')
parser.add_argument('-i','--inp',required=True)
parser.add_argument('-o','--out',required=True)
args=parser.parse_args()

Sam=classSam.Sam(f'{args.inp}')

s=config.s
w=config.w
p_count=config.p_count

def parse_sam_for_nonhom(Sam):

	classSam.Sam.nonhom_subgen(Sam)

def parse_sam_for_recombinants(Sam):

	for i in range(0,len(Sam.sam_sliding)):

		if Sam.sam_sliding.loc[i,'flag']!=0:
			Sam.p1_start_list=Sam.p1_start_list+["unmapped"]
			Sam.like_list=Sam.like_list+["unmapped"]
			Sam.classify_list=Sam.classify_list+["unmapped"]
			continue

		pos=classSam.Sam.init_pos(Sam,i)
		classSam.Sam.extend_cigar(Sam,pos,i)
		classSam.Sam.classify_subread(Sam,i)

parse_sam_for_recombinants(Sam)
parse_sam_for_nonhom(Sam)

p1_total = Sam.classify_list.count(config.parent1)
p2_total = Sam.classify_list.count(config.parent2)
if sum([(p1_total>=p_count),(p2_total>=p_count)])==2:
	Sam.recomb=1
	r=str(args.out)+"/recombinant_read_names.txt" 
	with open(r,'a') as recs:
		recs.write(Sam.read_name+"\n")

if Sam.nonhom==1:
	n=str(args.out)+"/nonhomologous_read_names.txt"
	with open(n,'a') as nonhoms:
		nonhoms.write(Sam.read_name+"\n")

if Sam.subgen==1:
	s=str(args.out)+"/subgenomic_read_names.txt"
	with open(s,'a') as subs:
		subs.write(Sam.read_name+"\n")	

outfile=str(args.out)+"/overall_summary.tsv"
with open(outfile,'a') as summ:
	summ.write("\t".join([Sam.read_name,str(Sam.recomb),str(Sam.nonhom),str(Sam.subgen)])+"\n")

report={
	'read_name': Sam.read_name,
	'position': Sam.p1_start_list,
	'likelihood': Sam.like_list,
	'classification': Sam.classify_list
}
jsonfile=f'{args.out}/reads/{Sam.read_name}.json'
with open(jsonfile, 'w') as out_f:
	json.dump(report, out_f, indent=4, default=int)