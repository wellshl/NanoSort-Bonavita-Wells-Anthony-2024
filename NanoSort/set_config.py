import argparse
from datetime import date

parser=argparse.ArgumentParser(description='HTS')
parser.add_argument('-i','--inputreads',required=True,type=str)
parser.add_argument('-p1','--parent1',required=True,type=str)
parser.add_argument('-p2','--parent2',required=True,type=str)
parser.add_argument('-pc','--p_count',required=False,default=41,type=int)
parser.add_argument('-s','--step',required=True,type=int)
parser.add_argument('-w','--window',required=True,type=int)
parser.add_argument('-wd','--workingdir',required=True,type=str)
parser.add_argument('-nd','--nanosortdir',required=True,type=str)
args=parser.parse_args()

with open(args.workingdir+"/config.py",'w') as in_handle:
	in_handle.write(
		f'input_reads = "{args.inputreads}"\n'
		f'date_of_run = "{date.today().strftime("%m/%d/%Y")}"\n'
		f'parent1 = "{args.parent1}"\n'
		f'parent2 = "{args.parent2}"\n'
		f's = {args.step}\n'
		f'w = {args.window}\n'
		f'p_count = {args.p_count}\n'
		f'fasta1 = "{args.nanosortdir}/data/singles/{args.parent1}/{args.parent1}.fasta"\n'
		f'fasta2 = "{args.nanosortdir}/data/singles/{args.parent2}/{args.parent2}.fasta"\n'
		f'count1 = "{args.nanosortdir}/data/singles/{args.parent1}/{args.parent1}_count_matrix.tsv"\n'
		f'count2 = "{args.nanosortdir}/data/singles/{args.parent2}/{args.parent2}_count_matrix.tsv"\n'
		f'pkl1 = "{args.nanosortdir}/data/singles/{args.parent1}/{args.parent1}.pkl"\n'
		f'pkl2 = "{args.nanosortdir}/data/singles/{args.parent2}/{args.parent2}.pkl"\n'
		f'bed1 = "{args.nanosortdir}/data/singles/{args.parent1}/{args.parent1}.bed"\n'
		f'bed2 = "{args.nanosortdir}/data/singles/{args.parent2}/{args.parent2}.bed"\n'
		f'cutoffs = "{args.nanosortdir}/data/pairs/{args.parent1}_{args.parent2}/{args.parent1}_{args.parent2}_cutoffs.csv"'
		)
	in_handle.close()
