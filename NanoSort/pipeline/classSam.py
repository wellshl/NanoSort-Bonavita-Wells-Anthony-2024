import config
import numpy as np
import pandas as pd
import pickle
from . import classAln

with open(config.pkl1, 'rb') as pkl:
	pos_to_p1 = pickle.load(pkl)
	p1_to_pos = {p1:pos for pos,p1 in pos_to_p1.items()}

with open(config.pkl2, 'rb') as pkl:
	pos_to_p2 = pickle.load(pkl)
	p2_to_pos = {p2:pos for pos,p2 in pos_to_p2.items()}

with open(config.bed1) as bed:
	b1 = pd.read_table(bed,delimiter='\t',header=0)
	b1_TRS = []
	b1_TRS_B_index = list(b1.index[b1["region"]=="TRS-B"])
	b1_TRS_L_index = list(b1.index[b1['region']=="TRS-L"])
	b1_TRS_index = b1_TRS_B_index + b1_TRS_L_index
	for i in b1_TRS_index:
		b1_TRS = b1_TRS + list(range(b1["start"][i],b1["end"][i]+1))

with open(config.bed2) as bed:
	b2 = pd.read_table(bed,delimiter='\t')
	b2_TRS = []
	b2_TRS_B_index = list(b2.index[b2["region"]=="TRS-B"])
	b2_TRS_L_index = list(b2.index[b2["region"]=="TRS-L"])
	b2_TRS_index = b2_TRS_B_index + b2_TRS_L_index
	for i in b2_TRS_index:
		b2_TRS = b2_TRS + list(range(b2["start"][i],b2["end"][i]+1))

with open(config.cutoffs) as csv:
	cuts = pd.read_table(csv, delimiter=',',dtype={0:int,1:float,1:float},header=0, names=['pos','p2_below','p1_above'])

s = config.s
w = config.w
p_count = config.p_count

class Sam:
	def __init__(self, file):
		self.filename = file
		colnames = ['read_name','flag','parent','pos','unk','cigar','star','first_0','second_0','sequence','quality','tag_NM','tag_ms','tag_AS','tag_nn','tag_tp','tag_cm','tag_s1','tag_s2','tag_de','tag_SA','tag_rl']
		self.sam_sliding = pd.read_table(f'{file}_sliding.sam',header=None,names=colnames,quoting=3,engine='python')
		self.sam = pd.read_table(f'{file}.sam',header=None,comment='@',names=colnames,quoting=3,engine='python')

		self.parent1 = config.parent1
		self.parent2 = config.parent2
		self.read_name = self.sam_sliding.loc[0,'read_name'].split("_sliding")[0]
		self.p1_start_list = []
		self.like_list = []
		self.classify_list = []
		self.recomb = 0
		self.subgen = 0
		self.nonhom = 0

	def init_pos(self, i):
	# initialize the starting position of the alignment
		self.mapped_parent = self.sam_sliding.loc[i,'parent']
		position = self.sam_sliding.loc[i,'pos']
		if self.mapped_parent==config.parent1:
			p1_pos = position
			self.p1_start_list = self.p1_start_list + [p1_pos]
			pos = p1_to_pos[p1_pos]
		if self.mapped_parent==config.parent2:
			p2_pos = position
			pos = p2_to_pos[p2_pos]
			if pos==0:
				p1_pos = 1
			else:
				p1_pos = pos_to_p1[pos] + 1
			self.p1_start_list = self.p1_start_list+[p1_pos]

		return(pos)

	def extend_cigar(self, pos, i):
	# extend the alignment based on the cigar string and calculate likelihood
		read = self.sam_sliding.loc[i,'sequence']
		cigar = self.sam_sliding.loc[i,'cigar']
		if self.mapped_parent==config.parent1:
			pos = p1_to_pos[self.sam_sliding.loc[i,'pos']]
		else:
			pos = p2_to_pos[self.sam_sliding.loc[i,'pos']]

		aligned = classAln.Aln(self.mapped_parent,read,cigar,pos)
		classAln.Aln.align(aligned)
		self.like_list = self.like_list + [aligned.loglike1-aligned.loglike2]

		return(aligned)

	def classify_subread(self, i):
	# classify subread as one parent or the other based on likelihood score 
		p1_min_cut = cuts.query(f'pos=={self.p1_start_list[i]}')['p1_above'].values
		p2_max_cut = cuts.query(f'pos=={self.p1_start_list[i]}')['p2_below'].values

		if self.like_list[i] > p1_min_cut:
			self.classify_list = self.classify_list + [config.parent1]

		elif self.like_list[i] < p2_max_cut:
			self.classify_list = self.classify_list + [config.parent2]

		else:
			self.classify_list = self.classify_list + ["unclassified"]

	def nonhom_subgen(self):
	# classify whole read as subgenomic or nonhomologous based on deletions
		primary = self.sam.loc[self.sam['flag']==0]
		if primary['parent'].values[0]==config.parent1:
			pos = p1_to_pos[primary['pos'].values[0]]
		else:
			pos = p2_to_pos[primary['pos'].values[0]]
		aligned = classAln.Aln(primary['parent'].values[0], primary['sequence'].values[0], primary['cigar'].values[0], pos)

		# primary starts at TRS-B but is trimmed ~90odd bp
		threshold = 20
		canonical_TRS_start = [el for i in b1_TRS_index[:-1] for el in [*range(b1['start'][i] - threshold, b1['start'][i] + threshold)]]
		possible_TRSL = [*range(b1['start'][b2_TRS_L_index].values[0] - threshold, b1['start'][b2_TRS_L_index].values[0] + threshold)]

		# reads start at TRS-B and leading soft/hard trim over TRS-L
		if pos_to_p1[pos] in canonical_TRS_start:
			if aligned.cigar_do[0] in ['S','H']:
				if aligned.cigar_count[0] in possible_TRSL:
					self.subgen = 1
		
		# leading soft trim over TRS-L
		elif aligned.cigar_do[0] in ['S','H']: #non-canonical
			if aligned.cigar_count[0] in possible_TRSL:
				self.subgen = 1
				
		# look for supplementary alignments starting at TRS-L with trailing hard trim
		supp = self.sam.loc[self.sam['flag']==2048]
		if not supp.empty:
			if primary['parent'].values[0]==config.parent1:
				supp_pos = p1_to_pos[supp['pos'].values[0]]
			else:
				supp_pos = p2_to_pos[supp['pos'].values[0]]
			supp_aligned = classAln.Aln(supp['parent'].values[0], primary['sequence'].values[0], supp['cigar'].values[0], supp_pos)

			if supp_aligned.cigar_do[-1]=='H':
				if supp_aligned.align()[-1] in possible_TRSL:
					if pos_to_p1[aligned.align()[1]] in canonical_TRS_start:
						self.subgen = 1
						
					else: #non-canonical
						self.subgen = 1

			# primary has a large soft trim and supp has a large hard trim 100 bp apart
			if self.subgen==0:
				if aligned.cigar_do[0]=='S':
					if aligned.cigar_count[0] > 100:
						if supp_aligned.cigar_do[-1]=='H':
							if supp_aligned.cigar_count[-1] > 100:
								self.nonhom = 1
								
			if aligned.cigar_do[-1]=='S':
				if aligned.cigar_count[-1] > 100:
					if supp_aligned.cigar_do[0]=='H':
						if supp_aligned.cigar_count[0] > 100:
							self.nonhom = 1

		## OR ##
		# deletion of size >100 bp
		dels=[i for i in range(len(aligned.cigar_do)) if aligned.cigar_do[i]=='D']
		where=np.argwhere(np.array(aligned.cigar_count)[dels]>100)
		if where:
			self.nonhom = 1


		






