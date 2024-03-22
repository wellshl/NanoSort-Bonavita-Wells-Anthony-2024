import config
import pickle
import re
from . import likelihood

with open(config.pkl1,'rb') as pkl:
	pos_to_p1 = pickle.load(pkl)
	p1_to_pos = {p1:pos for pos,p1 in pos_to_p1.items()}

with open(config.pkl2,'rb') as pkl:
	pos_to_p2 = pickle.load(pkl)
	p2_to_pos = {p2:pos for pos,p2 in pos_to_p2.items()}

def add_c_to_pos(pos, c, parent):
	if parent==config.parent1:
		try:
			return p1_to_pos[pos_to_p1[pos] + c]
		except:
			return max(p1_to_pos.values())
	else:
		try:
			return p2_to_pos[pos_to_p2[pos] + c]
		except:
			return max(p2_to_pos.values())

class Aln:

	def __init__(self,parent, read, cigar, pos):
		self.mapped_parent = parent
		self.read = read
		self.cigar = cigar
		self.cigar_do = re.findall('[A-Z]', self.cigar)
		self.cigar_count = [int(cc) for cc in re.findall('[0-9]+',self.cigar)]
		self.start_pos = pos
		self.pos = ['NA']*(len(self.cigar_do)+1)
		self.read_pos = ['NA']*(len(self.cigar_do)+1)
		self.loglike1 = 0.
		self.loglike2 = 0.

	def align(self):
		
		read = list(self.read)
		cigar_count = self.cigar_count
		cigar_do = self.cigar_do

		if self.start_pos==0:
			self.pos[0] = 1
			self.read_pos[0] = 1
			cigar_count[0] = cigar_count[0]-1
		else:
			self.pos[0] = self.start_pos
			self.read_pos[0] = 0

		for j in range(0,len(cigar_do)):

			# soft or hard trim: read position moves forward, parent position stays same
			if cigar_do[j]=="S":
				self.read_pos[j+1] = self.read_pos[j] + cigar_count[j]
				self.pos[j+1] = self.pos[j]

			if cigar_do[j]=="H":
				self.read_pos[j+1] = self.read_pos[j] + cigar_count[j]
				self.pos[j+1] = self.pos[j]

			# match: read position and parent position both move forward
			if cigar_do[j]=="M":
				for l in range(0,cigar_count[j]):
					new_loglike1, new_loglike2 = likelihood.get_ll(self.pos[j], l, self.read_pos[j], read,"M")
					self.loglike1 += new_loglike1
					self.loglike2 += new_loglike2

				self.read_pos[j+1] = self.read_pos[j] + cigar_count[j]
				self.pos[j+1] = add_c_to_pos(self.pos[j], cigar_count[j], self.mapped_parent)

			# insertion: read position moves forward, parent position stays same
			if cigar_do[j]=="I":
				for l in range(0,int(cigar_count[j])):	
					new_loglike1, new_loglike2 = likelihood.get_ll(self.pos[j], l, self.read_pos[j], read, "I")
					self.loglike1 += new_loglike1
					self.loglike2 += new_loglike2

				self.read_pos[j+1] = self.read_pos[j] + cigar_count[j]
				self.pos[j+1] = self.pos[j]

			# deletion: read position stays same, parent position moves forward
			if cigar_do[j]=="D":

				for l in range(0,int(cigar_count[j])):

					new_loglike1, new_loglike2=likelihood.get_ll(self.pos[j], l, self.read_pos[j], read, "D")				
					self.loglike1 += new_loglike1
					self.loglike2 += new_loglike2		

				self.read_pos[j+1] = self.read_pos[j]
				self.pos[j+1] = add_c_to_pos(self.pos[j], cigar_count[j], self.mapped_parent)
				
		return(self.pos)



