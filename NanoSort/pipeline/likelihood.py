import config
import numpy as np
import pandas as pd
import pickle

countmat1 = pd.read_table(config.count1, delimiter='\t', header=0, index_col=0)
countmat1 = countmat1.div(countmat1.sum(axis=1), axis=0)
countmat2 = pd.read_table(config.count2, delimiter='\t', header=0, index_col=0)
countmat2 = countmat2.div(countmat2.sum(axis=1), axis=0)

with open(config.fasta1) as fasta:
	s1 = list(fasta.read().split('\n')[1])
	l1 = len(s1)

with open(config.fasta2) as fasta:
	s2 = list(fasta.read().split('\n')[1])
	l2 = len(s2)

with open(config.pkl1,'rb') as pkl:
	p1 = pickle.load(pkl)

with open(config.pkl2,'rb') as pkl:
	p2 = pickle.load(pkl)

lim = min([l1,l2]) - 1
def get_ll(pos, l, read_pos, read, match_type):		

	if pos >= lim:
		kmer1 = "AAA"
		kmer2 = "AAA"

	else:
		if match_type != "I":
			kmer1 = "".join(s1[p1[pos]+l-2:p1[pos]+l+1])
			kmer2 = "".join(s2[p2[pos]+l-2:p2[pos]+l+1])
		else:
			kmer1 = "".join(s1[p1[pos]-2]+'-'+s1[p1[pos]+1])
			kmer2 = "".join(s2[p2[pos]-2]+'-'+s1[p2[pos]+1])

	if match_type != "D":
		read_nt = read[read_pos+l]
	else:
		read_nt = '-'
	
	if len(kmer1) < 3:
		loglike1 = 0
		loglike2 = 0
	if len(kmer2) < 3:
		loglike1 = 0
		loglike2 = 0
	else:
		loglike1 = np.log(countmat1.loc[kmer1,read_nt])
		loglike2 = np.log(countmat2.loc[kmer2,read_nt])

	return(loglike1,loglike2)

