# encoding: utf-8
# __author:"xiaoyaoyin"

# import the required packages
# import pandas as pd
import numpy as np
import random
import time
import argparse


# calculate the GC content of given sequence
def calcu_CG(seq):
	# seq = seq.upper()
	count_a = seq.count('A')
	count_t = seq.count('T')
	count_c = seq.count('C')
	count_g = seq.count('G')
	gc_content = (count_g + count_c) / (count_a + count_c + count_g + count_t)
	return gc_content


# get the reverse of a sequence
def reverse(seq):
	new_seq = list(reversed(seq))
	return "".join(new_seq)


# get the complement of a sequence
def complement(seq):
	new_seq = []
	for i in seq:
		if i=='A':
			new_seq.append('T')
		elif i=='T':
			new_seq.append('A')
		elif i=='C':
			new_seq.append('G')
		elif i=='G':
			new_seq.append('C')
		else:
			new_seq.append('N')
	return "".join(new_seq)


# calculate distance of two seqs
def calcu_distance(p, q):
	maxn=max(len(p),len(q))
	minn=min(len(p),len(q))
	l=0;
	if len(p)==maxn:
		for k in range(maxn+minn-1):
			if k<minn:
				s = 0
				p1 = p[0:k]
				q1 = q[minn-k:minn]
				for j in range(k):
					if p1[j]==q1[j]:
						s += 1
				l = max(l, s)
			elif (k>=minn) & (k<=maxn):
				s = 0
				p1 = p[k-minn:k]
				q1 = q[0:minn]
				for j in range(minn):
					if p1[j]==q1[j]:
						s += 1
				l = max(l, s)
			else:
				s = 0
				p1 = p[k-minn:maxn]
				q1 = q[0:maxn+minn-k]
				for j in range(maxn+minn-k):
					if p1[j]==q1[j]:
						s += 1
				l = max(l, s)
	else:
		for k in range(maxn+minn-1):
			if k<minn:
				s = 0
				q1 = q[0:k]
				p1 = p[minn-k:minn]
				for j in range(k):
					if p1[j]==q1[j]:
						s += 1
				l = max(l, s)
			elif (k>=minn) & (k<=maxn):
				s = 0
				q1 = q[k-minn:k]
				p1 = p[0:minn]
				for j in range(minn):
					if p1[j]==q1[j]:
						s += 1
				l = max(l, s)
			else:
				s = 0
				q1 = q[k-minn:maxn]
				p1 = p[0:maxn+minn-k]
				for j in range(maxn+minn-k):
					if p1[j]==q1[j]:
						s += 1
				l = max(l, s)
	return l


# find if subsequence of length 5 exists
def subsequence(N):
	n = len(N)
	x = reverse(N)
	y = complement(x)
	m = 0
	for i in range(n-5):
		for j in range(n-5):
			match_res = [(x[i+k]==y[j+k]) for k in range(5)]
			if sum(match_res)!=5:
				m += 1
	if m==(n-5)**2:
		p1 = 1
	else:
		p1 = 0
	return p1


# convert any decimal number to base x, x should be smaller than 16
def decimal_to_x(n, x):
	a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 'A', 'B', 'C', 'D', 'E', 'F']
	b = []
	while True:
		s = n // x  
		y = n % x
		b = b + [y]
		if s==0:
			break
		n = s
	b.reverse()
	dd = [str(i) for i in b]
	return "".join(dd)


# design sticky ends sequences given length and number
def Seq_Design(lengthsticky, numsticky, param1, param3, param4):
	N0 = []
	M = 4
	d = [i for i in range(M**lengthsticky)]
	NN = [decimal_to_x(i, M).zfill(lengthsticky) for i in d]
	for i in range(len(NN)):
		temp_list = []
		for j in NN[i]:
			if j=='0':
				temp_list.append('A')
			elif j=='1':
				temp_list.append('T')
			elif j=='2':
				temp_list.append('C')
			else:
				temp_list.append('G')
		NN[i] = "".join(temp_list)
	for i in range(len(NN)):
		N = NN[i]
		# if the max length of single word repeats is larger than param1, kick it
		N = Kick_RepeatOne(N,param1)
		# if the max length of consecutive sequence is larger than param4, kick it
		N = Kick_RepeatSeq(N,param4)
		# if the CG ratio isn't between 0.4 to 0.6, kick it
		if len(N)>0:
			GC_percent = calcu_CG(N)
			if (GC_percent<=0.4) | (GC_percent>=0.6):
				N = []
		# compare this new sequence with those we already have, if their distance ratio is closer than q, kick it
		if len(N)>0:
			k = 0
			for j in range(len(N0)):
				k += (calcu_distance(N, N0[j])>=param3*lengthsticky)
			if k>0:
				N = []
		# compare this new sequence with itself to prevent hairpin
		if len(N)>0:
			if subsequence(N)==0:
				N = []
		if len(N)>0:
			N0.append("".join(N))
	return N0



# KICKREPEAT gives a set of seq without repeat ones
def Kick_RepeatOne(N, k):
	temp = 0
	for i in range(len(N)-k):
		if calcu_distance(N[i:i+k], N[(i+1):len(N)])==k:
			temp += 1
	if temp>0:
		M = []
	else:
		M = N
	return M


# KICKREPEAT gives a set of seq without repeat ones.
def Kick_RepeatSeq(N, k):
	#N the former set, pp the later set, k the repeat num,at least 2
	pp = []
	q1 = [0 for i in range(k-1)]
	for i in range(len(N)):
		temp = 0
		for j in range(len(N[i])-k+1):
			for m in range(k-1):
				q1[m] = (N[i][j]==N[i][j+m])
			if min(q1)==1:
				temp += 1
		if temp==0:
			pp.append(N[i])
	return pp


# main
if __name__ == '__main__': 
	start_time = time.time()
	parser = argparse.ArgumentParser(description='manual to this script')
	parser.add_argument("--lengthsticky", type=int, default=6)
	parser.add_argument("--numsticky", type=int, default=1)
	parser.add_argument("--single_word_repeats", type=int, default=3)
	# parse.add_argument("--max_overlap", type=int, default=5)
	parser.add_argument("--max_ratio", type=float, default=0.5)
	parser.add_argument("--max_repeats", type=int, default=5)
	args = parser.parse_args()
	lengthsticky = args.lengthsticky
	numsticky = args.numsticky
	param1 = args.single_word_repeats # the max length of single word repeats
	# param2 = args.max_overlap # the total length of overlap between two frame sequences like a,b,c,d,e
	param3 = args.max_ratio # the max ratio of overlap between two sticky sequences
	param4 = args.max_repeats # the max length of repeats
	S0 = Seq_Design(lengthsticky, numsticky, param1, param3, param4)
	R1 = [reverse(i) for i in S0]
	S1 = [complement(i) for i in R1]
	print('sticky ends designed')
	for i in range(len(S0)):
		print(S0[i])
	end_time = time.time()
	print('time used ',end_time-start_time)
