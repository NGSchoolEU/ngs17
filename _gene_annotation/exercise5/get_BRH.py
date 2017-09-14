#!/usr/bin/env python

#Script used to obtain the BRH of a pair of reciprocal blast results

import argparse

def load_sequences(contigFile,delimiter):
	seqs = {}
	name = ""
	s = []
	for line in open(contigFile):
		line = line.strip()
		if ">" in line:
			if name != "":
				seqs[name] = "".join(s)
			if delimiter == "":
				name = line.replace(">","")
			else:
				name = line.replace(">","").split(delimiter)[0]
			s = []
		else:
			s.append(line.upper())
	if len(s) != 0:
		seqs[name] = "".join(s)
	return seqs

def filter_blast(BlastFile,seqs,threshold,evalue_thr=1e-05,overlap_thr=0.2):
	hits = {}
	stats = {}
	if evalue_thr > 0.1:
		print("WARNING: High e-value threshold!")
	for line in open(BlastFile):
		line = line.strip()
		dades = line.split("\t")
		if "#" not in line:
			try:
				overlap = (float(dades[7]) - float(dades[6])) / float(len(seqs[dades[0]]))
			except:
				overlap = 0.0
			if overlap > overlap_thr and float(dades[10]) < evalue_thr:
				if dades[0] not in hits:
					hits[dades[0]] = []
					stats[dades[0]] = {}
				if dades[1] not in hits[dades[0]] and len(hits[dades[0]]) < threshold:
					hits[dades[0]].append(dades[1])
					stats[dades[0]][dades[1]] = {}
					stats[dades[0]][dades[1]]["overlap"] = overlap
					stats[dades[0]][dades[1]]["eval"] = float(dades[10])
					stats[dades[0]][dades[1]]["ident"] = float(dades[2])
	return hits,stats

def get_BRH(hits1,hits2):
	BBH = []
	for c1 in hits1:
		h1 = hits1[c1][0]
		if h1 in hits2:
			if hits2[h1][0] == c1:
				if c1 > h1:
					pair = c1 +"-"+ h1
				else:
					pair = h1 + "-"+c1
				if pair not in BBH:
					BBH.append(pair)
	return BBH

parser = argparse.ArgumentParser(description="Obtains pairs of BRH")
parser.add_argument("-s1",dest="s1",action="store",required=True,help="Sequences of species 1")
parser.add_argument("-s2",dest="s2",action="store",required=True,help="Sequences of species 2")
parser.add_argument("-h1",dest="h1",action="store",required=True,help="Blast results from species 1 to species 2")
parser.add_argument("-h2",dest="h2",action="store",required=True,help="Blast results from species 2 to species 1")
parser.add_argument("-o",dest="outfileName",action="store",required=True,help="File where the results will be printed")
args = parser.parse_args()

seqs1 = load_sequences(args.s1," ")
seqs2 = load_sequences(args.s2," ")
hits1,stats = filter_blast(args.h1,seqs1,1)
hits2,stats = filter_blast(args.h2,seqs2,1)
pairs = get_BRH(hits1,hits2)
outfile = open(args.outfileName,"w")
for code in pairs:
	print >>outfile,code.replace("-","\t")
outfile.close()
