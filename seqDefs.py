#! /usr/bin/env python

import os, re, sys, math, time, stat, time, random, subprocess, decimal
from operator import add,mul
from string import *

from sio import *
import utilities as ut, xpermutations as perm


"""
This file contains methods pertaining to genomic sequences, io, manipulation and analysis.
"""
all_jobs = []

# ascii 33 to 126 - full range
fastqchars = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '\{', '\|', '\}', '~']

def getFastqchars():
	return fastqchars

# converting between phred and ascii and p values
A2Q = dict(zip(fastqchars, range(len(fastqchars))))
Q2P = lambda e: (e == 0 and 1-1e-16) or 10**(-e/10.) # avoid 1
A2P = dict( [ascii, 1-Q2P(A2Q[ascii])] for ascii in A2Q.keys() )

#
def hasnucs(astr, up=False):
	hasit = False
	if up: astr = astr.upper()
	for x in ['A','T','C','G']:
		if x in astr: hasit = True
	return hasit


# 4folds = ['G','V','P','T','A']
# 6folds = ['L','S','R']

def loadGeneticCode():
	"""Just does a dictionary lookup converting a nucleotide triplet to its amino acid sequence."""
	
	Z = {}
	Z['ATG'] = 'M'
	Z['TGG'] = 'W'
	
	Z['TGT'] = Z['TGC'] = 'C'
	Z['GAT'] = Z['GAC'] = 'D'
	Z['GAA'] = Z['GAG'] = 'E'
	Z['TTT'] = Z['TTC'] = 'F'
	Z['CAT'] = Z['CAC'] = 'H'
	Z['AAA'] = Z['AAG'] = 'K'
	Z['AAT'] = Z['AAC'] = 'N'
	Z['CAA'] = Z['CAG'] = 'Q'
	Z['TAT'] = Z['TAC'] = 'Y'
	
	Z['ATT'] = Z['ATC'] = Z['ATA'] = 'I'
	Z['TAA'] = Z['TAG'] = Z['TGA'] = '*'
	
	Z['GGT'] = Z['GGC'] = Z['GGA'] = Z['GGG'] = 'G'
	Z['GTT'] = Z['GTC'] = Z['GTA'] = Z['GTG'] = 'V'
	Z['CCT'] = Z['CCC'] = Z['CCA'] = Z['CCG'] = 'P'
	Z['ACT'] = Z['ACC'] = Z['ACA'] = Z['ACG'] = 'T'
	Z['GCT'] = Z['GCC'] = Z['GCA'] = Z['GCG'] = 'A'
	
	Z['TTA'] = Z['TTG'] = Z['CTT'] = Z['CTC'] = Z['CTA'] = Z['CTG'] = 'L'
	Z['TCT'] = Z['TCC'] = Z['TCA'] = Z['TCG'] = Z['AGT'] = Z['AGC'] = 'S'
	
	Z['CGT'] = Z['CGC'] = Z['CGA'] = Z['CGG'] = Z['AGA'] = Z['AGG'] = 'R'
	
	# also load reverse: from aa to codon
	Zi = {}
	Zi['M'] = ['ATG']
	Zi['W'] = ['TGG']
	
	Zi['C'] = ['TGT', 'TGC']
	Zi['D'] = ['GAT', 'GAC']
	Zi['E'] = ['GAA', 'GAG']
	Zi['F'] = ['TTT', 'TTC']
	Zi['H'] = ['CAT', 'CAC']
	Zi['K'] = ['AAA', 'AAG']
	Zi['N'] = ['AAT', 'AAC']
	Zi['Q'] = ['CAA', 'CAG']
	Zi['Y'] = ['TAT', 'TAC']
	
	Zi['I'] = ['ATT', 'ATC', 'ATA']
	Zi['*'] = ['TAA', 'TAG', 'TGA']
	
	Zi['G'] = ['GGT', 'GGC', 'GGA', 'GGG']
	Zi['V'] = ['GTT', 'GTC', 'GTA', 'GTG']
	Zi['P'] = ['CCT', 'CCC', 'CCA', 'CCG']
	Zi['T'] = ['ACT', 'ACC', 'ACA', 'ACG']
	Zi['A'] = ['GCT', 'GCC', 'GCA', 'GCG']
	
	Zi['L'] = ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG']
	Zi['S'] = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC']
	Zi['R'] = ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']
	
	return Z,Zi


Z,Zi = loadGeneticCode()
translateNT2AA = tr = nt2aa = lambda triple: 'N' in triple and 'X' or Z[triple]
aa2nt = rtr = lambda aa: aa in Zi and Zi[aa] or []

def translate(astr):
	if len(astr)%3 != 0: return ''
	
	aa = ''
	for i in range(0,len(astr),3):
		aa += tr(astr[i:i+3])
	return aa



def oligoHits(seq, oligo):
	"""Returns a list of indices where seq has a perfect 
	match to oligo. Assumes a 5'->3' sequence, so that 
	the index returned is the first and most 5' base in 
	the oligo. """
	dhits = dict()
	for i in range(len(seq)):
		index = seq.find(oligo, i)
		if index > -1:
			dhits[index] = 1
	# convert to list
	return dhits.keys()

def getSequence(db, id, start, stop):
	return 0

def iswellformed(seq='', strand='Plus'):
	"""Check whether an ORF looks OK based on multiple of 3 length and ATG...STOP"""
	wellformed = False
	if not len(seq): return False
	test = seq[:].upper()
	stdlow = strand.lower()
	if stdlow == 'plus' or stdlow == '+':
		# if test[0:3] == 'ATG':
		if test[0] == 'A':
			temp = test[len(seq)-3:len(seq)]
			if temp == 'TAG' or temp == 'TGA' or temp == 'TAA':
				wellformed = True
	else:
		# if test[len(seq)-3:len(seq)] == 'CAT':
		if test[len(seq)-1] == 'T':
			temp = test[0:3]
			if temp == 'CTA' or temp == 'TCA' or temp == 'TTA':
				wellformed = True
	return wellformed


def reverseComplementOld(seq='', uppercase=1):
	"""Translates a nucleotide sequence (including degenerate symbols) into
	its reverse complement. Always returns uppercase.
	"""
	seqlow = seq.upper()
	rcseq = ""
	for i in range(len(seq)-1, -1, -1):
		char = seqlow[i]
		# primary nucleotides
		if char == 'A': rcseq += 'T'
		elif char == 'T': rcseq += 'A'
		elif char == 'C': rcseq += 'G'
		elif char == 'G': rcseq += 'C'
		# now the degenerate ones
		elif char == '-': rcseq += '-'
		elif char == 'M': rcseq += 'K'	# m=(a/c) => k=(t/g)
		elif char == 'K': rcseq += 'M'
		elif char == 'R': rcseq += 'Y'	# r=(a/g) => y=(t/c)
		elif char == 'Y': rcseq += 'R'	
		elif char == 'S': rcseq += 'S'	# s=(c/g)
		elif char == 'V': rcseq += 'B'	# v=(a/g/c) => b=(t/g/c)
		elif char == 'B': rcseq += 'V'
		elif char == 'W': rcseq += 'W'	# w=(a/t)
		elif char == 'H': rcseq += 'D'	# h=(a/t/c) => d=(a/t/g)
		elif char == 'D': rcseq += 'H'
		else: rcseq += 'N'				# n=(a/t/g/c)
	if uppercase: return rcseq
	else: return rcseq.lower()

def reverseComplement(seq='', uppercase=0):
	"""Translates a nucleotide sequence (including degenerate symbols) into
	its reverse complement. Always returns uppercase.
	"""
	# seqlow = seq.upper()
	# rcseq = ""
	# for i in range(len(seq)-1, -1, -1):
		# c = seqlow[i]
	
	# m=(a/c) => k=(t/g)
	# r=(a/g) => y=(t/c)
	# s=(c/g)
	# v=(a/g/c) => b=(t/g/c)
	# w=(a/t)
	# h=(a/t/c) => d=(a/t/g)
	# n=(a/t/g/c)
	srev = lambda c:\
		(c=='A'and'T') or (c=='T'and'A') or (c=='C'and'G') or (c=='G'and'C') or (c=='U'and'A') or\
		(c=='a'and't') or (c=='t'and'a') or (c=='c'and'g') or (c=='g'and'c') or (c=='u'and'a') or\
		(c=='-'and'-') or (c=='M'and'K') or (c=='K'and'M') or (c=='R'and'Y') or (c=='Y'and'R') or\
		(c=='-'and'-') or (c=='m'and'k') or (c=='k'and'm') or (c=='r'and'y') or (c=='y'and'r') or\
		(c=='S'and'S') or (c=='V'and'B') or (c=='B'and'V') or (c=='W'and'W') or \
		(c=='s'and's') or (c=='v'and'b') or (c=='b'and'v') or (c=='w'and'w') or \
		(c=='H'and'D') or (c=='D'and'H') or 'N' or \
		(c=='h'and'd') or (c=='d'and'h') or 'n'
	
	rcseq = ''
	if len(seq) > 100e3:
		hlf = int(math.floor(len(seq)/3.0))
		# print 'mofo sequence', hlf, 2*hlf, len(seq)
		s1 = seq[0:hlf]; s2 = seq[hlf:2*hlf]; s3 = seq[2*hlf:len(seq)]
		rcs1 = map(srev,s1); rcs1.reverse()
		rcs2 = map(srev,s2); rcs2.reverse()
		rcs3 = map(srev,s3); rcs3.reverse()
		rcseq = ''.join(rcs3)+''.join(rcs2)+''.join(rcs1)
	else:
		rcseq = map(srev, seq)
		rcseq.reverse()
		rcseq = ''.join(rcseq)
	return (uppercase and rcseq) or rcseq

def rc(seq='', uppercase=1):
	return reverseComplement(seq,uppercase)

def reverse(string=''):
	"""reverses a string"""
	rev = ""
	for i in range(len(string)-1, -1, -1):
		rev += string[i]
	return rev


# need just a complement?


# NGS routines

def samToBed(line, pipe=sys.stdout):
	if line[0]=='@': return
	row = line[:-1].split('\t')
	if row[2] == '*': return
	strand = int(row[1]) & 16 and '-' or '+'
	# print 'test', row, row[3]
	# old
	# foobar = [row[0], '1', '%s%s'%(strand,row[2]), str(int(row[3])-1), str(int(row[3])+len(row[9])-1)]
	# new
	foobar = [row[2], str(int(row[3])-1), str(int(row[3])+len(row[9])-1),  row[0], '1', strand]
	# foobar = [row[0], '1', strand, row[2], str(int(row[3])-1), str(int(row[3])+len(row[9])-1)]
	# print 'foobar', foobar
	print >> pipe, '\t'.join(foobar)
	return




def nts2iupac(nts):
	X = {'A':'A', 'T':'T', 'C':'C', 'G':'G', 'N':'N', 'U':'U', '-':'-'}
	
	X['GT'] = 'K'
	X['AC'] = 'M'
	X['CGT'] = 'B'
	X['ACG'] = 'V'
	X['CG'] = 'S'
	X['AT'] = 'W'
	X['AGT'] = 'D'
	X['CT'] = 'Y'
	X['AG'] = 'R'
	X['ACT'] = 'H'
	
	# identities
	X['K'] = 'K'
	X['M'] = 'M'
	X['B'] = 'B'
	X['V'] = 'V'
	X['S'] = 'S'
	X['W'] = 'W'
	X['D'] = 'D'
	X['Y'] = 'Y'
	X['R'] = 'R'
	X['H'] = 'H'
	
	X['AA'] = 'A'
	X['TT'] = 'T'
	X['CC'] = 'C'
	X['GG'] = 'G'
	
	nts = ''.join(sorted(map(lambda x: x.upper(), nts)))
	return X[nts]


def iupac2nts(a,double=False):
	X = {'A':'A', 'T':'T', 'C':'C', 'G':'G', 'N':'N', 'U':'U', '-':'-'}
	
	X['K'] = 'GT'
	X['M'] = 'AC'
	X['B'] = 'CGT'
	X['V'] = 'ACG'
	X['S'] = 'CG'
	X['W'] = 'AT'
	X['D'] = 'AGT'
	X['Y'] = 'CT'
	X['R'] = 'AG'
	X['H'] = 'ACT'
	
	if double:
		X['A'] = 'AA'
		X['T'] = 'TT'
		X['C'] = 'CC'
		X['G'] = 'GG'
		X['N'] = 'NN'
		X['U'] = 'UU'
	
	return X[a.upper()]

def isIUPAC(astr):
	iup = ut.makeDict(['K','M','B','V','S','W','D','Y','R','H'])
	if len(filter(lambda x: x in iup, astr)) > 0: return True
	return False
def countIUPAC(astr):
	iup = ut.makeDict(['K','M','B','V','S','W','D','Y','R','H'])
	return len(filter(lambda x: x in iup, astr))

def nonnestedhamdist(s1,s2):
	assert len(s1) == len(s2)
	iup = ut.makeDict(['K','M','B','V','S','W','D','Y','R','H'])
	
	def nestednt(a,b):
		if len(ut.setIntersection(list(iupac2nts(a)), list(iupac2nts(b)))): return True
		return False
		
	return sum(not nestednt(a,b) for a,b in zip(s1,s2))


def readFasta(faname, usealtid=False, uniquenames=True, split='^>', idPat='', addins=[], include=[], TOUPPER=False, verbose=False,newline='\n', nentries=0):
	"""Parses a multiple fasta file into a dictionary object keyed by fasta identifier. 
	Option to key the dictionary by an alternative id found in the fasta header: altid=True|False.
	If directory of several fasta files is provided as faname, loads them into single dictionary.
	"""
	
	VERBOSE = verbose
	
	inkdct = dict( map(lambda x: (x,None), include) )
	
	fadict = {} # return dictionary of fasta entries
	faj = '' # load all information into single data string, faj
	if not os.access(faname, os.F_OK): print >> sys.stderr, 'readFasta(): File not found: %s'%(faname)
	files = getFiles(faname)
	
	for f in files:
		if VERBOSE: print '- Loading %s'%(f)
		fh = open(f)
		fa = fh.readlines()
		fh.close()
		# remove comment lines
		faj += ''.join(filter(lambda x: not re.match(' *#.*', x) and True, fa))
	if VERBOSE: print '- Parsing into blocks...'
	
	# parse by entry (>)
	getentries = re.compile(split, re.MULTILINE)
	faentries = re.split(getentries, faj)
	if VERBOSE: print '- Organizing blocks...'
	del faj
	
	# for some reason the first entry is the null string
	faentries.pop(0)
	
	repeatkeys = 0
	
	# parse individual entries
	while len(faentries):
		if nentries > 0 and len(fadict) >= nentries: break
	# for entry in faentries:
		entry = faentries.pop(0)
		# first  has format >name other-header-info
		# print entry
		(header, seq) = entry.split("%s"%(newline), 1)
		del entry
		
		# trim off the '>' character
		#header = header[0:len(header)]
		theid = ""
		altid = ""
		info = ""
		# further split up the info - for use of alternative identifier
		pattern = '^(\S+) (.*)$'
		if re.match(pattern, header):
			(crap, theid, info, crap) = re.split(pattern, header)
			info = info.strip(' ').strip('\t')
			if len(info) > 0:
				#print "\""+info+"\""
				(crap, altid, info, crap) = re.split('(\S+)(.*)', info)
				info = info.strip(' ').strip('\t')
		else:
			theid = header
			altid = ""
		
		if len(include) and theid not in inkdct: continue
		
		# remove newlines so the sequence data is contiguous
		seq = re.sub("\n", "", seq)
		# trim any carriage returns
		seq = re.sub("\r", "", seq)
		
		# remove terminal whitespace
		seq = seq.strip(' ')
		#seq = seq.replace('\n','')
		#print 'TEST!', seq
		
		if idPat != '':
			crap, theid, crap = re.split(idPat, theid)
			#print 'theid', theid
		
		for pre,suf,pat in addins:
			if re.match(pat, theid):
				# print 'theid', theid, 'pat', pat
				theid = pre+theid+suf
				# print 'new', theid
		
		if TOUPPER:
			seq = seq.upper()
		
		# print 'seq', seq[:200]
		
		# build the entry
		# add a little check to see if there are repeat entries
		if fadict.has_key(theid):
			if repeatkeys == 0 or verbose:
				sys.stderr.write('"%s": This key has already been used to store another fasta entry!'%(theid))
			repeatkeys += 1
			
			if uniquenames:
				theid += '.%s'%(repeatkeys)
				if repeatkeys == 1 or verbose:
					sys.stderr.write('Forcing unique id "%s"\n'%(theid))
			elif repeatkeys == 1 or verbose:
				sys.stderr.write('Skipping...\n')
			
		else:
			# build the entry for this id
			# allow choice of primary dict key
			if theid and altid and usealtid:
				fadict[altid] = {'id':theid, 'altid':altid, 'info': info, 'seq': seq}
			else:
				fadict[theid] = {'id':theid, 'altid':altid, 'info': info, 'seq': seq}
		
	if repeatkeys > 0:
		print >> sys.stderr, "A total of %s entries had shared keys"%(repeatkeys)
	
	return fadict

def fastaFromClustal(filename, verbose=False):
	"""Takes in a clustal multiple alignment file and returns a seqdict."""
	res = parseClustalOutput(filename, VERBOSE=verbose)
	return res['seq']


def readClustal(filename, verbose=False, keepgaps=True):
	return readFromClustal(filename, verbose, keepgaps)

def readFromClustal(filename, verbose=False, keepgaps=True):
	"""Takes in a clustal multiple alignment file and returns parsed alignment info."""
	return parseClustalOutput(filename, VERBOSE=verbose, keepgaps=keepgaps, rmTree=False)



def fasta2Phylip(dct={}, file='phylipi.txt', infile='', verbose=False):
	"""Convert a fasta file into Phylip interleaved format."""
	# print infile
	d = dct
	
	if infile: d = readFasta(infile)
	
	names = d.keys()
	# fix length
	maxlen = 10
	shorts = ['' for i in range(len(names))]
	for i in range(len(names)): 
		shorts[i] = ut.fixspacing(names[i], maxlen)+' '
	blank = ut.fixspacing('',maxlen+1)
	taxa = len(names)
	# print d
	# print 'test', d[d.keys()[0]]
	alilen = len(d[d.keys()[0]]['seq'])
	seqs = ['' for i in range(taxa)]
	
	if verbose: print shorts
	if verbose: print taxa, alilen
	if verbose: print map(len, seqs)
	enc = 'ascii'
	for i in range(taxa): seqs[i] = d[names[i]]['seq'].encode(enc)
	
	# print out an interleaved phylip formatted file
	pos = 0
	fh = open(file, 'w')
	# alilen = 715
	fh.write(str(taxa)+' '+str(alilen)+'\n'.encode(enc))
	
	for i in range(taxa):
		fh.write(shorts[i]+seqs[i][pos:pos+maxlen]+' '+seqs[i][pos+maxlen:pos+2*maxlen]+' '+seqs[i][pos+2*maxlen:pos+3*maxlen]+' '+seqs[i][pos+3*maxlen:pos+4*maxlen]+'\n'.encode(enc))
	pos = 2*maxlen
	fh.write('\n'.encode(enc))
	
	while pos < alilen:
		for i in range(taxa):
			fh.write(blank+seqs[i][pos:pos+maxlen]+' '+seqs[i][pos+maxlen:pos+2*maxlen]+' '+seqs[i][pos+2*maxlen:pos+3*maxlen]+' '+seqs[i][pos+3*maxlen:pos+4*maxlen]+'\n'.encode(enc))
		pos += 4*maxlen
		fh.write('\n'.encode(enc))
	
	fh.close()
	
	return 1


fasta2phylip = fasta2Phylip
printPhylip = fasta2Phylip

def isDiploid(fadct):
	"""Does this fasta file contain any diploid sequence?"""
	
	for k in fadct.keys():
		if isIUPAC(fadct[k]['seq']): return True
	return False


def readPhylip(filename, format='interleaved'):
	
	taxa = {}
	fh = open(filename)
	header = fh.readline().strip().split(' ')
	ntaxa = int(header[0])
	nbases = int(header[1])
	# read the first ntaxa lines for taxa names and first sequence
	eye2id = {}
	
	for i in range(ntaxa):
		line = fh.readline()
		row = line[:-1].split(' ')
		newid = row[0].strip()
		# print 'taxon', newid
		seq = ''.join(row[1:])
		taxa[newid] = {'seq':seq, 'info':'', 'id':newid}
		eye2id[i] = newid
	# now get remaining sequence
	crap = fh.readline()
	while crap:
		for i in range(ntaxa):
			line = fh.readline()
			row = line[:-1].strip().split(' ')
			seq = ''.join(row)
			taxa[eye2id[i]]['seq'] += seq
		crap = fh.readline()
	
	fh.close()
	
	return taxa


def readFastqReads(faname):
	"""General fastq formatted-file parser. Loads into a dictionary structure, keyed by ID"""
	
	D = {}
	fh = open(faname)
	lines = fh.readlines()
	fh.close()
	
	try:
		i = 0
		while i < len(lines):
			# entries have a 4-line structure:
			L1 = lines[i].strip()       # L1: @ID
			seq = lines[i+1].strip()    # L2: sequence
			L2 = lines[i+2].strip()     # L3: +ID
			Q = lines[i+3].strip()      # L4: quality sequence
		
			key = L1[1:len(L1)]
			key2 = L2[1:len(L2)]
			if key != key2:
				raise ValueError('Sequence and quality keys do not match:\n%s\n%s'%(key,key2))
			
			Qseq = map(lambda qj: A2Q[qj], Q)
			Pread = ut.prod(map(lambda qj: A2P[qj], Q)) # probability of no read error
			Qread = -10*math.log(1-Pread,10)
			D[key] = {'seq':seq, 'Q':Q, 'Qseq':Qseq, 'Pread':Pread, 'Qread':Qread, 'info':str(ut.r5(Qread)) }
		
			i += 4
	except IndexError: pass
	
	return D


def mergeFasta(target=[], reference=[], replchar='N', returnSNP=False):
	"""For each fasta entry in target and reference, update portions of target sequence which have "NNN" with corresponding reference sequence. Of course this assumes positional homology between target and reference."""
	
	print 'target keys:', target.keys()
	print 'reference keys:', reference.keys()
	
	SNP = {}
	
	new = {}
	for k in target.keys():
		if k in reference:
			print 'Match', k
			# find positions with replchar
			seq = target[k]['seq'].upper()
			rseq = reference[k]['seq'].upper()
			idx = seq.find(replchar)
			while idx != -1:
				# goofy string replace
				seq = seq[0:idx]+rseq[idx]+seq[idx+1:len(seq)]
				SNP[(k,idx,rseq[idx])] = 1
				# seq[idx] = rseq[idx]
				idx = seq.find(replchar)
			new[k] = target[k]
			new[k]['seq'] = seq
		else: # just keep old
			print 'No match', k
			new[k] = target[k]
	if returnSNP: return new, SNP
	return new


def diffFasta(target=[], reference=[], returnSNP=False):
	"""per-nucleotide diff between matching keys of 2 fasta dictionaries."""
	
	print 'target keys:', target.keys()
	print 'reference keys:', reference.keys()
	
	SNP = {}
	miscount = 0
	ncount = 0
	seqN = 0; rseqN = 0
	for k in target.keys():
		if k in reference:
			SNP[k] = []
			
			# find positions with replchar
			seq = target[k]['seq'].upper()
			rseq = reference[k]['seq'].upper()
			print 'Match', k, 'lseq/lrseq', len(seq), len(rseq)
			
			for i in range(len(seq)):
				try:
					if seq[i] == 'N': seqN += 1
					if rseq[i] == 'N': rseqN += 1
				
					if seq[i] != rseq[i] and (seq[i] != 'N' and rseq[i] != 'N'):
						miscount += 1
						# SNP[k] += [i,rseq[i],seq[i]]
						SNP[(k,i,seq[i],rseq[i])] = 1
					elif seq[i] != rseq[i] and seq[i] == 'N' or rseq[i] == 'N':
						ncount += 1
				except IndexError: pass
		else: # just keep old
			print 'No match', k
	if returnSNP: return {'miscount':miscount, 'ncount':ncount, 'arg1N':seqN, 'arg2N':rseqN}, SNP
	return miscount, ncount



def printFasta(fadict, file='', pipe=sys.stdout, key='seq', headerInfo=True, usealtid=False, maxwidth=None, addBlank=True, ioMethod='w', keyorder=None):
	"""Print a dictionary of fasta entries to a file in mfasta format"""
	
	keys = []
	if keyorder:
		keys = keyorder
	else:
		keys = sorted(fadict.keys())
	
	if file: pipe = open(file,ioMethod)
	
	if not maxwidth:
		for theid in keys:
			if 'altid' in fadict[theid] and not usealtid:
				header = '>'+str(theid)
				if headerInfo: header += ' '+str(fadict[theid]['altid'])+' '+str(fadict[theid]['info'])
			elif 'id' in fadict[theid]:
				header = '>'+str(theid)
				if headerInfo: header += ' '+str(fadict[theid]['id'])+' '+str(fadict[theid]['info'])
			else:
				header = '>'+str(theid)
				if headerInfo: header += ' '+str(fadict[theid]['info'])
			print >> pipe, header
			print >> pipe, fadict[theid][key]
			if addBlank: print >> pipe
	else:
		for theid in keys:
			if 'altid' in fadict[theid] and not usealtid:
				header = '>'+str(theid)
				if headerInfo: header += ' '+str(fadict[theid]['altid'])+' '+str(fadict[theid]['info'])
			elif 'id' in fadict[theid]:
				header = '>'+str(theid)
				if headerInfo: header += ' '+str(fadict[theid]['id'])+' '+str(fadict[theid]['info'])
			else:
				header = '>'+str(theid)
				if headerInfo: header += ' '+str(fadict[theid]['info'])
			print >> pipe, header
			
			# break up sequence into several lines each of length maxwidth characters
			tmp = fadict[theid][key]
			for zi in range(0,len(tmp),maxwidth):
				print >> pipe, tmp[zi:zi+maxwidth]
			if addBlank: print >> pipe
		
	if file: pipe.close()
	return 1

def printFastaEntry(fadict, id="noname", file='', pipe=sys.stdout, key='seq', ioMethod='w'):
	"""Print a dictionary to a file in fasta format"""
	if file: pipe = open(file,ioMethod)
	print >> pipe, ">"+id+" "+fadict['info']
	print >> pipe, fadict[key]
	if file: pipe.close()


def printSequential(fadict, file='', pipe=sys.stdout, key='seq', usealtid=False):
	"""Print a dictionary of fasta entries to a file in mfasta format"""
	keys = fadict.keys()
	keys.sort()
	
	if file: pipe = open(file,'w')
	
	ntaxa = len(keys)
	nseq = len(fadict[keys[0]]['seq'])
	print >> pipe, '%s %s'%(ntaxa, nseq)
	
	for theid in keys:
		print >> pipe, '%s  \n%s\n'%(fadict[theid]['id'], fadict[theid][key])
		# print >> pipe
	
	if file: pipe.close()


printPhylipSeq = printSequential

def list2fasta(seqs=[], ids=[], info=[]):
	dct = dict()
	uid = 0
	for i in range(len(seqs)):
		theid = str(uid)
		try: theid = ids[i]
		except IndexError: pass
		inf = ''
		try: inf = info[i]
		except IndexError: pass
		dct[theid] = {'id':theid, 'seq':seqs[i], 'info':inf}
		uid += 1
	return dct

def copyFasta(entry):
	"""Make a deep copy of a single fasta entry"""
	return {'id':entry['id'][:], \
	        'altid':entry['altid'][:], \
	        'info':entry['info'][:], \
	        'seq':entry['seq'][:]
	       }



def swapFastaIDs(fadict):
	"""Replaces altid with primary id and vice versa"""
	newfadict = {}
	for gob in fadict.keys():
		entry = copyFasta(fadict[gob])
		entry['id'] = entry['altid']
		entry['altid'] = gob
		newfadict[entry['id']] = entry
	return newfadict


def parseORFfromExtended(fadict):
	"""Given a fasta dictionary (e.g. returned by readFasta), trim both ends
	of each sequence that is lowercase and the CDS is uppercase
	"""
	newfadict = {}
	for gob in fadict.keys():
		entry = copyFasta(fadict[gob])
		seq = entry['seq']
		seq = ut.caseChangePos(seq,side='left',exclude=['N']) # handles 5' end
		seq = ut.caseChangePos(seq,side='right',exclude=['N']) # handles 3' end
		entry['seq'] = seq
		newfadict[gob] = entry
	return newfadict


def ORFsFromFile(filename, infoname, VERBOSE=False):
	""""""
	src = readFasta(filename)
	info = readDict(infoname)
	return ORFsFromGenome(src, info, VERBOSE=VERBOSE)

def ORFsFromGenome(fadict, info, VERBOSE=False):
	"""Given a raw genome (as dict) and corresponding exon track, 
	return the parsed ORF genome.
	
	WORK IN PROGRESS
	"""
	ORFeome = dict()
	notintrack = list()
	# get gene names from info table
	qkeys = info.keys()
	qkeys.sort()
	for qid in qkeys:
		# create object for qid
		fasta = {'id':qid,
		      'altid': 'NA',
		       'info': 'NA',
		        'seq': 'NA'}
		
		# the track may not list all IDs - this is problem
		if qid in info:
			# parse exon indices from track entry
			exonCount = int(info[qid]["exonCount"])
			exonStarts = info[qid]["exonStarts"].strip("\"").split(",")
			exonStarts.pop() # last entry is garbage
			exonStops = info[qid]["exonEnds"].strip("\"").split(",")
			exonStops.pop() # last entry is garbage
			exonStarts = map(int, exonStarts)
			exonStops = map(int, exonStops)
			
			# get which fasta/supercontig entry in target file
			
			# do a strand adjustment
			# relative/absolute indices
			# trim sequence by start and stop
			start = exonStarts[0]
			stop = exonStops[len(exonStops)-1]
			fasta['seq'] = fadict[qid]['seq'][start][stop]
		else:
			# do something smart here
			notintrack.append(qid)
			# just assume entry has one exon
			#continue
			exonCount = 1
		
		# store the entry
		ORFeome[qid] = fasta
	if VERBOSE:
		print >> sys.stderr, "ORFsFromGenome:", len(notintrack), "genes not in track."
	return ORFeome

def loadGenomes(files, verbose=False):
	# each file is mfasta - need to read them all in and group entries
	master = {}
	indepCount = 0
	comCount = 0
	singleCount = 0
	
	for f in files:
		genome = readFasta(f)
		if verbose: print 'read', len(genome.keys()), 'genome keys', f
		parts = f.split('/')
		species = parts[len(parts)-1].strip('.fasta')
		indepCount += len(genome.keys())
		for gid in genome.keys():
			try: master[gid][species] = genome[gid]
			except KeyError: master[gid] = {species:genome[gid]}
	comCount = len(master.keys())
	mk = master.keys()#; mk.sort()
	for gid in mk:
		num = len(master[gid].keys())
		if num == 1: singleCount += 1
		
	if verbose: print 'read total of', indepCount, 'sequences'
	if verbose: print singleCount, 'singletons remain'
	
	return master



def readUTRs(utrfile, idfile, default=ut.nastr, utr5='utr5', utr3='utr3', verbose=False):
	if verbose: print 'Loading UTR annotations...'
	# name2info = readDict(utrfile,key=True)
	(minfo, mrows, mcols) = readMatrix(utrfile, key=True) # read the table containing the UTR info
	# some rows have multiple entries, use the single best
	bestinfo = []
	bestrows = []
	i = 0
	while i < len(minfo):
		cid = mrows[i]
		qual = minfo[i][2]
		# print i, 'orig', cid, qual
		# are there more rows of same ID?
		j = i+1
		highestj = j
		while j < len(minfo):
			nid = mrows[j]
			nqual = minfo[j][2]
			highestj = j
			# print nid, '==', cid, '?'
			if nid == cid:
				# print '-- comp', j, nid, nqual
				# 'quality' column is transcript overlap - rank these
				if (nqual == '>=50%' and qual == '<50%') or \
				   (nqual == 'complete' and qual != 'complete'):
					cid = nid; qual = nqual
					i = j
					# print 'better', i
				j += 1
			else: j = len(minfo)
		# print 'keep', i
		bestinfo.append(minfo[i]); bestrows.append(cid)
		i = highestj
	
	# print 'len utr', len(bestinfo), len(bestrows)
	# sys.exit()
	name2info = dictFromMatrix(bestinfo, bestrows, mcols) # convert to dictionary
	name2id = readDict(idfile,key=True)
	
	genenames = name2info.keys(); genenames.sort()
	# print 'genenamnes', len(genenames)
	name2idfailures = []; utrpartialfailures = []; utrfullfailures = []
	workinggenes = {}; successcount = 0
	
	for cgname in genenames:
		#print gene, info[gene][heading]
		# the first thing we need is a conversion between the id of the gene 
		# in the fasta file and its id in the info table
		
		# split any compound gnames
		gnames = cgname.strip('\"').split(',')
		# print 'working', cgname
		# associate a sgd ID for this gene
		sgdid = ut.nastr
		gnotfound = True
		# first try the utr table itself
		featnames = name2info[cgname]['overlapFeatAll'].strip('\"').split(',')
		for anid in featnames:
			# name here is a sgdid... find it in name2id
			name = ut.deepLookup(name2id, anid, "id")
			# print 'anid', anid, 'short name', name
			# does this name match the name in gnames?
			idx = 0
			while idx < len(gnames):
				# print '-- comp', name, 'vs', gnames[idx]
				if gnames[idx] == name:
					sgdid = anid
					gnotfound = False
					idx = len(gnames)
				else: idx += 1
		# if not gnotfound: print 'match', sgdid, cgname, name
		idx = 0; name = ut.nastr
		
		while gnotfound and idx < len(gnames):
			name = gnames[idx]
			if name in name2id:
				sgdid = name2id[name]['id']
				if len(sgdid) > 0:
					gnotfound = False
				else:
					# in a few cases the registry table has an systematic id
					# but it falls under a different alias
					# perform a reverse lookup
					alias = ut.deepLookup(name2id, name, "alias")
					# print "found name:", alias
					gnames.append(alias)
			idx += 1
		if gnotfound: 
			name2idfailures.append(cgname)
			continue # move to next item
		
		# get the utr indices from the info table
		utr5len = name2info[cgname][utr5]
		utr3len = name2info[cgname][utr3]
		# print cgname, utr5, utr5len, utr3, utr3len
		
		# now filter the successes by which actually have listed UTR values
		if utr5len != "NA" and utr3len != "NA":
			workinggenes[sgdid] = {'utr5':int(utr5len), 'utr3':int(utr3len)}
			successcount += 1
		elif (utr5len != "NA" and utr3len == "NA"):
			workinggenes[sgdid] = {'utr5':int(utr5len), 'utr3':default}
			successcount += 1
		elif (utr5len == "NA" and utr3len != "NA"):
			workinggenes[sgdid] = {'utr5':default, 'utr3':int(utr3len)}
			successcount += 1
		else:
			utrfullfailures.append(cgname)
	
	if verbose:
		print "--- DATA SUMMARY ---"
		print "SUCCESSES, cleanly parsing", successcount, "genes"
		print "NAME2ID LOOKUP FAILURES:", len(name2idfailures)
		# print "FASTA DATA LOOKUP FAILURES:", len(datafailures)
		print "UTR PARTIAL FAILURES:", len(utrpartialfailures)
		print "UTR FULL FAILURES:", len(utrfullfailures)
		print
		print "WORKING SET OF GENES WITH UTR DATA:", len(workinggenes.keys())
	
	return workinggenes



################ MISC ################

def fastqToFasta(f, header, newfile):
	outdir = ut.getPath(newfile)
	fastafile1 = outdir+ut.getLabel(newfile)+'.fasta'
	fastafile2 = outdir+ut.getLabel(newfile)+'.2.fasta'
	fastafile3 = outdir+ut.getLabel(newfile)+'.fin.fasta'
	
	# how to convert fastq file to fasta for blast analysis
	call = 'grep -A 1 "@%s" "%s" > "%s"'%(header, f, fastafile1)
	print '  Running:', call
	os.system(call)
	call = 'sed -e "s/@%s/>%s/" "%s" > "%s"'%(header, header, fastafile1,fastafile2)
	print '  Running:', call
	os.system(call)
	call = 'sed -e "s/--//" "%s" > "%s"'%(fastafile2, fastafile3)
	print '  Running:', call
	os.system(call)
	call = 'mv "%s" "%s"'%(fastafile3, fastafile1)
	print '  Renaming:', call
	os.system(call)
	print '  Removing temporary files...'
	os.system('rm "%s"'%(fastafile2))
	# os.system('rm "%s"'%(fastafile3))
	



############ ALIGNMENTS ############
def subseqHelper(seq, lengths, a):
	# replace -1 by full length sequence
	s1 = seq[:]
	ls = lengths[:]
	
	if ls[0] == '-1': ls[0] = len(s1)
	if ls[1] == '-1': ls[1] = len(s1)
	# replace negative number by full length minus that number
	if ls[0] < 0: ls[0] = len(s1) + ls[0]
	if ls[1] < 0: ls[1] = len(s1) + ls[1]
	# hackish
	if a == '5primeInner':
		# get subseq from 500 to end of UTR
		ls[0] = max(len(s1)-ls[0], 0)
		ls[1] = ls[1]
	elif a == '5primeOuter':
		# get subseq from seqlen to 500
		ls[1] = max(len(s1)-ls[1], 0)
	elif a == '3primeInner':
		# get subseq from end of UTR to 500
		ls[1] = min(len(s1), ls[1])
	elif a == '3primeOuter':
		# get subseq from 500 to seqlen
		ls[0] = min(len(s1), ls[0])
	return s1[ls[0]:ls[1]]


def blast2Seqs(seq1="", seq2="", params='-F F -p blastn', savepath="", name='gene', DELETE=True, VERBOSE=False):
	"""Given two strings to be aligned, run the pairwise local alignment program
	BLAST from NCBI, parse the resulting output file, and return 
	a dictionary, where key 'qseq' returns the query alignment, 'tseq' the target
	alignment, and 'ngaps', 'lenwithgaps', and 'niden' alignment statistics."""
	
	#savepath = savepath.strip('/')
	if len(savepath) and savepath[len(savepath)-1] != '/': 
		#savepath = savepath[0:len(savepath)-2]
		savepath += '/'
	s1file = savepath+name+'_query.fasta'
	s2file = savepath+name+'_target.fasta'
	bl2file = savepath+name+'_bl2seq.txt'
	
	s1 = open(s1file, "w")
	q = {"info":"", "seq":seq1}
	printFastaEntry(q, name+'_query', pipe=s1)
	s1.close()
	
	s2 = open(s2file, "w")
	t = {"info":"", "seq":seq2}
	printFastaEntry(t, name+'_target', pipe=s2)
	s2.close()
	
	blcall = "bl2seq -i "+s1file+" -j "+s2file+" "
	blcall += params
	blcall += " > \'"+bl2file+"\'"
	if VERBOSE: print 'Call:', blcall
	os.system(blcall)
	
	# save the alignment sequences for each gene
	#[qseq, tseq] = parseBlastOutput(bl2file)
	blresults = parseBlastOutput(bl2file)
	if VERBOSE: print 'Results:', blresults
	
	if DELETE: 
		os.remove(s1file)
		os.remove(s2file)
		os.remove(bl2file)
	
	# return the results dictionary
	return blresults


def blastall(queryfile, dbfile, parms='-F F -p blastn -b 1 -m 8 -M BLOSUM80', expect=10.0, savepath='', verbose=True, reformatdb=False, redo=True, returnFilename=False, blastresultsfile='', ispeptide=False, nprocs=1):
	
	if blastresultsfile == '':
		blastresultsfile = savepath+'blastall_table.txt'
	
	if not os.access(savepath+'formatdb.log', os.F_OK) or not os.access(dbfile+'.nin', os.F_OK) or reformatdb:
		# print 'blast redoing', dbfile+'.nin'
		# format the target file as a blast database
		formatcall = 'formatdb -i \''+dbfile+'\' -p %s -l \''%(ispeptide and 'T' or 'F')
		formatcall += savepath+'formatdb.log\''
		if verbose: 
			print "Formatting blast db..."
			print formatcall
		os.system(formatcall)
		if verbose: print "formatdb.log created."
	
	if not os.access(blastresultsfile, os.F_OK) or redo:
		# run a blastall against the target database
		if verbose: print "Running blastall against target database..."
		blastcall = 'blastall -d \''+dbfile+'\' -i \''+queryfile
		blastcall += '\' '+parms+' -e %s'%(expect)
		blastcall += ' > \''+blastresultsfile+'\''
		os.system(blastcall)
		if verbose: print 'Call:', blastcall
	
	if not returnFilename:
		if verbose: print "Parsing Blast file..."
		# headings = ['qid', 'tid', 'pct', 'hitlen', 'mismatches', 'gapopens', 'qstart', 'qstop', 'tstart', 'tstop', 'evalue', 'bitscore']
		Qgenes = parseBlastFile(blastresultsfile)
		if verbose: print "finished."
		return Qgenes
	else:
		return blastresultsfile


def launchProcesses(jobs, nsubprocs=1, sleep=1, verbose=0):
	# global all_jobs, running_jobs, nlaunced, ncomplete, nkilled
	
	global all_jobs
	all_jobs += jobs
	
	# multiprocessing queue
	# ---------------------
	njobs = len(all_jobs) # initial job count
	nlaunched = 0
	running_jobs = []
	ncomplete = 0
	nkilled = 0
	
	while len(all_jobs):
		# as long as there are more jobs to be processed
		if nlaunched - ncomplete - nkilled < nsubprocs:
			if verbose: print >> sys.stderr, ' - Scheduling job', nlaunched+1, '(%s completed, %s killed | x%s processors)'%(ncomplete, nkilled, nsubprocs)
			# and we can launch another job
			# fetch the next job
			call = all_jobs.pop(0)
			# and execute it
			job = subprocess.Popen(call, shell=True)
			running_jobs += [job]
			nlaunched += 1
			time.sleep(sleep)
		else:
			# pipeit('waiting...')
			# check if existing job has been killed
			nkilled = 0
			marked = []
			for ji in xrange(len(running_jobs)): 
				if running_jobs[ji].poll() != None: 
					# print 'poll killed', rj[ji].poll(), len(running_jobs)
					ncomplete += 1
					marked += [running_jobs[ji]]
			for m in marked:
				running_jobs.pop(running_jobs.index(m))
					
			# if nkilled > 0: print '%s jobs killed'%(nkilled)
			# check to see if an existing job has finished
			# ncomplete = len(getFiles(completedir, include='likelihood'))
			time.sleep(sleep)
		# print 'tet', nlaunched
	for job in running_jobs: job.wait() # halt code until all jobs are finished.
	return 1


def batchGlobalAlignment(seq1, fadct, program='stretcher', nprocs=1, matrix='', params='', savepath='./', name='gene', rmAln=True, rmAlt=True, refresh=True, hardCut=10, indels=False, terminals='-', keepgaps=True, verbose=False):
	"""Run many alignment calls in parallel
	seq1 is a string
	fadct is a dict of additional sequences
	"""

	calls = []

	if savepath[len(savepath)-1] != '/': savepath += '/'
	jobs = []
	counter = 1
	outfiles = []
	for seq2lab in fadct.keys():
		myname = name+'vs.%s.%s'%(seq2lab,counter)
		counter += 1

		s1file = savepath+myname+'_query.fasta'
		s2file = savepath+myname+'_target.fasta'
		outfile = savepath+myname+'_stretcher.txt'
		
		seq2 = fadct[seq2lab]['seq']

		RUNWITHIT = True
		if hardCut != None and (len(seq1) > hardCut*len(seq2) or len(seq2) > hardCut*len(seq1)): 
			print >> sys.stderr, 'Error: sequence too long for stretcher'
			RUNWITHIT = False
		
		# if forced alignment or not forced, but can't access results, redo alignment
		if RUNWITHIT and (refresh or (not refresh and not os.access(outfile, os.F_OK))):
			s1 = open(s1file, "w")
			q = {"info":"", "seq":seq1}
			printFastaEntry(q, myname+'_query', pipe=s1)
			s1.close()
		
			s2 = open(s2file, "w")
			t = {"info":"", "seq":seq2}
			printFastaEntry(t, myname+'_target', pipe=s2)
			s2.close()
			alcall = ''
			if program == 'stretcher':
				alcall += 'stretcher -asequence "'+s1file+'" -bsequence "'+s2file+'" '
				if params:
					alcall += params
				else:
					alcall += ' -gapopen 16 -gapextend 4 '
			else:
				alcall += 'needle -asequence "'+s1file+'" -bsequence "'+s2file+'" '
				if params:
					alcall += params
				else:
					alcall += ' -gapopen 10.0 -gapextend 0.5 '
			if matrix:
				alcall += '-datafile '+matrix+' '
			alcall += ' -outfile "'+outfile+'" '
			alcall += "&> /dev/null"
			
			if verbose: print 'Call:', alcall
			jobs += [alcall]
			# os.system(alcall)
			outfiles += [(s1file,s2file,outfile,seq2lab)]

	launchProcesses(jobs,nprocs, verbose=verbose)

	# read results
	resset = {}
	for s1file,s2file,outfile,seq2lab in outfiles:
		results = parseSNAlignOutput(outfile, indels=indels, terminals=terminals)
		resset[seq2lab] = results
		if rmAlt and os.access(s1file,os.F_OK): os.remove(s1file)
		if rmAlt and os.access(s2file,os.F_OK): os.remove(s2file)
		if rmAln and os.access(outfile,os.F_OK): os.remove(outfile)
		
	# return the results dictionary
	return resset




def globalAlignment(seq1="", seq2="", program='stretcher', matrix='', params='', savepath='./', name='gene', rmAln=True, rmAlt=True, refresh=True, hardCut=10, indels=False, terminals='-', keepgaps=True, verbose=False):
	"""Given two strings to be aligned, run the pairwise global alignment program
	stretcher from the EMBOSS package, parse the resulting output file, and return 
	a dictionary, where key 'qseq' returns the query alignment, 'tseq' the target
	alignment, and 'ngaps', 'lenwithgaps', and 'niden' alignment statistics.
	
	program can be either stretcher or needle. Use stretcher for long sequences as it saves space, but needle saves time by a factor of 2.
	"""
	
	if program == 'clustal' or program=='clustalw':
		return multipleAlignment({'q'+name:{'seq':seq1, 'info':''}, 
		                            't'+name:{'seq':seq2, 'info':''}}, keepgaps=keepgaps,
		                            tempdir=savepath, rmAln=rmAln, rmAlt=rmAlt, 
		                            indels=indels, terminals=terminals)
	
	#savepath = savepath.strip('/')
	if savepath[len(savepath)-1] != '/': 
		#savepath = savepath[0:len(savepath)-2]
		savepath += '/'
	s1file = savepath+name+'_query.fasta'
	s2file = savepath+name+'_target.fasta'
	outfile = savepath+name+'_stretcher.txt'
	
	RUNWITHIT = True
	if hardCut != None and (len(seq1) > hardCut*len(seq2) or len(seq2) > hardCut*len(seq1)): 
		print >> sys.stderr, 'Error: sequence too long for stretcher'
		RUNWITHIT = False
	
	# if forced alignment or not forced, but can't access results, redo alignment
	if RUNWITHIT and (refresh or (not refresh and not os.access(outfile, os.F_OK))):
		s1 = open(s1file, "w")
		q = {"info":"", "seq":seq1}
		printFastaEntry(q, name+'_query', pipe=s1)
		s1.close()
	
		s2 = open(s2file, "w")
		t = {"info":"", "seq":seq2}
		printFastaEntry(t, name+'_target', pipe=s2)
		s2.close()
		alcall = ''
		if program == 'stretcher':
			alcall += 'stretcher -asequence "'+s1file+'" -bsequence "'+s2file+'" '
			if params:
				alcall += params
			else:
				alcall += ' -gapopen 16 -gapextend 4 '
		else:
			alcall += 'needle -asequence "'+s1file+'" -bsequence "'+s2file+'" '
			if params:
				alcall += params
			else:
				alcall += ' -gapopen 10.0 -gapextend 0.5 '
		if matrix:
			alcall += '-datafile '+matrix+' '
		alcall += ' -outfile "'+outfile+'" '
		alcall += "&> /dev/null"
		
		if verbose: print 'Call:', alcall
		os.system(alcall)
	
	results = parseSNAlignOutput(outfile, indels=indels, terminals=terminals)
	# print 'score', results['score']
	# print 'globalAli', 'removing?', rmAln, outfile, savepath
	if rmAlt and os.access(s1file,os.F_OK): os.remove(s1file)
	if rmAlt and os.access(s2file,os.F_OK): os.remove(s2file)
	if rmAln and os.access(outfile,os.F_OK): os.remove(outfile)
		
	# return the results dictionary
	return results

def efGlobalAlignment(seq1="", seq2="", program='clustalw', matrix='', params='', savepath='./', name='gene', rmAln=True, rmAlt=True, refresh=True, hardCut=10, indels=False, terminals='-'):
	"""Pseudo-end free global alignment"""
	results = {}
	symbol = '*'
	if program == 'stretcher':
		results = globalAlignment(seq1,seq2,program,matrix,params,savepath,
					name,rmAln,rmAlt,refresh,hardCut, 
					indels=indels, terminals=terminals)
		symbol = ':'
	else:
		results = multipleAlignment({'q'+name:{'seq':seq1, 'info':''}, 
		                               't'+name:{'seq':seq2, 'info':''}}, 
		              tempdir=savepath, rmAln=rmAln, rmAlt=rmAlt, 
		              indels=indels, terminals=terminals)
	try:
		# just want from first match to last, without outer gaps
		first = results['idenseq'].index(symbol)
		ids = map(lambda x:x, results['idenseq'])
		ids.reverse(); ids = ''.join(ids); last = len(ids) - ids.index(symbol)
		# print 'first and last', first, last
		if program == 'stretcher':
			results['qseq'] = results['qseq'][first:last]
			results['tseq'] = results['tseq'][first:last]
		else:
			results['qseq'] = results['seq']['q'+name]['seq'][first:last]
			results['tseq'] = results['seq']['t'+name]['seq'][first:last]
			del results['seq']
		
		results['idenseq'] = results['idenseq'][first:last]
		results['iden'] = len(ut.filterstr(' ')(map(lambda x:x,results['idenseq'])))
		results['lenwithgaps'] = len(results['idenseq'])
		results['niden'] = results['lenwithgaps'] - results['iden']
	except ValueError: 
		results['qseq'] = ''
		results['tseq'] = ''
		pass
	return results



def multipleAlignment(mfDct, fn='clustalfile', program='clustalw', rmAln=True, rmAlt=True, tempdir='', refresh=True, indels=False, terminals='-', keepgaps=False, readTree=False):
	
	if not tempdir: tempdir = os.getcwd()
	if tempdir[len(tempdir)-1] != '/': tempdir += '/'
	outfile = tempdir+fn
	results = {}
	nseqs = len(mfDct.keys())
	
	if refresh or (not refresh and not os.access(outfile+'.aln', os.F_OK)):
		createdir(tempdir)
		oh = open(outfile+'.fasta', 'w')
		printFasta(mfDct, pipe=oh)
		oh.close()
		
		results = {}
		if program == 'clustalw':
			call = 'clustalw -align -infile='+outfile+'.fasta'+' &> /dev/null'
			os.system(call)
			# print 'will trim gap', terminals
			results = parseClustalOutput(outfile+'.aln', indels=indels, terminals=terminals,keepgaps=keepgaps, readTree=readTree, rmTree=rmAlt)
			
		if rmAln: os.system('rm "'+outfile+'.aln"')
		if rmAlt: 
			os.system('rm "'+outfile+'.dnd"')
			os.system('rm "'+outfile+'.fasta"')
	else: 
		results = parseClustalOutput(outfile+'.aln', indels=indels, terminals=terminals, readTree=readTree, rmTree=rmAlt)
		
	return results


def multipleAlignmentWithLengths(mfDct, fn='clustalfile', program='clustalw', rmAln=True, rmAlt=True, tempdir='', refresh=True, indels=False, terminals='-', lengths={}, returnSeq=False, subSeqOnly=False, readTree=False):
	"""Similar to allPairsAlignment given optional argument lengths, this 
	 performs several multiple alignments on subsequences of the given input 
	 sequences, partitioned according to [start,end] pairs in the list lengths
	"""
	
	if not tempdir:
		tempdir = os.getcwd()
	if tempdir[len(tempdir)-1] != '/': tempdir += '/'
	outfile = tempdir+fn
	
	keys = mfDct.keys()
	
	idenseqs = []
	# save all alignments into subdirectory (fn must be specified for this to work)
	if not rmAln:
		tempdir += fn
		createdir(tempdir)
	
	if not len(lengths.keys()): lengths['full'] = [0,-1,'-']
	
	results = {}
	for a in lengths.keys():
		results[a] = {}
		results[a]['pi'] = []
		results[a]['score'] = []
		results[a]['iden'] = []
		results[a]['snp'] = []
		results[a]['ngaps'] = []
		results[a]['lenwithgaps'] = []
	
	# global align this pair for all sub-parts
	for a in lengths.keys():
		# print a, fn, '-----------------------------'
		# create subseq dict, making sure that lengths match sequence properly
		# python internally handles when n > len(s)
		subDct = {}
		for k in keys:
			entry = copyFasta(mfDct[k])
			entry['seq'] = subseqHelper(entry['seq'], lengths[a], a)
			subDct[k] = entry
		
		if not subSeqOnly:
			# determine whether to trim or count terminal gaps
			terminals = lengths[a][2]
			# print 'terminals?', a, terminals
			gresults = multipleAlignment(subDct, tempdir=tempdir, fn=a+'_'+fn, 
			  program=program, rmAln=rmAln, rmAlt=rmAlt, refresh=refresh, 
			  indels=indels, terminals=terminals, readTree=readTree)
		
			# figure out results for subsequence
			results[a]['pi'] = gresults['pi']
			results[a]['score'] = gresults['score']
			results[a]['iden'] = gresults['iden']
			results[a]['snp'] = gresults['niden']
			results[a]['ngaps'] = gresults['ngaps']
			results[a]['lenwithgaps'] = gresults['lenwithgaps']
			# off by default to conserve memory
			if returnSeq: results[a]['seq'] = gresults['seq']
		else:
			results[a]['seq'] = subDct
	return results


def subSeqs(mfDct,lengths,a):
	# if not len(lengths.keys()): lengths['full'] = [0,-1,'-']
	
	results = {}
	# for a in lengths.keys(): results[a] = {}
	
	# for k in mfDct.keys():
	entry = copyFasta(mfDct)
	# for a in lengths.keys():
		# create subseq dict, making sure that lengths match sequence properly
		# python internally handles when n > len(s)
	entry['seq'] = subseqHelper(entry['seq'], lengths, a)
	# results[a][k] = entry
	return entry


def allPairsAlignment(mfDct, fn='clustalfile', program='stretcher', rmAln=True, rmAlt=True, tempdir='', refresh=True, indels=False, terminals='-', lengths={}):
	"""Lengths is optional, but will return all pairs results for each key:(start,end) item
	"""
	
	if not tempdir:
		tempdir = os.getcwd()
	if tempdir[len(tempdir)-1] != '/': tempdir += '/'
	outfile = tempdir+fn
	
	# want to do global alignment for all unique pairs of sequences to obtain
	# an appropriate overall distance value
	keys = mfDct.keys()
	
	idenseqs = []
	# save all alignments into subdirectory (fn must be specified for this to work)
	if not rmAln:
		tempdir += fn
		createdir(tempdir)
	
	if not len(lengths.keys()): lengths['full'] = [0,-1,'-']
	nkeys = len(lengths.keys())
	
	results = {}
	for a in lengths.keys():
		results[a] = {}
		results[a]['pi'] = []
		results[a]['V(pi)'] = ut.nastr
		results[a]['score'] = []
		results[a]['iden'] = []
		results[a]['snp'] = []
		results[a]['V(snp)'] = ut.nastr
		results[a]['ngaps'] = []
		results[a]['lenwithgaps'] = []
	
	combos = ut.makeCombinations(len(keys), 2, start_idx=0)
	for c in combos:
		k1 = keys[c[0]]; k2 = keys[c[1]]
		s1 = mfDct[k1]['seq']; s2 = mfDct[k2]['seq']
		
		# global align this pair for all sub-parts
		for a in lengths.keys():
			pref = (nkeys > 1 and a+'_') or ''
			# get subseqs, making sure that lengths match sequence properly
			# python internally handles when n > len(s)
			subs1 = subseqHelper(s1, lengths[a], a)
			subs2 = subseqHelper(s2, lengths[a], a)
			# determine whether to trim or count terminal gaps
			terminals = lengths[a][2]
			# print 'terminals?', a, terminals
			gresults = globalAlignment(subs1, subs2, savepath=tempdir, name=pref+k1+'_'+k2,
			  program=program, rmAln=rmAln, rmAlt=rmAlt, refresh=refresh, 
			  indels=indels, terminals=terminals)
			
			aliLen = gresults['lenwithgaps']
			nIden = gresults['iden']
			nNiden = gresults['niden']
			# nQgaps = aliLen - len(''.join(sqseq.split('-')))
			# nTgaps = aliLen - len(''.join(stseq.split('-')))
			
			results[a]['lenwithgaps'].append(aliLen)
			results[a]['iden'].append(nIden)
			results[a]['snp'].append(nNiden)
			results[a]['ngaps'].append(gresults['ngaps'])
			results[a]['score'].append(gresults['score'])
			results[a]['pi'].append(gresults['pi'])
	
	for a in lengths.keys():
		# figure out results for subsequence
		pis = results[a]['pi']
		results[a]['pi'] = ut.avg(pis)
		results[a]['V(pi)'] = ut.variance(pis)
		results[a]['score'] = ut.avg(results[a]['score'])
		results[a]['iden'] = ut.avg(results[a]['iden'])
		snps = results[a]['snp']
		results[a]['snp'] = ut.avg(snps)
		results[a]['V(snp)'] = ut.variance(snps)
		results[a]['ngaps'] = ut.avg(results[a]['ngaps'])
		results[a]['lenwithgaps'] = ut.avg(results[a]['lenwithgaps'])
	
	return results


### PARSING ALIGNMENT OUTPUT FILES ###
def fixIndelHelper(seqs, trimTerminal='-'):
	# treat all indel runs as a single mismatch / SNP
	# for each position, need to check each character and
	# see whether it is a gap, and build a subseq of positions
	# with a gap
	lenwithgaps = 0
	try:
		lenwithgaps = len(seqs[0]) # all sequences assumed to have same length - already aligned
	except IndexError: pass
	nseqs = len(seqs)
	pointGaps = 0
	groupGaps = 0
	gapruns = [] # store start,end indices for each gap run
	start = 0
	end = 0
	
	while start < lenwithgaps:
		end = start
		# extend a gap (if exists)
		gapRun = True    # will enter into an attempted gap run
		gapInRun = False # have found a gap in this attempt
		while gapRun and end < lenwithgaps:
			# check for a gap in this column
			gapInCol = False
			for i in range(nseqs):
				if seqs[i][end] == '-': gapInCol = True
			if gapInCol: 
				pointGaps += 1
				end += 1
				gapInRun = True
			else: gapRun = False
		if gapInRun: groupGaps += 1
		# end += 1
		oldstart = start
		# if end > start: start += end-start
		if end > start+1: start += end-start-1
		gapruns.append([oldstart, start+1])
		start += 1
	
	# print 'finished at start', start
	
	# trim terminal gap runs if desired
	if trimTerminal != '-' and lenwithgaps > 0 and groupGaps > 0:
		# print 'TRIMMING TERMINALS', len(gapruns), lenwithgaps
		newgapruns = gapruns[:]
		# find a gap run with start = 0 or end = lenwithgaps
		hits = 0
		if trimTerminal.count('5') > 0:
			if gapruns[0][0] == 0: newgapruns.pop(0); hits += 1
		if trimTerminal.count('3') > 0:
			if gapruns[len(gapruns)-1][1] == lenwithgaps: 
				newgapruns.pop(len(newgapruns)-1); hits += 1
		groupGaps -= hits
		
		# if hits > 0: 
			# if newgapruns[0][0] != 0: print 'found', hits, 'terminal runs start'
			# if newgapruns[len(newgapruns)-1][1] != gapruns[len(gapruns)-1][1]: print 'found', hits, 'terminal runs end'
			# print 'wtf', newgapruns[len(newgapruns)-1][1], gapruns[len(gapruns)-1][1]
			# print 'orig', gapruns[0:2], gapruns[len(gapruns)-2:len(gapruns)]
			# print 'now', newgapruns[0:2], newgapruns[len(newgapruns)-2:len(newgapruns)]
		# else:
			# print 'no hits'
		
		gapruns = newgapruns
		# print 'final length', len(gapruns)
	
	return {'groupGaps':groupGaps, 'lenwithgaps':len(gapruns), 'pointGaps':pointGaps}


# parsing global pairwise and multiple alignments
# -----------------------------------------------
def parseSNAlignOutput(filename, VERBOSE=False, indels=False, terminals='-', missing='N'):
	"""Given the path to a stretcher alignment results file, this
	method will read the file, and return a dictionary, 
	
	NB these sequences will be the actual alignment results, not the input 
	sequence given to bl2seq."""
	# print 'got indels', indels
	# results dictionary
	results = {'qseq':'', 'tseq':'','idenseq':'', 'pi':ut.nastr, 'iden':ut.nastr,
	           'niden':ut.nastr, 'lenwithgaps':ut.nastr, 'ngaps':ut.nastr, 
	           'score':ut.nastr, 'file':filename}
	
	# get the alignment file data
	data = []
	try:
		file = open(filename)
		data = file.readlines()
		file.close()
	except IOError:
		if VERBOSE:
			parts = filename.split('/')
			print '-- parseSNAlignOutput: Input not found,', parts[len(parts)-1]
		return results
	
	qseq = ''; tseq = ''; idenseq = ''
	# read the alignment line by line, ignoring lines that don't
	# have a query-sbjct pattern as in stretcher output format
	cl = 0
	nnl = 0	# number of new lines seen
	pat = '^ *[\S\_]+ *\S? ([a-zA-Z-]+).*$'
	# print 'RUNNING------------------------------------'
	while cl < len(data):
		# get current line
		# print "line:", data[cl].strip("\n")
		# check if it starts with name of seq1 header
		if re.match('#.*', data[cl]) or re.match('^( *\d+)+$', data[cl]):
			# print 'TOSS:', data[cl]
			pass
		# elif re.match('^ *\S+.*$', data[cl]):
		elif re.match('^ *[\S\_]+ *\S? [a-zA-Z-]+.*$', data[cl]):
			# first parse out the query line
			# print "qline:", re.split(pat, data[cl])
			[crap1, qsubseq, crap2] = re.split(pat, data[cl])
			# print "qsub:", crap1, qsubseq
			qseq += qsubseq
			
			# how long is name part?
			rowheadlen = len(data[cl]) - len(qsubseq) - 1
			# print 'rowheadlen', rowheadlen
			
			# now get the alignment indicator line
			cl += 1
			line = data[cl].strip("\n")
			# strip leading spaces
			line = line[rowheadlen:len(line)]
			# print 'line \|-'+line+'-\"'
			idenseq += line
			# print 'len idenseq', len(idenseq)
			# then get the subject line
			cl += 1
			line = data[cl].strip("\n")
			# print "t curr:", line
			if re.match('^ *[\S\_]+ *\S? [a-zA-Z-]+.*$', data[cl]):
				[crap1, tsubseq, crap2] = re.split(pat, data[cl])
				#print "t got line:", crap1, tsubseq
				tseq += tsubseq
				# print 'tsub:', crap1, tubseq
				# print 'seqs \|-'+tsubseq+'-\"'
				cl += 3
			else:
				if VERBOSE: print '-- parseSNAlignOutput: Cannot parse output '+filename
				return results
				
		cl += 1
	
	# now we should have both full length sequences, resulting
	# from the alignment
	results['qseq'] = qseq
	results['tseq'] = tseq
	results['idenseq'] = idenseq
	
	# print 'lqs', len(qseq), 'lts', len(tseq), 'ids', len(idenseq)
	# print '-------------------------------------------------------'
	# store header/results info regarding the alignment
	#results = data
	
	foundtheline = False
	cl = 0
	parsed = list()
	while not foundtheline and cl < len(data):
		#parlen = re.split('# Length: *(\d+).*\n', results[line])
		parsed = re.split('# Identity: *(\d+)/(\d+)(.*)', data[cl])
		if len(parsed) > 1:
			foundtheline = True
		else:
			cl += 1
	if cl < len(data) - 2:
		(crap, iden, lenwithgaps, crap2, crap3) = parsed
		# store results
		results['iden'] = int(iden) # in stretcher this does not count N vs N
		results['lenwithgaps'] = int(lenwithgaps)
		results['niden'] = int(results['lenwithgaps'] - results['iden'])
		results['pi'] = results['niden'] / float(results['lenwithgaps'])
		
		# now parse gaps from later line
		cl += 2
		anygaps = re.split('# Gaps: *(\d+)/(\d+)(.*)', data[cl])
		if len(anygaps) > 1:
			(crap, ngaps, crap2, crap3, crap4) = anygaps
			results['ngaps'] = int(ngaps)
			#print "ngaps", ngaps
		cl += 1
		#print 'line', data[cl]
		(crap,score,crap) = re.split('# Score: *(-?\d+)*', data[cl])
		#print >> sys.stderr, 'score!', score
		results['score'] = int(score)
	
	# now correct lenwithgaps and ngaps if not indels
	if not indels:
		nres = fixIndelHelper([qseq, tseq], terminals)
		results['lenwithgaps'] = nres['lenwithgaps']
		results['ngaps'] = nres['groupGaps']
		if VERBOSE: print 'initial sub', nres['lenwithgaps'], 'by', results['iden']
		results['niden'] = nres['lenwithgaps'] - results['iden']
		results['pi'] = ut.divna(results['niden'], results['lenwithgaps'])
	
	# now correct for missing nucleotides indicated by the missing='x' argument
	if missing != '':
		missingindices = {}
		# which positions are missing values?
		qs = results['qseq'].upper()
		for i in range(len(qs)):
			if qs[i]==missing: missingindices[i] = 1
		# repeat for tseq
		ts = results['tseq'].upper()
		for i in range(len(ts)):
			if ts[i]==missing: missingindices[i] = 1
		missingidxlist = missingindices.keys()
		
		# go through idenseq to determine how many missing positions are matches (N vs N)
		isq = results['idenseq']
		missingiden = 0  
		missingniden = 0 # how many are mismatches (N vs A)
		missinggap = 0   # how many are gaps (N vs -)
		for i in range(len(isq)):
			if i in missingidxlist: 
				# is this position an identity?
				if isq[i] == ':':
					# so this must reflect A vs A; let remain
					pass
				elif isq[i] == ' ' and (ts[i]==missing and qs[i]==missing):
					missingiden += 1
				elif isq[i] == ' ' and (ts[i]=='-'and qs[i]==missing or qs[i]=='-'and ts[i]==missing):
					# N vs gap
					# print 'N vs gap', ts[i], qs[i]
					missinggap += 1
				elif isq[i] == ' ' and (ts[i]==missing or qs[i]==missing):
					# important that we do not discount N vs gap, since N is a real nt,
					# we just don't know the particular identity
					# print 'N vs nt', ts[i], qs[i]
					# N vs nucleotide
					missingniden += 1
		
		# apply correction
		# results['lenwithgaps'] -= len(missingidxlist)
		# results['ngaps'] -= missinggap
		# results['iden'] = isq.count(':')
		correction = missingniden
		if VERBOSE: print 'correct for %s N vs nt mismatches'%(correction)
		# results['iden'] -= missinggap # DO count N vs - as a mismatch
		# results['niden'] = results['lenwithgaps'] - results['iden']
		results['niden'] -= correction # don't count N vs A mismatches
		results['pi'] = results['niden'] / float(results['lenwithgaps'])
		
	
	if VERBOSE: print 'indels', indels, 'lenwithgaps', results['lenwithgaps'], \
		'ngaps', results['ngaps'], 'iden', results['iden'], 'niden', results['niden']
	return results

def parseClustalOutput(filename, VERBOSE=False, indels=False, terminals='-', missing='N', keepgaps=True, readTree=False, rmTree=True):
	# print 'got indels', indels
	# results dictionary
	results = {'seq':{}, 'idenseq':'', 'pi':ut.nastr, 'iden':ut.nastr, 'pairwisepid':ut.nastr,
	           'niden':ut.nastr, 'lenwithgaps':ut.nastr, 'ngaps':ut.nastr, 
	           'score':ut.nastr, 'file':filename}
	
	# get the alignment file data
	data = []
	try:
		file = open(filename)
		data = file.readlines()
		file.close()
	except IOError: return results
	
	if len(data): data.pop(0)
	
	# determine number of sequences in alignment
	dstr = ''.join(data)
	i = 0
	while dstr[i:i+1] == '\n':
		dstr = dstr[i+2:len(dstr)]
		i += 1
	
	groups = dstr.split('\n\n')
	tmp = groups[:]
	map(lambda x: x.strip('\n'), tmp)
	
	sg = groups[0].split('\n')
	# potential null entry
	if sg[len(sg)-1] == '': sg.pop(len(sg)-1)
	nseqs = max(0,len(sg)-1) # consensus sequence
	# print 'numseqs', nseqs
	# print 'stuff', sg
	# parse out labels
	ids = ['' for i in range(nseqs)]
	pat = '(^.+) +(\S+)$'
	for i in range(nseqs):
		(crap, ids[i], crap, crap) = re.split(pat, sg[i])
		ids[i] = ids[i].strip(' ')
	# print 'ids', ids
	
	cl = 0
	iden = 0
	lenwithgaps = 0
	seqs = ['' for i in range(nseqs)]
	consensus = ''
	# the key to clustal parsing is the *** lines
	# print 'len data', len(data), 'nseq', nseqs
	lheader = 0
	while nseqs > 1 and cl < len(data)-nseqs:
		pat = '^.+ +(\S+)$'
		if re.match(pat, data[cl]):
			# print 'data line', data[cl]
			
			# argh = data[cl].split(' ')
			
			if not lheader:
				(crap, seq, crap) = re.split('^(.+ +)[ATCGatcgnN\-]+\n$', data[cl])
				# print seq
				lheader = len(seq)
				# lheader = len(data[cl])-len(seq)-1
			
			for i in range(nseqs):
				(crap, seq, crap) = re.split(pat, data[cl])
				seqs[i] += seq
				# print 'seq', seq
				cl += 1
				
			stars = data[cl][lheader:-1]
			# print 'star line "%s"'%(stars)
			consensus += stars
			cl += 1
			
		cl += 1
	
	# print 'check', len(consensus), map(len,seqs)
	# for zz in range(len(consensus)):
	# 	if consensus[zz] == ' ': print 'line', zz
	
	# except IndexError:
	# 	print 'index error', len(data), nseqs, cl
	# 	print ids
	# 	data
	# 	sys.exit()
	
	seqdct = {}
	for i in range(nseqs):
		tmp = ''
		if keepgaps: tmp = seqs[i]
		else: tmp = seqs[i].replace('-','')
		# print 'seq', tmp
		fa = {'id':ids[i], 'altid':'', 'info': '', 'seq': tmp}  # no gaps
		seqdct[ids[i]] = fa
	results['seq'] = seqdct
	results['idenseq'] = consensus
	
	# print 'len consensus', consensus
	# this needs to change for indels
	nres = fixIndelHelper(seqs, terminals)
	lenwithgaps = len(results['idenseq'])
	ngaps = nres['pointGaps']
	if not indels:
		lenwithgaps = nres['lenwithgaps']
		ngaps = nres['groupGaps']
	
	iden = consensus.count('*')
	
	
	# ---------------------------------------------------------
	# if readTree:
	# 	ofile = filename[0:len(filename)-4]
	# 	if not os.access(ofile+'.dnd', os.F_OK):
	# 		os.system('clustalw '+filename+' &> /dev/null')
	# 	results['tree'] = readTree(ofile+'.dnd')
	# 	if rmTree: os.system('rm "'+ofile+'.dnd"')
	# ---------------------------------------------------------
	
	if len(ids)>0: results['ancestor'] = ids[0] # top row of alignment is most ancestral
	else: results['ancestor'] = ''
	
	
	
	results['iden'] = iden
	results['lenwithgaps'] = lenwithgaps
	results['niden'] = ut.subna(lenwithgaps, iden)
	results['file'] = filename
	results['ngaps'] = ngaps
	results['pi'] = ut.divna(results['niden'], float(results['lenwithgaps']))
	
	# also compute average pairwise identity
	sk = seqdct.keys()
	ds = []
	for i,j in ut.makeCombinations(len(sk),2):
		ds += [ut.subna(1,ut.divna(ut.hamdist(seqdct[sk[i]]['seq'], seqdct[sk[j]]['seq']), lenwithgaps))]
	results['pairwisepid'] = ut.avg(ds)
	
	if missing != '':
		# correct for N-N identities and N vs nt mismatches
		# however this is a multiple alignment. if we have a column e.g. AATNNN
		# which involves 6 taxa, although there are missing data, we KNOW there is
		# a true mismatch, so this position should be considered a segretating site. 
		# it is only in pairwise comparisons where we have N vs A or N vs T 
		# where Pi will be affected
		# this means that we only want to correct niden at positions that have
		# a single allele and Ns - only 2 characters. do NOT consider this segregating.
		
		missingindices = {}
		# find positions with missing value
		for k in seqdct.keys():
			tmp = seqdct[k]['seq']
			for i in range(len(tmp)):
				if tmp[i] == missing: missingindices[i] = 1
				
		# for each position with missing data, how many involve a single nucleotide and N?
		missingniden = 0
		missingidxlist = missingindices.keys()
		for i in missingidxlist:
			pool = {}
			for k in seqdct.keys(): 
				# print 'test', k, seqdct[k]['seq']
				try:
					pool[seqdct[k]['seq'][i]] = 1
				except IndexError:
					pass
			ks = pool.keys()
			# important that we do not discount N vs gap, since N is a real nt,
			# we just don't know the particular identity
			if len(ks) == 2 and missing in ks and '-' not in ks: missingniden += 1
		
		# apply correction
		correction = missingniden
		results['niden'] -= correction
		results['pi'] = ut.divna(results['niden'], float(results['lenwithgaps']))
		
	
	if VERBOSE: print 'indels', indels, 'lenwithgaps', lenwithgaps, 'ngaps', ngaps, 'iden', iden, 'niden', results['niden']
	return results

# -----------------------------------------------
# -----------------------------------------------

def sliceAlignment(chrom,start,stop,maindir, alnnfile='mavid.mfa'):
	
	res = readFasta(ut.slash(maindir)+ut.slash(str(chrom))+alnnfile)
	
	tab = []
	for k in res.keys():
		tab += [[k,res[k]['seq'][start:stop]]]
	
	def constar(lst):
		setnt = lst[0]
		for i in range(1,len(lst)):
			if lst[i]!=setnt: return ' '
		return '*'
	
	
	def consensus(lst):
		af = {'A':0,'T':0,'C':0,'G':0,'N':0,'-':0}
		for nt in lst: af[nt.upper()]+=1
		foo = [(af[nt],nt) for nt in af.keys()]
		foo.sort()
		return foo[-1][1]
	
	
	# add star
	tt = ut.transpose(tab)
	try:
		foo = zip(*tt[1:][0])
		star = ''.join(map(constar, foo))
		cons = ''
		try: cons = ''.join(map(consensus, foo))
		except IndexError: pass
		tab += [['match',star]]
		tab += [['cons ',cons]]
	except IndexError:
		pass
	return tab


def parseBlastFile(filename, headings=['qid', 'tid', 'pct', 'hitlen', 'mismatches', 'gapopens', 'qstart', 'qstop', 'tstart', 'tstop', 'evalue', 'bitscore'], upto=None):
	Qgenes = dict()
	
	# read in the blast results file
	file = open(filename)
	count = 0 
	for line in file:
		if upto:
			if count >= upto[1]: break
			elif count < upto[0]: 
				count += 1
				continue
		# if count % 100000 == 0: print ' - reached %s'%(count)
		
		rowlist = line[:-1].strip('\t').split('\t')
		
		# print 'testing'
		# print line[:-1].split('\t')
		# print rowlist
		
		
		# create a rowdict from rowlist
		rowdict = dict()
		try:
			rowdict[headings[0]] = rowlist[0]
			rowdict[headings[1]] = rowlist[1]
			rowdict[headings[2]] = float(rowlist[2])
			rowdict[headings[3]] = int(rowlist[3])
			rowdict[headings[4]] = int(rowlist[4])
			rowdict[headings[5]] = int(rowlist[5])
			rowdict[headings[6]] = int(rowlist[6])
			rowdict[headings[7]] = int(rowlist[7])
			rowdict[headings[8]] = int(rowlist[8])
			rowdict[headings[9]] = int(rowlist[9])
			rowdict[headings[10]] = float(rowlist[10])
			rowdict[headings[11]] = float(rowlist[11])
		except ValueError: print '- parsing error:', rowlist[:1]
		
		#printDict(rowdict)
		#print
		
		qid = rowdict[headings[0]]
		#print "Current gene", qid
		tid = rowlist[1]
		
		# append the rowdict onto a list for the query gene
		# qid = rowdict[headings[0]]
		# if Qgenes.has_key(qid): Qgenes[qid] = Qgenes[qid]+[rowdict]
		# else: Qgenes[qid] = [rowdict]
		
		
		# append the rowdict onto a list for the query gene
		if Qgenes.has_key(qid):
			temp = Qgenes[qid]
			# temp.append(rowdict)
			temp += [rowdict]
			Qgenes[qid] = temp
		else:
			Qgenes[qid] = [rowdict]
		
		count += 1
	
	file.close()
	return Qgenes


def parseBlastTable(table, headings=['qid', 'tid', 'pct', 'hitlen', 'mismatches', 'gapopens', 'qstart', 'qstop', 'tstart', 'tstop', 'evalue', 'bitscore']):
	
	Qgenes = dict()
	for rowstr in table:
		rowlist = rowstr.strip("\n").strip("\t").split("\t")
		# create a rowdict from rowlist
		rowdict = dict()
		rowdict[headings[0]] = rowlist[0]
		rowdict[headings[1]] = rowlist[1]
		rowdict[headings[2]] = float(rowlist[2])
		rowdict[headings[3]] = int(rowlist[3])
		rowdict[headings[4]] = int(rowlist[4])
		rowdict[headings[5]] = int(rowlist[5])
		rowdict[headings[6]] = int(rowlist[6])
		rowdict[headings[7]] = int(rowlist[7])
		rowdict[headings[8]] = int(rowlist[8])
		rowdict[headings[9]] = int(rowlist[9])
		rowdict[headings[10]] = float(rowlist[10])
		rowdict[headings[11]] = float(rowlist[11])
		
		#printDict(rowdict)
		#print
		
		qid = rowdict[headings[0]]
		#print "Current gene", qid
		tid = rowlist[1]
		
		# append the rowdict onto a list for the query gene
		if Qgenes.has_key(qid):
			temp = Qgenes[qid]
			temp.append(rowdict)
			Qgenes[qid] = temp
		else:
			Qgenes[qid] = [rowdict]
			#print Qgenes['YGR102C']
			#print
			
	return Qgenes

def parseBlastOutput_old(filename):
	"""Given the path to a bl2seq alignment results file, this
	method will read the file, and return a double, the first entry
	is the query sequence alignment and the second the target sequence.
	
	NB these sequences will be the actual blast sequence, not the input sequence given to bl2seq."""
	
	# get the alignment file data
	file = open(filename)
	data = file.readlines()
	file.close()
	
	# results dictionary
	blresults = {'qseq':"", 'tseq':"", 'niden':-1, 'lenwithgaps':-1, 'ngaps':-1, 'file':filename}
	
	qseq = ""
	tseq = ""
	# read the alignment line by line, ignoring lines that don't
	# have a query-sbjct pattern as in a typical blast output format
	cl = 0
	nnl = 0	# number of new lines seen
	while cl < len(data) and nnl < 3:
		# get current line
		line = data[cl]
		# check if it starts with "Query:"
		# also check that it's the first alignment in the file
		# as an indicator, three newlines in a row is used - find three and we're done
		if re.match('Query: .*', line):
			# first parse out the query line
			[crap1, qsubseq, crap2] = re.split('Query: \d+ +(.*) \d+', line)
			qseq += qsubseq
			# then get the subject line
			cl += 2
			line = data[cl]
			if re.match('Sbjct: .*', line):
				[crap1, tsubseq, crap2] = re.split('Sbjct: \d+ +(.*) \d+', line)
				tseq += tsubseq
			else:
				print "ERROR!"
			# reset nnl
			nnl = 0
		elif cl > 10 and line == '\n':
			nnl += 1
		
		cl += 1
	# now we should have both full length sequences, resulting
	# from the alignment
	
	blresults['qseq'] = qseq
	blresults['tseq'] = tseq
	
	# store header/results info regarding the alignment
	results = data
	
	foundtheline = False
	line = 0
	parsed = list()
	while not foundtheline and line < len(results):
		parsed = re.split(' Identities = (\d+)/(\d+)(.*)', results[line], 1)
		if len(parsed) > 1:
			foundtheline = True
		else:
			line += 1
	if not line == len(results):
		(crap, niden, lenwithgaps, checkforgaps, crap2) = parsed
		# store results
		blresults['niden'] = int(niden)
		blresults['lenwithgaps'] = int(lenwithgaps)
		# see if there were gaps
		blresults['ngaps'] = 0
		anygaps = re.split('Gaps = (\d+)/.*', checkforgaps)
		if len(anygaps) > 1:
			blresults['ngaps'] = int(anygaps[1])
	
	#return [qseq, tseq]
	return blresults

def parseBlastOutput(filename, method='bl2seq'):
	"""Given the path to a bl2seq alignment results file, this
	 method will read the file, and return a double, the first entry
	 is the query sequence alignment and the second the target sequence.
	 
	 NB these sequences will be the actual blast sequence, 
	 not the input sequence given to bl2seq."""
	
	
	# get the alignment file data
	file = open(filename)
	data = file.readlines()
	file.close()
	
	# pattern for matching nonidentities
	panpat = re.compile('[ +.]')
	
	qseq = ''
	tseq = ''
	editdistance = list() # |vector| = alignment length, binary valued
	# read the alignment line by line, ignoring lines that don't
	# have a query-sbjct pattern as in a typical blast output format
	
	# results dictionary
	blresults = {'qpos':ut.nastr, 'tpos':ut.nastr, 'qseq':'', 'tseq':'', 
	           'idenseq':'', 'pi':ut.nastr, 'iden':ut.nastr,
	           'niden':ut.nastr, 'lenwithgaps':ut.nastr, 'ngaps':ut.nastr, 
	           'score':ut.nastr, 'file':filename, 'expect':ut.nastr, 'qstrand':'', 'tstrand':''}
	
	cl = 0
	nnl = 0 # number of new lines seen
	while cl < len(data) and nnl < 3:
		# get current line
		line = data[cl]
		# print 'X',line
		# check if it starts with "Query:"
		# also check that it's the first alignment in the file
		# as an indicator, three newlines in a row is used - find three and we're done
		
		if re.match('.*Score.*', line):
			crap, score, expect, crap = re.split('.*Score += +(\d+).+bits.+,.+Expect +=(.+)$', line)
			blresults['score'] = int(score)
			try:
				blresults['expect'] = float(expect)
			except ValueError:
				blresults['expect'] = float('1'+expect.strip(' '))
		
		if re.match('.*Strand.*', line):
			crap, std1, std2, crap = re.split('.*Strand = (.+) \/ (.+)$', line)
			blresults['qstrand'] = std1
			blresults['tstrand'] = std2
		
		if re.match('Query: .*', line):
			# first parse out the query line
			# print 'line', line, re.split('Query: (\d+)( +)(.*) (\d+)', line)
			# print 'line', line, re.split('Query: (\d+)( +)(.*) (\d+)', line)
			
			# [crap,qpos,qsubseq,tpos] = line.strip('\n').split(' ')
			
			[crap1, qpos, spacer, qsubseq, tpos,crap2] = re.split('Query: (\d+)( +)(.*) (\d+)', line)
			
			qseq += qsubseq
			
			# if blresults['qpos'] == ut.nastr: blresults['qpos'] = int(qpos)
			
			# now get the number of identities in this alignment subsequence
			cl += 1
			line = data[cl][len(spacer):len(line)].strip('\n')
			# print "line:", line
			# print "spacer cl", cl, len(spacer), len(qsubseq), len(line), line,
			# for each position, did the 2 sequences align?
			for i in range(0,len(line)):
				if re.match(panpat, line[i]): editdistance.append(0)
				else: editdistance.append(1)
			# then get the subject line
			cl += 1
			line = data[cl]
			if re.match('Sbjct: .*', line):
				# print 'test!!!', re.split('Sbjct: (\d+) +(.*) \d+', line)
				[crap1, tpos, tsubseq, crap2] = re.split('Sbjct: (\d+) +(.*) \d+', line)
				tseq += tsubseq
				# if blresults['tpos'] == ut.nastr: blresults['tpos'] = int(tpos)
			else:
				print "ERROR!"
			
			blresults['qpos'] = int(qpos)
			blresults['tpos'] = int(tpos)
			
			# reset nnl
			nnl = 0
		elif method == 'bl2seq' and cl > 10 and line == '\n':
			nnl += 1
		
		cl += 1
	# now we should have both full length sequences, resulting
	# from the alignment
	
	blresults['qseq'] = qseq
	blresults['tseq'] = tseq
	blresults['editdistance'] = editdistance
	
	# store header/results info regarding the alignment
	results = data
	
	foundtheline = False
	line = 0
	parsed = list()
	while not foundtheline and line < len(results):
		parsed = re.split(' Identities = (\d+)/(\d+)(.*)', results[line], 1)
		if len(parsed) > 1:
			foundtheline = True
		else:
			line += 1
	if not line == len(results):
		(crap, iden, lenwithgaps, checkforgaps, crap2) = parsed
		# store results
		blresults['iden'] = int(iden)
		blresults['lenwithgaps'] = int(lenwithgaps)
		blresults['niden'] = int(blresults['lenwithgaps'] - blresults['iden'])
		blresults['pi'] = blresults['niden'] / float(blresults['lenwithgaps'])
		
		blresults['ngaps'] = 0
		anygaps = re.split('Gaps = (\d+)/.*', checkforgaps)
		if len(anygaps) > 1:
			blresults['ngaps'] = int(anygaps[1])
	
	#return [qseq, tseq]
	return blresults


def properStrand(gdict, outfile=sys.stdout, ISCDS=False, key='seq', plexy=False, critLen=1e6, minLen=90): 
	"""Some sequences may not reflect the actual CDS sequence, but
	are the reverse complement. this detects start/stop codons for each
	entry in an mfasta file and takes the reverse complement if needed.
	"""
	
	# indices are gene ids corresponding to improperly formed gene sequences
	badentries = list()
		
	print >> outfile, "Determining which sequences have",
	print >> outfile, "improper start and stop codons..."
	
	nstart = 0
	nstop = 0
	nflip = 0
	nmult3 = 0
	print >> outfile
	print >> outfile, "gene\treason"
	keys = gdict.keys()
	keys.sort()
	
	for k in keys:
		noteflip = 0
		#print k, gdict[k]
		s = gdict[k][key]
		
		first3 = s[0:3].lower()
		last3 = s[len(s)-3:len(s)].lower()
		
		if not plexy:
			if first3 == 'atg':
				nstart += 1
			else:
				print >> outfile, k,"\t> missing 5' seq - first 3 not start codon but:", first3
				noteflip += 1
				badentries.append(k)
		
			if last3 == 'tga' or last3 == 'tag' or last3 == 'taa':
				nstop += 1
			else:
				print >> outfile, k,"\t< missing 3' seq - last 3 not stop codon but:", last3
				noteflip += 1
				badentries.append(k)
		
			# check for multiple of 3
			if ISCDS:
				if len(s) % 3 == 0:
					nmult3 += 1
				else:
					print >> outfile, k, "\t< sequence not a multiple of 3"
					badentries.append(k)
		
			if len(s) > critLen:
				print >> outfile, k, '\t> sequence too long:', len(s)
				badentries.append(k)
			elif noteflip == 2:
				nflip += 1
				gdict[k][key] = reverseComplement(s).upper()
		else:
			# we won't use start/stop codon criteria to find a gene
			# rather just use minimum length
			if len(s) < minLen:
				print >> outfile, k, '\t> sequence too short:', len(s)
				badentries += [k]
				
		
	print >> outfile
	print >> outfile, "number of starts:", nstart
	print >> outfile, "number of stops:", nstop
	print >> outfile, "number of flips:", nflip
	if ISCDS:
		print >> outfile, "number of seqs multiple of 3:", nmult3
	
	print >> outfile, "total seqs:", len(gdict.keys())
	print >> outfile, "\n"
	
	return badentries

def properStrandFromFile(filename, savefilename="tossme.txt", ISCDS=False, key='seq'):
	"""Some sequences may not reflect the actual CDS sequence, but
	are the reverse complement. this detects start/stop codons for each
	entry in an mfasta file and takes the reverse complement if needed."""
	
	# get the alignment file data
	gdict = readFasta(filename)
	
	properStrand(gdict, savefilename, ISCDS, key)


def getExonIndices(cdsseq="", genseq="", buflen=15):
	
	""" 
	Given two sequences, presumably corresponding to the same protein coding
	locus, return a dictionary with 3 keys: exonCount (number), exonStarts (list),
	exonEnds (list)
	
	 you can recover exon indices by iterating through the cds sequence by nt,
	comparing the current base against that in the genomic sequence. Since these 2
	sequence are identical except for introns, the first cds base that does not
	match the current genomic base is the start of an intron. Increment the
	genomic index until you match the current cds sequence base again. repeat sic
	until end of cds sequence.
	
	You will quickly realize that there is a good chance (1/4) that any pair of
	bases will be the same, and so one cannot determine the number of exons by
	counting how many mismatches there are between the two sequences. I take a
	probabilistic approach, whereby we want to make sure we obtain all
	exons, at the cost of increasing the exonCount incorrectly. There will be a
	buffer of genomic sequence until you are pretty sure that the current
	subsequence is the next exon. Buffer length is specified as an argument. 
	The default is	a reasonable minimum length of an intron 
	(in mouse this is 20bp). 10 bases is used here as a default for yeast.
	(false exonCount increment every 1/4^10 exons)
	"""
	
	cds = cdsseq.upper()
	gen = genseq.upper()
	
	exonCount = 0
	exonStarts = list()
	exonEnds = list()
	
	NEWEXON = True
	NEWINTRON = True
	
	i = 1
	j = i
	while i <= len(cds) and j <= len(gen):
		#print "comparing:",i, cds[i-1], "-", gen[j-1], "...",
		
		if cds[i-1] == gen[j-1]:
			if NEWEXON:
		#		print "new exon"
				exonStarts.append(j-1)
				NEWEXON = False
		#	else:
		#		print "OLD exon"
			i += 1
			j += 1
			#print "\texon"
		elif cds[i-1] != gen[j-1]:
			if NEWINTRON:
		#		print "new intron"
				exonEnds.append(j-1)
				exonCount += 1
				NEWINTRON = False
		#	else:
		#		print
			
			#print "\t", i, "intron"
			
			# found an intron; now go through gen until we match cds again
			buffer = 0	# start a buffer
			itry = i	# pretend like we are incrementing i also
			jtry = i+1	# current intron index
			jbak = jtry	# backup
			JUSTFAILED = False
			while buffer < buflen and jtry <= len(gen):
				# maintain current buffer until mismatch
		#		print "\tbuffer length:", buffer, "/", buflen, "intron:", gen[jtry-1], "vs exon", cds[itry-1], "itry:", itry, "jtry:", jtry
				if gen[jtry-1] == cds[itry-1]:
					buffer += 1
					itry += 1
					jtry += 1
					jbak += 1
					JUSTFAILED = False
				else:
					buffer = 0
					itry = i
					if JUSTFAILED:
						jtry += 1
						#jbak += 1
					else:
						JUSTFAILED = True
						jtry = jbak
					
			# exit while loop either with a full buffer, or having reached
			# end of sequence
			# now we are sure the buffer is good, so drop back to exon before buffer
			i = itry - buffer
			# but we will be continuing where we just left off in gen seq
			j = jtry - buffer
		#	print "finished intron, i is ", i, "j is ", j, "len", len(gen)
			NEWEXON = True
			NEWINTRON = True
			if j > len(gen):
				#print "returning?"
				# special case
				#exonStarts.append(j-1)
		#		print "in the intron at end"
				return {'exonCount':exonCount, 
						  'exonStarts':exonStarts, 
						  'exonEnds':exonEnds}
		#	else:
		#		print "not returning..."
		#else:
		#	print
		#print "finishing current iteration with i:", i
	#print "out of while loop"
	# we've finished looping through entire sequence, ending with an exon
	# so return the annotations
	exonCount += 1
	exonEnds.append(len(gen))
	return {'exonCount':exonCount, 
			  'exonStarts':exonStarts, 
			  'exonEnds':exonEnds}

def printExonsInBED(dct={}):
	"""Takes a dictionary keyed by gene mapping to an exon dict of 3 keys: 
	exonCount (integer), exonStarts (list), exonEnds (list). Prints out
	each row as a gene followed by these three values.
	"""
	
	keys = dct.keys()
	keys.sort()
	# print header
	print "#name\texonCount\texonStarts\texonEnds"
	for g in keys:
		exonCount = dct[g]['exonCount']
		exonStarts = dct[g]['exonStarts']
		exonEnds = dct[g]['exonEnds']
		print g+"\t",
		print str(exonCount)+"\t",
		for i in exonStarts:
			print str(i)+",",
		print "\t",
		for i in exonEnds:
			print str(i)+",",
		print


def trimFasta(data, lengths={}, lengths3={}, trimDef=0, side='both', verbose=False, nclen=2000):
	"""Given a dictionary of fasta files, trim len bases from one or both sides.
	If a lengths dictionary is provided, trim each entry by the amount specified by
	lengths[name], otherwise trim by default value
	"""
	keys = data.keys(); keys.sort()
	for e in keys:
		seq = data[e]['seq']
		newseq = ''
		
		if side=='both': 
			try: trim5 = max(nclen - lengths[e], 0)
			except KeyError: trim5 = trimDef
			
			try: trim3 = max(len(seq) - lengths3[e], 0)
			except KeyError: trim3 = trimDef
			newseq = seq[trim5:len(seq)-trim3-nclen]
			if verbose: print e, '5:', trim5, '3:', len(seq)-trim3
		elif side=='5': 
			try: trim5 = max(len(seq) - lengths[e], 0)
			except KeyError: trim5 = trimDef
			newseq = seq[trim5:len(seq)]
			if verbose: print e, '5:', trim5
		elif side=='3': 
			try: trim3 = max(len(seq) - lengths[e], 0)
			except KeyError: trim3 = trimDef
			newseq = seq[0:len(seq)-trim3]
			if verbose: print e, '3:', len(seq)-trim3
		data[e]['seq'] = newseq
	return data
	

def trimFastaFile(faname, lengths={}, lengths3={}, trimDef=0, side='both', verbose=False, nclen=2000):
	"""given a dictionary of fasta files, trim len bases from both ends."""
	return trimFasta(readFasta(faname), lengths, lengths3, trimDef, side, verbose, nclen=nclen)



def blasthist(faname, column=0, delim='\t', nlines=-1, parseGenome=True, alnlen=None, alnpct=None, retlast=False):
	"""Make a categorical frequency historgram of a blast all results table for specified column."""
	
	dct = {}
	fh = open(faname)
	lineno = 1
	lastid = None
	for line in fh:
		if nlines > -1 and lineno >= nlines: break
		
		row = line[:-1].split(delim)
		lastid = row[0]
		
		geninfo = row[column]
		if parseGenome: 
			geninfo = geninfo.split('|')
			geninfo = geninfo[0]+'|'+geninfo[1]
			# print 'info', geninfo
		
		apct = ut.floatna(row[2])
		alen = ut.intna(row[3])
		
		BOTH = False
		if alnlen != None and alnpct != None: BOTH = True
		
		if (BOTH and alen >= alnlen and apct >= alnpct):
			# print 'test', row[:5]
			try: dct[geninfo] += 1
			except KeyError: dct[geninfo] = 1
			lineno += 1
		elif not BOTH and (alnlen == None and alnpct == None) or (alnlen != None and alnlen >= alen) or (alnpct != None and alnpct >= apct):
			try: dct[geninfo] += 1
			except KeyError: dct[geninfo] = 1
			lineno += 1
			
	fh.close()
	
	if retlast: return [dct.items(),lastid]
	else: return dct.items()
	# return ut.histogram(categorical=1, dat=dct.items())


# Statistical tests of selection / neutral theory
# -----------------------------------------------

def tajimaDHelper(K,S,n, verbose=False):
	# D statistic
	a1 = float(sum(map(ut.inv(1), range(1, n))))
	a2 = float(sum(map(ut.inv(2), range(1, n))))
	b1 = (n+1)/float(3*(n-1))
	b2 = (2*(n**2 + n + 3))/float(9*n*(n - 1))
	c1 = b1 - 1/a1
	c2 = b2 - (n+2)/(a1*n) + a2/a1**2
	e1 = c1/a1
	e2 = c2/(a1**2 + a2)
	
	# positive d => excess of common SNPs: 
	#   positive/directional selection increasing minor allele freq.
	# negative d => excess of rare (deleterious) SNPs: 
	#   purifying selection, bottleneck, etc
	d = K - S/a1 # difference in 2 estimates of Theta
	Vd = e1*S + e2*S*(S - 1)
	
	D = 0.0; pval = ut.nastr
	if Vd >= 0.0: 
		try: D = d/math.sqrt(Vd)
		except ZeroDivisionError: pass
	
	if D != ut.nastr: pval = ut.TajimaDistPDF(D, n)
	
	if verbose: print
	if verbose: print '%s = %s - %s (d = khat - S/a1)'%(d, K, S/a1)
	if verbose: print "  a1/a2", a1, a2, "b1/b2", b1, b2
	if verbose: print "  c1/c2", c1, c2, "e1/e2", e1, e2
	if verbose: print "  Vd", Vd
	if verbose: print '\n  Results:' 'D', D, 'p-value', pval
	
	return {'D':D, 'P':pval, 'a1':a1, 'Theta':S/a1}



def tajimaD(seqs, S=ut.nastr, K=ut.nastr, verbose=False, tempdir='', DELETE=True, REFRESH=True, name='gene', indels=False, haltIfMissing=False, terminals='-'):
	"""Returns the value of the Tajima D statistic for the input 
	list of sequences. Optionally takes S, the number of segregating
	sites from an alignment of the sequences. This is computed if
	not provided, using Clustalw.
	 
	If S and K are given, assumes they are across the locus, not per base!
	"""
	# verbose = 1
	if verbose: print '# Tajima D:', name
	
	# accept either fasta dictionary or list of sequences
	seqDct = seqs
	try: seqs.keys()
	except AttributeError: seqDct = list2fasta(seqs)
	
	n = len(seqs)
	if n < 2 or (haltIfMissing and (S==ut.nastr or K==ut.nastr)): 
		return {'D':ut.nastr, 'S':S, 'K':K, 'Var(K)':ut.nastr, 'n':n, 'Lk':ut.nastr, 'Ls':ut.nastr, 'p':ut.nastr, 'Pi_K':ut.nastr, 'Pi_S':ut.nastr, 'Theta':ut.nastr}
	
	# S: Total number of segregating sites
	Sf = S
	Ls = ut.nastr
	if S == ut.nastr:
		if verbose: print "\n  Generating multiple alignment of "+str(n)+" seqs with ClustalW..."
		results = multipleAlignment(seqDct, tempdir=tempdir, rmAln=DELETE, refresh=REFRESH, fn=name, 
		                            indels=indels, terminals=terminals)
		Sf = results['niden']
		Ls = results['lenwithgaps']
		if verbose: 
			print '  S', Sf, "segregating sites over length", Ls
			try:
				del results['qseq']; del results['tseq']; del results['idenseq']; del results['file']
			except KeyError: pass
			# for r in results.keys(): print r, results[r]
		if Sf == ut.nastr: return {'D':ut.nastr, 'S':Sf, 'K':K, 'Var(K)':ut.nastr, 'n':n, 'Lk':ut.nastr, 'Ls':Ls, 'p':ut.nastr, 'Pi_K':ut.nastr, 'Pi_S':ut.nastr, 'Theta':ut.nastr}
	if verbose: print '  S', Sf, 'seg sites over length', Ls
	
	# K: Average pairwise number of segregating sites
	khat = K; Lhat = ut.nastr; varK = ut.nastr
	if K == ut.nastr:
		if verbose: print '\n  Running pairwise alignments...'
		results = allPairsAlignment(seqDct, program='stretcher', tempdir=tempdir, 
		                            indels=indels, terminals=terminals)
		khat = results['full']['snp']
		Lhat = results['full']['lenwithgaps']
		varK = results['full']['V(snp)']
		if khat == ut.nastr:
			return {'D':ut.nastr, 'S':Sf, 'K':khat, 'Var(K)':varK, 'n':n, 'Lk':Lhat, 'Ls':Ls, 'p':ut.nastr, 'Pi_K':ut.nastr, 'Pi_S':ut.nastr, 'Theta':ut.nastr}
	if verbose: print '  Khat', khat, 'avg pairwise differences over length', Lhat
	
	R = tajimaDHelper(K=khat,S=Sf,n=n)
	D = R['D']; pval = R['P']; theta = R['Theta']; a1 = R['a1']
	
	return {'D':D, 'S':Sf, 'Theta':theta, 'K':khat, 'Var(K)':varK, 'n':n, \
	        'Lk':Lhat, 'Ls':Ls, 'p':pval, \
	        'Pi_K':ut.divna(khat,Lhat), 'Pi_S':ut.divna(Sf, Ls) \
	       }

def hkaHelper(name, seqA={}, seqB={}, Sa={}, Sb={}, D={}, verbose=False, tempdir='', DELETE=True, REFRESH=True, indels=False, useSegSites=True):
	"""Returns segregating sites for single locus of two species, and interspecies divergence."""
	import random
	
	locDat = {'Sa':ut.nastr, 'Sb':ut.nastr, 'D':ut.nastr, 
	          'La':ut.nastr, 'Lb':ut.nastr, 'Ld':ut.nastr,
	          'Na':ut.nastr, 'Nb':ut.nastr, 'Nd':ut.nastr}
	
	# accept either fasta dictionary or list of sequences
	# ---------------------------------------------------
	seqADct = seqA
	try: seqA.keys()
	except AttributeError: seqADct = list2fasta(seqA)
	seqBDct = seqB
	try: seqB.keys()
	except AttributeError: seqBDct = list2fasta(seqB)
	
	# Total number of segregating sites for species A
	# -----------------------------------------------
	Sfa = ut.nastr
	if useSegSites: Sfa = Sa[name]['S']
	else: Sfa = Sa[name]['K']
	
	Na = Sa[name]['Taxa']
	La = Sa[name]['Ls']
	if Sfa == ut.nastr and len(seqA):
		Na = len(seqA)
		if verbose: print "Generating multiple alignment of "+str(Na)+" seqs with ClustalW..."
		results = multipleAlignment(seqADct, tempdir=tempdir, rmAln=DELETE, refresh=REFRESH, fn=name, indels=indels)
		Sfa = results['niden']
		La = results['lenwithgaps']
		if verbose: 
			# print "There are", Sf, "polymorphic sites for", len(seqs), "sequences, with alignment."
			del results['qseq']; del results['tseq']; del results['idenseq']; del results['file']
			for r in results.keys(): print r, results[r]
		if Sfa == ut.nastr: return locDat
	if La == 0 or La == ut.nastr: Sfa = ut.nastr # must do this or will use 0 as number of seg sites
	
	# Total number of segregating sites for species B
	# -----------------------------------------------
	Sfb = ut.nastr
	if useSegSites: Sfb = Sb[name]['S']
	else: Sfb = Sb[name]['K']
	Nb = Sb[name]['Taxa']
	Lb = Sb[name]['Ls']
	if Sfb == ut.nastr and len(seqB):
		Nb = len(seqB)
		if verbose: print "Generating multiple alignment of "+str(Nb)+" seqs with ClustalW..."
		results = multipleAlignment(seqBDct, tempdir=tempdir, rmAln=DELETE, refresh=REFRESH, fn=name, indels=indels)
		Sfb = results['niden']
		Lb = results['lenwithgaps']
		if verbose: 
			# print "There are", Sf, "polymorphic sites for", len(seqs), "sequences, with alignment."
			del results['qseq']; del results['tseq']; del results['idenseq']; del results['file']
			for r in results.keys(): print r, results[r]
		if Sfb == ut.nastr: return locDat
	if Lb == 0 or Lb == ut.nastr: Sfb = ut.nastr
	
	# Divergence estimate between A and B
	# -----------------------------------
	# Randomly choose a sequence from each species
	Df = D[name]['S']
	Nd = 2
	Ld = D[name]['L']
	if Df == ut.nastr and len(seqA) and len(seqB):
		Nd = 2
		As = seqADct.values(); As.sort()
		Bs = seqBDct.values(); Bs.sort()
		randA = random.randrange(0, Na)
		randB = random.randrange(0, Nb)
		divResults = globalAlignment(As[randA], Bs[randB], program='clustalw', indels=indels)
		Df = divResults['niden']
		Ld = divResults['lenwithgaps']
	if Ld == 0 or Ld == ut.nastr or Ld == '': Df = ut.nastr; Ld = ut.nastr
	
	locDat['Sa'] = ut.floatna(Sfa); locDat['La'] = ut.floatna(La)
	locDat['Sb'] = ut.floatna(Sfb); locDat['Lb'] = ut.floatna(Lb)
	locDat['D'] = ut.floatna(Df); locDat['Ld'] = ut.floatna(Ld)
	locDat['Na'] = int(Na); locDat['Nb'] = int(Nb)
	locDat['Nd'] = 2
	return locDat

def HKA(specA={}, specB={}, Sa={}, Sb={}, D={}, verbose=False, tempdir='', DELETE=True, REFRESH=True, indels=False, useSegSites=True, twoTerm=False):
	"""
	HKA test of correlation between polymorphism and divergence
	for a single locus sampled in two species.
	
	Required input: dictionaries A, B, keyed by gene name, containing 
	either mfasta dictionaries or lists of sampled sequences at locus 
	'name', for two species A and B. 
	-------------------------------------------------------------------
	
	Polymorphism here is the number of segregating sites among a
	multiple alignment of samples of one species (intraspecific).
	
	Divergence is number of segregating sites between one randomly
	chosen sequence from each species (interspecific).
	"""
	import maths
	
	C = lambda n: float(sum(map(ut.inv(1), range(1, n))))
	C2 = lambda n: float(sum(map(ut.inv(2), range(1, n))))
	
	# Obtain statistics: Sa, Sb, and D for each common locus
	# ------------------------------------------------------
	namesA = specA.keys(); namesB = specB.keys()
	loci = ut.setIntersection(namesA, namesB)
	L = len(loci)
	
	# if no sequences, then parse loci from divergence data
	HAVESEQ = True
	if L == 0:
		# load gene names from stats files
		# in particular the divergence genes are all that's required
		loci = D.keys()
		L = len(D)
		HAVESEQ = False
		if L < 2: 
			return {'X2':ut.nastr, 'dof':ut.nastr, 'p':ut.nastr, \
			        'T':ut.nastr, 'f':ut.nastr, 'loci':{}}
	loci.sort()
	
	if verbose: print 'Computing HKA for', L, 'loci...'
	stats = [{} for i in range(L)]
	for i in range(L):
		locus = loci[i]
		# dict/list of sequences in species A/B for this locus
		locA = {}; locB = {}
		if HAVESEQ: locA = specA[locus]; locB = specB[locus]
		stats[i] = hkaHelper(locus, locA, locB, Sa, Sb, D, verbose, tempdir, DELETE, REFRESH, indels, useSegSites=useSegSites)
	
	# Estimate parameters Theta, T, and f
	# -----------------------------------
	# system of L + 2 equations
	Sa = [stats[i]['Sa'] for i in range(L)]
	Sb = [stats[i]['Sb'] for i in range(L)]
	Na = [stats[i]['Na'] for i in range(L)]
	Nb = [stats[i]['Nb'] for i in range(L)]
	D =  [stats[i]['D']  for i in range(L)]
	allDat = [loci, Sa,Na,Sb,Nb,D] # backup for saving info
	loci, Sa, Na, Sb, Nb, D = ut.filterStrFromLists(allDat, ut.nastr)
	# print Sa[0:10], Na[0:10], Sb[0:10]
	# internal hook for testing against the drosophila ADH locus data
	if False:
		loci = ['5p', 'adh']
		Sa = [30, 20]
		Sb = [0,0]
		Na = [11,11]
		Nb = [11,11]
		D = [78, 16]
	
	# revise L based on missing values
	L = len(loci)
	if verbose: print 'Corrected Number of loci:', L
	if L < 2: 
		if verbose: print '... need at least 2 loci to compute test statistic,'
		return {'X2':ut.nastr, 'dof':ut.nastr, 'p':ut.nastr, \
		        'T':ut.nastr, 'f':ut.nastr, 'loci':{}}
	
	# assume each loci has sample number of samples, Na
	# using the average number of samples over all loci
	CNa = [C(int(math.floor(Na[i]))) for i in range(L)]
	CNaHat = ut.avg(CNa)
	CNb = [C(int(math.floor(Nb[i]))) for i in range(L)] #C(int(math.floor(ut.avg(Nb))))
	CNbHat = ut.avg(CNb)
	# estimate Sum of thetas
	ThetaSum = sum(Sa)/CNaHat
	# estimate fhat
	fHat = sum(Sb)/(ThetaSum*CNbHat)
	if twoTerm: fHat = 1.0
	# estimate That
	That = (sum(D)/ThetaSum) - (1 + fHat)/2.0
	# estimate L-1 Theta_i values
	Theta = [0 for i in range(L)]
	for i in range(L):
		Theta[i] = (Sa[i] + Sb[i] + D[i]) / float(That + (1+fHat)/2.0 + CNa[i] + fHat*CNb[i])
	if twoTerm:
		for i in range(L):
			Theta[i] = (Sa[i] + D[i]) / float(That + (1+fHat)/2.0 + CNa[i])
	# compute Theta_L
	ThetaAlt = ThetaSum - sum(Theta)
	
	# check
	tolerance = .05 # error within 5%
	# print 'ThetaSum', ThetaSum, '?=', 'SumTheta', sum(Theta), '=', ThetaSum-sum(Theta)
	# print 'NOT TOLERANT', abs(ThetaSum-sum(Theta)),'>', tolerance*ThetaSum
	
	if abs(ThetaSum-sum(Theta)) > tolerance*ThetaSum or Theta[L-1] < 0.0:
		if verbose: 
			print 'ThetaSum', ThetaSum, '?=', 'SumTheta', sum(Theta)
			print 'NOT TOLERANT', abs(ThetaSum-sum(Theta)),'>', tolerance*ThetaSum
			print 'Theta[L-1] =', Theta[L-1], 'cf avg =', ut.avg(Theta)
		return {'X2':ut.nastr, 'dof':ut.nastr, 'p':ut.nastr, \
		        'T':ut.nastr, 'f':ut.nastr, 'loci':{}}
	
	# Compute Expectations
	# --------------------
	
	CNa2 = [C2(int(math.floor(Na[i]))) for i in range(L)] #C2(int(math.floor(ut.avg(Na))))
	CNb2 = [C2(int(math.floor(Nb[i]))) for i in range(L)] #C2(int(math.floor(ut.avg(Nb))))
	CNa2Hat = ut.avg(CNa2); CNb2Hat = ut.avg(CNb2)
	
	ESa = [Theta[i]*CNa[i] for i in range(L)]
	VSa = [ESa[i] + (CNa2[i]*Theta[i])**2 for i in range(L)]
	
	ESb = [fHat*Theta[i]*CNb[i] for i in range(L)]
	VSb = [ESb[i] + (CNb2[i]*Theta[i])**2 for i in range(L)]
	
	ED  = [Theta[i]*(That + (1+fHat)/2.0) for i in range(L)]
	VD  = [ED[i] + (Theta[i]*(1+fHat/2.0))**2 for i in range(L)]
	if twoTerm:
		ED  = [Theta[i]*That for i in range(L)]
		VD  = [ED[i] + Theta[i]**2 for i in range(L)]
	
	# Compute X2 random variable
	# --------------------------
	X2a = ut.sumna( [ ut.divna((Sa[i]-ESa[i])**2, VSa[i]) for i in range(L)] )
	X2b = ut.sumna( [ ut.divna((Sb[i]-ESb[i])**2, VSb[i]) for i in range(L)] )
	X2d = ut.sumna( [ ut.divna((D[i]-ED[i])**2, VD[i]) for i in range(L)] )
	X2 = X2a + X2b + X2d
	dof = 2*L-2
	if twoTerm: 
		X2 = X2a + X2d
		dof = L-1
	
	# print 'here:',X2,dof, X2a, X2b, X2d
	# print Theta
	pval = maths.ChiSquareCDF(X2,dof)
	
	if pval > .5: pval = 1.0-pval
	if verbose: print 'Result:', 'X2', X2, 'dof', dof, 'p-value', pval, 'T', That, 'f', fHat
	
	# organize statistics by locus
	# [loci, Sa,Na,Sb,Nb,D] = allDat
	results = {} # per locus
	for i in range(len(loci)):
		# print 'test', Sa[i], Sb[i], Na[i], Nb[i]
		results[loci[i]] = {
			'Sa':Sa[i], 'Na':Na[i], \
			'Sb':Sb[i], 'Nb':Nb[i], \
			'D':D[i], \
			'Theta':Theta[i], 'ESa':ESa[i], 'ESb':ESb[i], \
			'ED':ED[i]
		}
	
	return {'X2':X2, 'dof':dof, 'p':pval, 'T':That, 'f':fHat, 'loci':results, 'Sa':sum(Sa), 'Sb':sum(Sb), 'D':sum(D),
	        'X2a':X2a, 'X2b':X2b, 'X2d':X2d}


def divergenceGens(K=12.46, L=276.61, mu=1.84e-10, filename='', verbose=False, useK=True, useD=True, Dlow=-1.0, Dhigh=1.0):
	"""Returns the number of generations since the taxa separated evolutionarily.
	 filename is a matrix of genes x khat values
	 mu is mutation rate per base per generation
	"""
	
	low = Dlow
	hig = Dhigh
	d = K/float(L)
	diffs = []
	ells = []
	if filename:
		d = readDict(filename, key=True)
		for g in d.keys():
			K = ut.nastr
			L = ut.nastr
			if useK: 
				K = ut.floatna(d[g]['K'])
				L = ut.floatna(d[g]['Lk'])
			else: 
				K = ut.floatna(d[g]['S'])
				try:
					L = ut.floatna(d[g]['Ls'])
				except KeyError:
					L = ut.floatna(d[g]['L'])
			if useD:
				D = ut.floatna(d[g]['D'])
				if D < hig and D > low:
					diffs.append( ut.divna(K,L) )
			else: diffs.append( ut.divna(K,L) )
			ells.append(L)
		d = ut.avg(diffs)
	if verbose:
		if filename: print len(ut.filterstr(ut.nastr)(diffs)), 'loci'
		print 'Expected per site difference of', d, 'using K', K, 'L', L, 'mu', mu
		print 'average length L', ut.avg(ells)
	
	denom = float(2*mu)
	gens = ut.filterstr(ut.nastr)([ut.divna(diffs[i],denom) for i in range(len(diffs))])
	# avggens = d / denom # years
	# also want the variance
	# avggens = ut.avg(gens); 
	# vargens = ut.variance(gens)
	# print 'avggens', avggens, 'avg dist', avggens, 'var', vargens
	return ut.avg(gens), gens
	



def divergenceTime(n, r=2920, filename='', verbose=False):
	"""Only parameters are average K/L and mu and r
	r is generations per year
	"""
	
	if filename: n,alln = divergenceGens(K,L,mu,filename,verbose)
	if verbose: print 'Generations since divergence:', n
	return n/float(r)



def getSNS(filename, VERBOSE=False):
	M = readFromClustal(filename)
	
	# for each position, determine morphisms
	
	nonsynon = 0
	synon = 0
	
	gapflag = 0
	
	i = 0
	while i < len(M['idenseq']):
		# always remain in frame
		# print 'position', i
		# morph = {}
		codon = {}
		for j,so in M['seq'].items():
			seq = so['seq'].upper()
			# morph[seq[i]] = 1
		
			triplet = seq[i:i+3]
			# if frame == 0: 
			# elif frame == 1: triplet = seq[i-1:i+2]
			# elif frame == 2: triplet = seq[i-2:i+1]
		
			# this will need to be controlled for missing values, eg ANA
			if 'N' not in triplet and '-' not in triplet and len(triplet)==3: codon[triplet] = tr(triplet)
			elif '-' in triplet: 
				print 'problematic codon', triplet
				gapflag = 1
			
		# lmb = len(morph.keys())
		# if 'N' in morph: lmb -= 1
	
		# is it synon or nonsynon?
		# we have a multiple alignment, so for now just assess all codons and do an all pairs comparison
	
		pairs = codon.items() # codons with aas
		
		if VERBOSE: print i,'pairs', pairs
	
		for j in range(len(pairs)):
			for k in range(j):
				# first check amino acid identity
				if pairs[j][1] != pairs[k][1]:
					nonsynon += 1
					if VERBOSE: print i, 'NS: %s->%s, %s->%s'%(pairs[j][0], pairs[j][1], pairs[k][0], pairs[k][1])
				elif pairs[j][1] == pairs[k][1] and pairs[j][0] != pairs[k][0]:
					# then it is synonymous
					synon += 1
					if VERBOSE: print i, 'S: %s->%s, %s->%s'%(pairs[j][0], pairs[j][1], pairs[k][0], pairs[k][1])
	
		i += 3
	
	nonsynon += gapflag # count as a single event
	
	return {'S':synon, 'NS':nonsynon}





def loadRestrictionEnzymes(filename='/Users/simola/bin/simolacode/restrictionEnzymes.txt'):
	
	E = {}
	fh = open(filename)
	for line in fh:
		if line[0] == '#': continue
		row = line[:-1]
		# AatI         AGG/CCT                                  
		parts = filter(lambda x: len(x), map(lambda x: x.strip(), row.split(' ')))
		if not len(parts): continue
		
		name = parts[0]
		motifs = filter(lambda x: x != 'and', parts[1:])
		
		fm = []
		for m in motifs:
			m = m.replace('/','')
			
			# idx = m.find(' ')
			# if idx:
			# 	m = m[:idx]
			# 	# print 'here', m
			
			if re.match('.*\(\d+\).*', m):
				crap,keep,keep2,crap = re.split('(.*)\(\d+\)(.*)', m)
				m2 = ''
				if '(' not in keep: m2 += keep
				if '(' not in keep2: m2 += keep2
				m = m2
			if len(m) and '(' not in m and ')' not in m and len(m)>1:
				# print name, m
				fm += [m]
		E[name] = fm
	fh.close()
	
	return E


def enumMotif(seq,i=0,lst=['']):
	if i == len(seq)-1:
		# print 'final', i, seq[i]
		tmp2 = []
		for item in lst:
			for nt in iupac2nts(seq[i]):
				tmp2 += [item+nt]
		return tmp2
	else: 
		# print 'else'
		
		tmp2 = []
		for item in lst:
			for nt in iupac2nts(seq[i]):
				tmp2 += [item+nt]
		# print 'tmppy', tmp2
		return enumMotif(seq, i+1, tmp2)


def digest(seq, E={}):
	def check(seq,elist):
		ret = False
		for e in elist:
			# this should be improved!!!
			# for iue in enumMotif(e):
			iue = ''.join(map(iupac2nts,e))
			if len(e) == len(iue):
				pat = '.*(%s).*'%(iue)
				# print 'test', pat, seq, rc(seq)
				if re.match(pat,seq) or re.match(pat,rc(seq)): ret = True; break
			# else: print 'skipping', 
		return ret
	
	# print 'seq', seq
	return filter(lambda e: check(seq,E[e]), E.keys())


if __name__ == "__main__":
	print 'hello'
	
	res = parseClustalOutput('/Users/simola/Documents/Publications/paradoxus/work/alignSNPsToGenomes/alignments/CBS432.chr01_1968_A-G_ma.aln')
	sys.exit()
	
	
	print blast2Seqs('CCAAAGTTTTGTCAAGAACTAACCTCTGTTACCCAACGCCGACACTGACGCTGACGCTACCAA', 'CCAAAGTTTTGTCAAGAACTAACCTCTGTTACCCAACGCCGACACTGACGCTGACAAAAAAA', savepath='/Users/simola/', name='test', DELETE=False, params='-F F -p blastn -E 100000000000000')
	
	sys.exit()
	
	fn = '/Users/simola/Desktop/3primeUTR_YAL040C.aln'
	seq = fastaFromClustal(fn, terminals=True)
	print 'seq'
	for s in seq.keys():
		print s, '\n', seq[s]['seq']
	
	
	# fn = '/Volumes/Supra/working/Sparadoxus_genomes2/aligned_intron/YBR189W/YBR189W-Q32_3-intron_YBR189W-CBS432-intron_stretcher.txt'
	# parseSNAlignOutput(fn, indels=1, VERBOSE=1)
	# parseSNAlignOutput(fn, indels=0, VERBOSE=1)
	
	# fn = '/Volumes/Supra/working/Sparadoxus_genomes2/maligned_intron/YGR001C.aln'
	# parseClustalOutput(fn, nseqs=2, indels=1)
	# parseClustalOutput(fn, nseqs=2, indels=0)
	
	
	# switch = sys.argv.pop(1)
	# if switch == '-align':
	# 	file1 = readFasta(sys.argv.pop(1))
	# 	file2 = readFasta(sys.argv.pop(1))
	# 	key1 = file1.keys()[0]
	# 	key2 = file2.keys()[0]
	# 
	# 	stuff = stretcher(file1[key1]['seq'], file2[key2]['seq'])
	# 	print "ngaps", stuff['ngaps']
	# 	print 'niden', stuff['niden']
	# 	print 'lenwithgaps', stuff['lenwithgaps']
	# elif switch == '-orf':
	# 	filename = sys.argv.pop(1)
	# 	info = sys.argv.pop(1)
	# 	orfeome = ORFsFromFile(filename, info)
	# 	printFasta(orfeome)
	
	print 'Nothing to test 8('