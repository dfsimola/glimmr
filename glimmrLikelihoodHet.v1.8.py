#!/usr/bin/env python

GLIMMR_VERSION = 1.7

# NEW:
# C code for statePriS
# maxdepth default changed to 100
# default bufferlimit changed to 10000
# improved runtime for -d 1 (UNIQUE) only analysis

"""GlimmrLikelihood.py: likelihood computational backend for glimmrHet.py
 
 Version: %s
 
 Daniel F. Simola (simola@upenn.edu)
 Laboratory of Shelley L. Berger
 University of Pennsylvania
 October 2015
 
 Copyright (c) 2015, Daniel F. Simola and Shelley L. Berger, University of
 Pennsylvania.  All Rights Reserved.
 
 You may not use this file except in compliance with the terms of our
 License. You may obtain a copy of the License at XXX
 
 Unless required by applicable law or agreed to in writing, this
 software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 CONDITIONS OF ANY KIND, either express or implied.  See the License
 for the specific language governing permissions and limitations
 under the License.
"""%(GLIMMR_VERSION)



import os, re, sys, math, time, stat, time, random, subprocess
from operator import add,mul
from string import *

# GLOBAL VARIABLES AND FUNCTIONS
# -----------------------------------------------------------------
multichromtoken = '.'      # separator for --rest <chr1.chr2.....chrn> separation
nastr = '"NA"'             # represents missing value
SGr = {}                   # global current read pileup dictionary
PRELIK = {}                
URB = {'SE':[] , 'PE':[]}  # unique read buffer
DRB = {'SE':[], 'PE':[]}   # degenerate index read buffer
bufferlimit = 10000       # number of lines to buffer in memory (10MB)
PCR_HASH_WIDTH = PCR_HASH_NMAX = 5 # DO NOT CHANGE; FIXED IN FUNCTION



# helper function for pileup updating
def delete(x):
	global SGr, PRELIK
	del SGr[x]
	del PRELIK[x]
	return 1


def fetchLines(handler,list_of_byte_len):
	"Random access to lines in the read map files"
	for nbytes,nline in sorted(list_of_byte_len):
		handler.seek(nbytes)
		for x in xrange(nline): yield handler.next() # return a full line
		# yield handler.read(nline)


def pipeit(s, nl=0, pipe=sys.stdout):
	pipe.write(s+(nl==1 and '\n' or '')); sys.stdout.flush()

def isnan(x):
    return str(x) == str(1e400*0)

def unique(lst):
	"""Returns a subset of lst containing the unique elements."""
	d = {}
	for i in lst:
		if i not in d: d[i] = 0
		else: d[i] += 1
	e = d.keys(); e.sort()
	return e

def maxlenval(tup):
	mlv = 0
	for i in tup:
		x = len(str(i))
		if x>mlv:mlv = x
	return mlv

def fixspacing(item,nsp):
	v = str(item)
	if len(v) == nsp: return v
	delt = nsp-len(v)
	if delt < 0:
		v = v[0:nsp]
	else:
		for i in range(delt): v+=' '
	return v

def transpose(mat):
	"""Given a 2d array of arrays, will invert the dimensions."""
	newm = []
	for c in range(len(mat[0])):
		newrow = []
		for r in range(len(mat)):
			newrow.append(mat[r][c])
		newm.append(newrow)
	return newm


fixNull = lambda nastr: lambda x: x=='' and nastr or x
def flatten(sequence):
	def rflat(seq2):
		seq = []
		for entry in seq2:
			if seqin([entry]):
				seq.extend([i for i in entry])
			else:
				seq.append(entry)
		return seq
	
	def seqin(sequence):
		for i in sequence:
			## all sequences have '__contains__' in their dir()
			## parentheses present to aid commenting mid-condition
			if ('__contains__' in dir(i) and type(i) != str and type(i) != dict):
				return True
		return False
	
	seq = [sequence][:] # in case parameter isn't already a sequence
	while seqin(seq): seq = rflat(seq)
	return seq


def createdir(dirname):
	if not dirname: return 0
	if not os.access(dirname, os.F_OK):
		try: os.mkdir(dirname); return 1
		except OSError: pass
	return 0

def record(astr, filename, ioMethod='a'):
	if ioMethod == 'a':
		assert os.access(filename, os.F_OK), 'Cannot access '+filename
	
	pmrfh = open(filename, ioMethod)
	print >>pmrfh, astr
	pmrfh.close()
	
	return 1


def trimComments(dat):
	# parse comments
	ndat = []
	for line in dat:
		# only care about non commented lines
		if not re.match('(^\s*#.*$)|(^\s$)', line):
			# parse any commentary on this line
			cleanme = re.split('^(.*)#.*$', line)
			if len(cleanme) > 1: cleanme.pop(0)
			ndat.append(cleanme[0].strip('\n'))
	return ndat

def getFiles(indirfile, include=[], exclude=[]):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	return getContents(indirfile, typ='files', include=include, exclude=exclude)

def getDirectories(indirfile, include=[], exclude=[]):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	files = getContents(indirfile, typ='dirs', include=include, exclude=exclude)
	return map(slash, files)

def getContents(indirfile, typ='all', include=[], exclude=[]):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	keepfiles = []
	if len(include) == 0: include = ['']
	try:
		mode = os.stat(indirfile)[stat.ST_MODE]
		INISADIR = False
		if stat.S_ISDIR(mode):
			INISADIR = True
			indirfile = slash(indirfile)
		files = []
		if INISADIR and os.access(indirfile, os.F_OK):
			for f in os.listdir(indirfile):
				# print indirfile, 'vs', f
				# if f[0] != '.' and f[0] != '~' and f[len(f)-3:len(f)] not in exclude:
				if f[0] != '.' and f[0] != '~':
					# do includes
					for item in include:
						if re.match('.*%s.*'%(item), f):
							# confirm not a directory
							mode = os.stat(indirfile+f)[stat.ST_MODE]
							if typ=='all': 
								if stat.S_ISDIR(mode): files.append(indirfile+f)
								else: files.append(indirfile+f)
							elif typ=='files' and not stat.S_ISDIR(mode):
								files.append(indirfile+f)
							elif typ=='dirs' and stat.S_ISDIR(mode):
								files.append(indirfile+f)
							# else: print indirfile+f, 'is a directory'
			for f in files:
				inexclude = False
				for item in exclude:
					if re.match('.*%s.*'%(item), f): inexclude = True
				if inexclude == False:
					# confirm not a directory
					keepfiles.append(f)
					
					# try:
					# 	print 'testing', indirfile+f
					# 	mode = os.stat(indirfile+f)[stat.ST_MODE]
					# 	print 'mode', mode
					# except OSError:
					# 	print 'why?'
					# if typ=='all': 
					# 	if stat.S_ISDIR(mode): keepfiles.append(indirfile+f+'/')
					# 	else: keepfiles.append(indirfile+f)
					# elif typ=='files' and not stat.S_ISDIR(mode):
					# 	keepfiles.append(indirfile+f)
					# elif typ=='dirs' and stat.S_ISDIR(mode):
					# 	print 'got a dir'
					# 	keepfiles.append(indirfile+f+'/')
					# else:
					# 	print 'failure', stat.S_ISDIR(mode)
		elif not INISADIR: keepfiles.append(indirfile)
		else: sys.exit('cannot access dir '+indirfile)
	except OSError: pass
	return keepfiles

def readList(fname, header=False):
	crap,thelist,crap = readTable(fname, header=header)
	return flatten(thelist)

def printList(lst, title="", delim="\n", pipe=sys.stdout, file='', sort=False, newline='\n', ioMethod='w'):
	if file: pipe = open(file,ioMethod)
	if sort: lst.sort()
	if title: 
		lst = [title]+lst#print >> pipe, title
	for i in range(len(lst)):
		if i == len(lst)-1:
			print >> pipe, str(lst[i]),
		else:
			print >> pipe, str(lst[i])+delim,
	if newline: print >> pipe, newline,
	if file: pipe.close()

def readTable(filename, header=True, rownames=True, delim="\t", comments=True, keepKey=False):
	fh = open(filename)
	table = fh.readlines()
	fh.close()
	if comments: table = trimComments(table)
	
	tabHeader = []; rowIDs = []; data = []
	if header: 
		tabHeader = table.pop(0).strip('\n').split(delim)
		if not keepKey: tabHeader.pop(0)
		# else: tabHeader = table.pop(0).strip('\n').split(delim)
	for line in table:
		sline = line.strip('\n').split(delim)
		if rownames:
			rowIDs.append(sline.pop(0))
		data += [sline]
	return data, rowIDs, tabHeader

def printTable(tab, header=[], rows=[], delim="\t", newline='\n', file='', pipe=sys.stdout, ioMethod='w'):
	
	if file: pipe = open(file,ioMethod)
	
	if header:
		printList(header, delim=delim, pipe=pipe, newline=newline)
	if len(rows):
		assert len(rows) == len(tab), "Rows and tab must have same length."
		for ri in range(len(tab)):
			r = [rows[ri]]+tab[ri]
			printList(r, delim=delim, pipe=pipe, newline=newline)
	else:
		for r in tab:
			printList(r, delim=delim, pipe=pipe, newline=newline)
	
	if file: pipe.close()

def printFormattedTable(tab, header=[], rows=[], delim=" ", newline='\n', file='', pipe=sys.stdout, nastr='"NA"', ioMethod='w', colSpacing=[]):
	"""not good for tables > say 100 MB"""
	if file: pipe = open(file,ioMethod)
	if tab == [] and file: pipe.close(); return
	tab2 = tab[:]
	if rows:
		for i in range(len(tab)): tab2[i] = [rows[i]]+tab[i]
	
	# need to determine for each column the max length value and space to that
	collens = colSpacing[:]
	if not len(collens):
		# print 'printformattedtable: fixing col spacing'
		temp = transpose(tab2)
		if len(header):
			for i in range(len(temp)):
				r = temp[i]
				r.append(header[i])
				collens.append(maxlenval(r))
		else:
			for r in temp: collens.append(maxlenval(r))
	
	# now create a new tuple with proper spacing
	temp = []; sheadings = []
	if len(header):
		for i in range(len(header)):
			sheadings.append(fixspacing(header[i],collens[i]))
	
	for row in tab2:
		newrow = []
		for i in range(len(row)):
			val = row[i]==nastr and '.' or row[i]
			#if val==nastr: val = '.' # SAS missing value
			newrow.append(fixspacing(val,collens[i]))
		temp.append(newrow)
		
	# pipe = open(file,ioMethod)
	printTable(temp,header=sheadings,pipe=pipe,delim=delim,newline=newline)
	if file: pipe.close()

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

def readFasta(faname, usealtid=False, split='^>', idPat='', addins=[], TOUPPER=False, VERBOSE=False):
	"""Parses a multiple fasta file into a dictionary object keyed by fasta identifier. 
	Option to key the dictionary by an alternative id found in the fasta header: altid=True|False.
	If directory of several fasta files is provided as faname, loads them into single dictionary.
	"""
	fadict = {} # return dictionary of fasta entries
	faj = '' # load all information into single data string, faj
	files = getFiles(faname)
	
	for f in files:
		if VERBOSE: print '- Loading %s'%(f)
		fh = open(f)
		fa = fh.readlines()
		fh.close()
		# remove comment lines
		faj += ''.join(filter(lambda x: not re.match(' *#.*', x) and True, fa))
	
	
	# parse by entry (>)
	getentries = re.compile(split, re.MULTILINE)
	faentries = re.split(getentries, faj)
	
	# for some reason the first entry is the null string
	faentries.pop(0)
	
	# parse individual entries
	for entry in faentries:
		# first  has format >name other-header-info
		(header, seq) = entry.split("\n", 1)
		# trim off the '>' character
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
		
		# remove newlines so the sequence data is contiguous
		seq = re.sub("\n", "", seq)
		# remove terminal whitespace
		seq = seq.strip(' ')
		
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
		
		# build the entry
		# add a little check to see if there are repeat entries
		if fadict.has_key(theid):
			print theid, ": This key has already been used to store another fasta entry!"
		else:
			# build the entry for this id
			# allow choice of primary dict key
			if theid and altid and usealtid:
				#print "OK", theid, altid
				fadict[altid] = {'id':theid, 'altid':altid, 'info': info, 'seq': seq}
			else:
				fadict[theid] = {'id':theid, 'altid':altid, 'info': info, 'seq': seq}
	return fadict

def list2fasta(seqs=[], ids=[]):
	dct = dict()
	uid = 0
	for i in range(len(seqs)):
		theid = str(uid)
		try: theid = ids[i]
		except IndexError: pass
		dct[theid] = {'id':theid, 'seq':seqs[i], 'info':''}
		uid += 1
	return dct

def printFasta(fadict, file='', pipe=sys.stdout, key='seq', usealtid=False):
	"""Print a dictionary of fasta entries to a file in mfasta format"""
	keys = fadict.keys()
	keys.sort()
	
	if file: pipe = open(file,'w')
	
	for theid in keys:
		if 'altid' in fadict[theid] and not usealtid:
			header = '>'+theid+' '+fadict[theid]['altid']+' '+fadict[theid]['info']
		elif 'id' in fadict[theid]:
			header = '>'+theid+' '+fadict[theid]['id']+' '+fadict[theid]['info']
		else:
			header = '>'+theid+' '+fadict[theid]['info']
		print >> pipe, header
		print >> pipe, fadict[theid][key]
		print >> pipe
	
	if file: pipe.close()

def fasta2Phylip(dct={}, file='phylipi.txt', infile='', verbose=False):
	"""Convert a fasta file into Phylip interleaved format."""
	# print infile
	d = dct
	
	if infile: d = readFasta(infile)
	
	names = d.keys()
	# fix length
	maxlen = 10+1
	shorts = ['' for i in range(len(names))]
	for i in range(len(names)): 
		shorts[i] = fixspacing(names[i], maxlen)
	blank = fixspacing('',maxlen)
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


eq = lambda s: lambda i: str(s)!=str(i)
filterstr = lambda s: lambda v: filter(eq(s), v)

def argmax(v):
	"""Return argmax of a list"""
	mi = 0
	mv = v[0]
	for i in range(1,len(v)):
		if float(v[i]) > float(mv):
			mv = float(v[i])
			mi = i
	return mi

def slash(astr):
	if astr[len(astr)-1] !='/': astr+='/'
	return astr

def makeContigs(intervals):
	intervals.sort()
	# print intervals[0:100]
	# combine overlaps
	coverage = 0 # total bases mapped
	contigs = []
	curr = 1
	while curr < len(intervals):
		bak = curr
		start = intervals[curr-1][0]; stop = intervals[curr-1][1]
		
		while curr < len(intervals) and intervals[curr][0] <= stop:
			if intervals[curr][1] > stop: stop = intervals[curr][1]
			curr += 1
		
		contigs += [(start,stop)]
		coverage += abs(start-stop)
		curr += 1
	# print 'contigs', contigs
	return {'contigs':contigs, 'coverage':coverage}


roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))
floatna = lambda x: x == nastr and nastr or float(x)
intna = lambda x: (x == nastr or x == '' or x == 'NA' and nastr) or int(x)

def factit(x):
	"""Iterative implementation of factorial function."""
	if x <= 1: return 1
	a = [0 for i in range(x+1)]
	a[0] = 1; a[1] = 1
	for i in range(2,x+1): a[i] = i*a[i-1]
	return a[x]

nchoosek = lambda n,k: factit(n)/(factit(k)*factit(n-k)) # binomial coefficient
prod = lambda lst: reduce(mul, lst) # product function

def avg(lst):
	try: return reduce(add, lst)/float(len(lst))
	except TypeError: return nastr

def avg0(lst):
	try: return reduce(add, lst)/float(len(lst))
	except TypeError: return 0

def wavg(v, w):
	if len(v) == 0: return ret
	if len(w) > 0 and len(w) != len(v): return ret
	return sum([v[i]*float(w[i]) for i in range(len(v))])

def harmonicMean(v, minval=1e-20):
	if len(v) == 0: return 0
	# prevent 0 division
	v = map(lambda x: x==0 and minval or x, v)
	return len(v) / sum(map(lambda x: 1./x, v))

def wHarmonicMean(v, w, minval=1e-20):
	if len(v) == 0: return 0
	if len(w) > 0 and len(w) != len(v): return 0
	# prevent 0 division
	v = map(lambda x: x==0 and minval or x, v)
	# return len(v) / sum(map(lambda x: 1./x, v))
	return sum(w) / sum(map(lambda x,y: y/x, v,w))

def geometricMean(lst):
	try: return reduce(mul, lst)**1/float(len(lst))
	except TypeError: return 0


ssd = lambda lst, mu: reduce(add, map(lambda x: (x-mu)*(x-mu), lst))
def variance(lst):
	try: return ssd(lst,avg(lst))/float((len(lst)-1))
	except TypeError: return nastr

def sd(lst):
	try: return math.sqrt(variance(lst))
	except TypeError: return nastr

def percentile(v=[], k=0.5):
	"""
	Return the value of array that corresponds to the kth percentile of the data.
	"""
	if len(v) <= 0: return 0.0
	temp = map(float, v)
	temp.sort()
	n = len(temp)
	idx = float(k)*(n+1)
	lidx = math.floor(idx)
	hidx = math.ceil(idx)
	if idx < 1.0:
		return float(temp[0])
	elif idx > n-1:
		return float(temp[n-1])
	elif not (idx == lidx or idx == hidx):
		return float(temp[int(lidx)] + temp[int(hidx)])/2.0
	else:
		return float(temp[int(idx)])

median = lambda v: percentile(v,.5)
getQuants = lambda dat: ','.join([str(q)+':'+str(round(percentile(dat, q),4)) for q in [0.05, 0.25, 0.5, 0.75, 0.95]])

# list operations
head = lambda x: x[0]
tail = lambda x: x[len(x)-1]
car = lambda x: x[0]
cdr = lambda x: x[len(x)-1]

# arithmetic operations with missing values
subna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) - float(v)
multna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) * float(v)
divna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) / float(v)
roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))

# ref string may contain SNPs. Choose among possible allelic strings the one with fewest mismatches to obs string.
# This is a greedy approach, so we can do in linear time.
# --------------------------------------------------------

# this is slower
# def bestRefSeqHelper(rnt,ont):
# 	# if length == 1 then this is monoallelic
# 	if len(rnt)==1: return rnt
# 	for nt in rnt:
# 		if nt == ont: return nt
# 	# if we don't find a matching allele, choose the first (i.e. randomly)
# 	return rnt[0]
# 	
# bestRefSeq = lambda ref,obs: ''.join([bestRefSeqHelper(iupac2nts(ref[i]), obs[i]) for i in range(len(ref))])

# faster version
def bestRefSeq(ref,obs):
	bestref = ''
	
	ref = len(ref) < len(obs) and ref+obs[len(ref):] or ref
	
	if len(ref) != len(obs):
		print >> sys.stderr, 'ERROR!', len(ref), len(obs)
		print >> sys.stderr, ref
		print >> sys.stderr, obs
		sys.exit()
	
	for i in range(len(ref)):
		# bestref += bestRefSeqHelper(iupac2nts(ref[i]), obs[i])
		
		bestnt = False
		nts = iupac2nts(ref[i])
		for nt in nts:
			if nt == obs[i]: 
				bestref += nt; bestnt = True; break
		# if we don't find a matching allele, choose the first (i.e. randomly)
		if not bestnt: bestref += nts[0]
	if len(bestref) != len(ref): 
		print 'ERROR: RETURNING INCORRECT REF SEQUENCE!', 'original', ref, 'new', bestref
	return bestref
# # --------------------------------------------------------

def makeRange(xmin,xmax,interval=None,n=None):
	if n == None: n = int(abs(xmax-xmin)/float(interval))
	else: 
		interval = (xmax-xmin)/float(n)
		interval += interval/float(n-1)
	return [xmin+interval*i for i in range(n)]


# dictionary keyed by ascii character returning phred q-score
# ascii 33 to 126 - full range
fastqchars = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '\{', '\|', '\}', '~']
A2Q = dict(zip(fastqchars, range(len(fastqchars))))

# string conversion to integer (potentially a long)
encodeNt = lambda x: x=='A'and'1' or x=='T'and'2' or x=='C'and'3' or x=='G'and'4' or x=='N'and'5'
encodeSeq = lambda seq: int(''.join(map(encodeNt, seq)))
decodeNt = lambda x: x=='1'and'A' or x=='2'and'T' or x=='3'and'C' or x=='4'and'G' or x=='5'and'N'
decodeSeq = lambda val: ''.join(map(decodeNt, str(val)))

encodeQ = dict(zip(fastqchars, map(str,range(10,len(fastqchars)+10)))) # every character is 2 integer digits
decodeQ = dict(zip(map(str,range(10,len(fastqchars)+10)), fastqchars))
encodeQual = lambda qual: int(''.join(map(lambda x: encodeQ[x], qual)))
pairUP = lambda qual: map(lambda x: qual[x-1]+qual[x] ,range(1,len(qual),2))
decodeQual = lambda qual: ''.join(map( lambda x: decodeQ[x], pairUP(str(qual)) ))

# getStrand = lambda flag: (int(flag) & (0x10)) and '-' or '+'

# Core process for likelihood computation
def glimmrLikelihood(readfiles, GR, outfile, expPhred=30, maxreadlen=150, tmpdir='', Lambda=.5, verbose=False, bounds={}, ploidy=2, libraryErrorFiles={}, LIPieces=20):
	"""Computes the per nucleotide likelihood for occupancy of the protein generating the observed data (readfiles)"""
	
	USECCODE = False
	try: import snipercore; USECCODE = True; pipeit(' - Sniper using snipercore.so for faster likelihood calculations.',1)
	except ImportError: pipeit(' - NB, Sniper did not find snipercore.so module.',1)
	
	
	# count total genomic positions
	flg = sum(map(lambda x: len(x['seq']), GR.values()))
	
	readmaph = {}
	for k in readfiles.keys(): 
		try: readmaph[k] = open(readfiles[k])
		except IOError: pass
	
	
	# Local methods
	# ===================================================================
	ln = lambda x: math.log(x, math.e)
	log = lambda x: math.log(x, 10)
	lowestQ = int(-10*log(1e-16))
	# convert between phred score and p-value
	Q2P = lambda e: (e == 0 and 1-1e-16) or 10**(-e/10.) # avoid 1
	P2Q = lambda p: -10*log(p)
	
	A2P = dict( [ascii, 1-Q2P(A2Q[ascii])] for ascii in A2Q.keys() ) # Probability of no error
	invA2P = dict( [ascii, Q2P(A2Q[ascii])] for ascii in A2Q.keys() ) # Probability of error
	
	# compute mismatches between 2 sequences - this must be FAST
	derror = lambda a,b: (a == b and 0) or sum(map(lambda x,y: x != y and 1 or 0, a,b))
	
	def intoPieces(x,width,nbins):
		out = []
		i = 0
		npieces = 0
		while i < len(x) and npieces < nbins-1:
			out += [x[i:i+width]]
			npieces += 1
			i += width
		out += [x[i:]]
		# fill out to nbins for consistency
		while len(out) < nbins: out += ['']
		return out
	
	
	def cfPCRduplicates(D,astr):
		global PCR_HASH_WIDTH, PCR_HASH_NMAX
		p = intoPieces(astr, PCR_HASH_WIDTH, PCR_HASH_NMAX)
		# expectation is exception, so make faster using if else
		
		if p[0] in D:
			if p[1] in D[p[0]]:
				if p[2] in D[p[0]][p[1]]:
					if p[3] in D[p[0]][p[1]][p[2]]:
						if p[4] in D[p[0]][p[1]][p[2]][p[3]]:
							return True
						else: return False
					else: return False
				else: return False
			else: return False
		else: return False
		
		# try: 
		# 	D[p[0]][p[1]][p[2]][p[3]][p[4]]
		# 	return True
		# except KeyError: return False
	
	
	def loadPCRduplicates(filename):
		# load the read ids to ignore (due to PCR duplicates)
		global PCR_HASH_WIDTH, PCR_HASH_NMAX
		D = {}
		
		def loadIntoHash(count, astr):
			p = intoPieces(astr, PCR_HASH_WIDTH, PCR_HASH_NMAX)
			try: D[p[0]][p[1]][p[2]][p[3]][p[4]] = count
			except KeyError:
				try: D[p[0]][p[1]][p[2]][p[3]] = {p[4]: count}
				except KeyError:
					try: D[p[0]][p[1]][p[2]] = {p[3]: {p[4]: count}}
					except KeyError:
						try: D[p[0]][p[1]] = {p[2]: {p[3]: {p[4]: count}}}
						except KeyError:
							D[p[0]] = {p[1]: {p[2]: {p[3]: {p[4]: count}}}}
							
			# control on accurate storage
			# if ''.join(pieces) != astr: 
			# 	print "HEY!"
			# 	print 'input ', astr, ''.join(pieces) == astr
			# 	print 'output', ''.join(pieces)
			# 	if ''.join(pieces) != astr: print pieces
			# 	print
			
			return
		
		
		try:
			lh = open(filename)
			[loadIntoHash(*line[:-1].split('\t')[:2]) for line in lh]
			lh.close()
		except IOError: pipeit('Warning: Cannot access PCR dup file %s'%(filename),1)
		return D
	
	
	def loadPCRduplicatesOLD(filename):
		# load the read ids to ignore (due to PCR duplicates)
		badIDs = {}
		try:
			lh = open(filename)
			for line in lh:
				count, badid = line[:-1].split('\t')[:2]
				badIDs[hash(badid)] = None
			lh.close()
		except IOError: pipeit('Warning: Cannot access PCR dup file %s'%(filename),1)
		return badIDs
	
	
	def loadUniqueReadsOLD(tmpfile, fromChromosomes, filterids={}, subsample=None):
		# initialize D with chrom names
		D = dict( map(lambda x: (x,{}), fromChromosomes) )
		
		if not os.access(tmpfile, os.F_OK): 
			pipeit('Warning: Cannot access %s'%(tmpfile))
			return D
		fh = open(tmpfile)
		count = 1; nbytes = 0
		for line in fh:
			if count % 1000000 == 0: sys.stdout.write('%s,'%(count)); sys.stdout.flush()
			count+=1; ll = len(line)
			
			if subsample and random.random() > subsample: nbytes += ll; continue
			
			try:
				lst = line[:-1].split('\t')
				chrom = lst[2]; hhid = hash(lst[0])
				if chrom == '*' or (chrom not in D) or hhid in filterids: 
					# print 'FILTERING UNI', hhid
					nbytes += ll; continue
				
				# if chrom == '*' or (chrom not in D): nbytes += ll; continue
				# else: 
				# 	sign = str(hhid)[0] == '-' and 1 or 0 # 0 is positive
				# 	p1 = sign==1 and int(str(hhid)[1]) or int(str(hhid)[0])
				# 	p2 = sign==1 and int(str(hhid)[2]) or int(str(hhid)[1])
				# 	if hhid in filterids[sign][p1][p2]: 
				# 		print 'FILTERING UNI', hhid
				# 		nbytes += ll; continue
				# 	else: print "OK", hhid, count
				
				pos = int(lst[3])-1
				
				# store bytes in file to this line
				if pos in D[chrom]: D[chrom][pos] += [(nbytes,ll)]
				else: D[chrom][pos] = [(nbytes,ll)]
				
			# exception if a line beginning with @
			except IndexError: pass
			except ValueError: pass
			nbytes += ll
		fh.close()
		return D
	
	
	def loadDegReads(readfiles, chromosomesOfInterest, filterids={}, subsample=None):
		"""Load the degenerate portion of the read map for any reads aligning to this chromosome."""
		
		# convert to dict
		chint = dict( map(lambda x: (x,None), chromosomesOfInterest) )
		
		# gather reads of interest wrt this chrom group
		# Since read map is organized by read, really we just need to save
		# information for the first alignment along with line number in the file and number of alignments
		RdxChr = {'PE':{}, 'SE':{}}
		# chrKeep = {}
		for chrom in GR.keys(): # needs all chrom keys
			# chrKeep[chrom] = None
			RdxChr['PE'][chrom] = {}
			RdxChr['SE'][chrom] = {}
		
		# now collect all reads with read names of interest for this portion of genome
		for pese in ['PE', 'SE']:
			degfile = readfiles[pese+'deg']
			if os.access(degfile, os.F_OK):
				# temporary storage to summarize info for all deg alignments for a given hitid
				curHitid = None; curPese = None; curInfo = []; curByte = []
				# curByteStart = 0; curNBytes = 0
				count = 1; nbytes = 0 # update line by line
				fh = open(degfile)
				for line in fh:
					# if count == 10000: break
					if count % 1000000 == 0: pipeit('%.1e,'%(count))
					count+=1; ll = len(line) # each character is 1 byte in python
					
					# allow coverage subsampling
					if subsample and random.random() > subsample: nbytes += ll; continue
					
					try:
						lst = line[:-1].split('\t')
						chrom = lst[2]; hitid = lst[0]
						
						if chrom == '*': 
							# if hitid == curHitid: curNBytes += ll
							nbytes += ll; continue
						
						pos=int(lst[3])-1; seq=lst[9]; qual=lst[10]; rl = len(seq)
						
						# must record all instances of multiply-mapped reads (relevant to bounds)
						# but for each entry store first line for a new readid and number of aligments
						# this allows us to avoid the otherhits dictionary
						if curHitid == None:
							# So this is the first read we encounter in the file
							curHitid = hitid # name of the multiply-mapped read
							curPese = pese # store in PE or SE part of dictionary
							curInfo = [(chrom,pos,rl)] # minimal info for each alignment
							# curByteStart = nbytes; curNBytes = ll # pointer to block of lines in read map
							curByte = [(nbytes,ll)]
						elif hitid == curHitid:
							# continuing with multiple alignments of the same read
							curInfo += [(chrom,pos,rl)] # add this alignment
							# curNBytes += ll # the total number of bytes to read
							curByte += [(nbytes,ll)]
						else:
							# then we moved on to a new read. Save previous information
							# to the non-unique read pileup... if it is relevant
							
							# print '=>', fetchLines(open(degfile), curByte)
							# print
							
							# Currently, only want alignments for reads that occur on chromosomes of interest
							# so if none of these alignments meets this criterion, toss them all
							meetsCriterion = False
							# if curHitid not in filterids:
							if not cfPCRduplicates(filterids, curHitid):
								for cc,cp,clen in curInfo:
									if cc in chint: meetsCriterion = True; break
							
							if meetsCriterion:
								# to recover all of the alignments for this read
								# we only need the starting byte number and total bytes to read
								# stuff = (curByteStart,curNBytes) # create object here so all entries point to it
								
								# only store first maxnumalign alignments
								for cc,cp,clen in curInfo: # store each alignment, pointing them to the block of bytes
									if cp in RdxChr[curPese][cc]: RdxChr[curPese][cc][cp] += [curByte[:maxnumalign]]
									else: RdxChr[curPese][cc][cp] = [curByte[:maxnumalign]]
							
							# update
							curHitid = hitid
							curPese = pese
							curInfo = [(chrom,pos,rl)]
							curByte = [(nbytes,ll)]
					
					# exceptions if line begins with @
					except IndexError: pass
					except ValueError: pass
					nbytes += ll
				fh.close()
			else:
				pipeit('Warning: Cannot access %s'%(degfile))
		return RdxChr
	
	
	def loadDegIndex(idxhandle, chromosomesOfInterest, bounds={}, subsample=None, maxaln=1e10):
		"""Load the degenerate portion of the read map for any reads aligning to this chromosome."""
		
		# convert to dict
		chint = dict( map(lambda x: (x,None), chromosomesOfInterest) )
		
		Bb = {} # set up dictionary based on bounds
		for k in bounds.keys():
			for low,high,flag in filter(lambda x: nastr not in x, sorted(bounds[k])):
				if flag == '+':
					try: Bb[k] += [(low,high)]
					except KeyError: Bb[k] = [(low,high)]
					
		# gather reads of interest wrt this chrom group
		# Since read map is organized by read, really we just need to save
		# information for the first alignment along with line number in the file and number of alignments
		RdxChr = {}
		for chrom in chromosomesOfInterest: RdxChr[chrom] = {}
		
		# now collect all reads with read names of interest for this portion of genome
		count = 1
		try:
			line = idxhandle.next()
			while line:
				count % 1000000 == 0 and pipeit('%.1e,'%(count))
				count +=1 
				
				chrom,pos,startBytes,nlines = line[:-1].split('\t')
				pos = int(pos); startBytes = int(startBytes); nlines = int(nlines)
				if chrom in chint:
					# if len(Bb) and not poverlaps(Bb[chrom], pos): continue
					if subsample and random.random() > subsample: line = idxhandle.next(); continue
					
					# load up to specified number of alignments
					parsedAlns = [[startBytes,nlines>maxaln and maxaln or nlines]]
					
					# only store the first maxnumalign alignments
					# trick is that we don't need to store this at all locations...but pileup below should handle it
					try: RdxChr[chrom][pos] += [parsedAlns]
					except KeyError: RdxChr[chrom][pos] = [parsedAlns]
				else: # file is sorted so just return
					return RdxChr
				line = idxhandle.next()
		except StopIteration: return RdxChr
		return RdxChr
	
	
	
	def reloadBuffer(target, handle, limit):
		try:
			for i in xrange(limit): target += [handle.next()]
		except StopIteration: return -1
		return 1
	
	
	
	# Variable initialization
	# ===================================================================
	nucleotides = ['A', 'T', 'C', 'G']
	
	genotypes = None; genoUni = None
	if ploidy == 2:
		genotypes = [('A','A'),('A','T'),('A','C'),('A','G'),('T','T'),('T','C'),('T','G'),('C','C'),('C','G'),('G','G')]
		genoUni = {('A','A'):'AA',('A','T'):'AT',('A','C'):'AC',('A','G'):'AG',
				   ('T','T'):'TT',('T','C'):'CT',('T','G'):'GT',('C','C'):'CC',
				   ('C','G'):'CG',('G','G'):'GG'}
	else:
		genotypes = ['A', 'T', 'C', 'G']
		genoUni = {'A':'A', 'T':'T', 'C':'C', 'G':'G'}
		
	rangelengenotypes = range(len(genotypes))
	GT = genotypes # alias
	
	invGT = None
	if ploidy == 2:
		invGT = dict( map(lambda i: ((GT[i][0],GT[i][1]),i) , rangelengenotypes) )
		invGT2 = dict( map(lambda i: ((GT[i][1],GT[i][0]),i) , rangelengenotypes) )
		invGT.update(invGT2)
	else:
		invGT = dict( map(lambda i: (GT[i],i) , rangelengenotypes) )
	
	states = None
	stateUni = None
	if ploidy == 2:
		# FUTURE WORK: allele-specific protein binding
		# states = [('N','N'), ('N','Y'), ('Y','Y')]
		# stateUni = {('N','N'):'NN', ('N','Y'):'YN', ('Y','Y'):'YY'}
		
		# states = [('N','N'), ('Y','Y')]
		# stateUni = {('N','N'):'NN', ('Y','Y'):'YY'}
		
		states = [('Y','Y')]
		stateUni = {('Y','Y'):'YY'}
		
	else:
		states = ['Y']
		stateUni = {'Y':'Y'}
	
	rangelenstates = range(len(states))
	Ystate = states[0]
	# ===================================================================
	
	
	# Precompute nt probability scores
	# =========================================================================
	
	# PriS = pr. sequencing read ri from template sequence S
	# --------------------------------------------------------------
	# generalize this for all read lengths - 0 to maxreadlen
	binCoef = [[nchoosek(rl,d) for d in range(rl+1)] for rl in range(maxreadlen+1)]
	# call: binomialP(readlength)(mismatches, Perror)
	binomialP = lambda rl: lambda d,p: binCoef[rl][d] * p**d * (1-p)**(rl-d)
	
	# tune between read/position-specific error vs a global/position specific error
	expError = 10**(-expPhred/10.)
	PriSbin = [[Lambda * binomialP(rl)(d, expError) for d in range(rl+1)] for rl in range(maxreadlen+1)]
	
	# index by [rl]
	offsetRange = map(lambda x: range(x), xrange(maxreadlen+1))
	opLambda = 1. - Lambda
	statePris = None
	
	# binomial + quality model
	if Lambda >=0 and Lambda < 1.0:
		print ' - Using tuned binomial+quality read likelihood model for P(ri | s, G) with lambda=%s'%(Lambda)
		
		def statePriS(S, R, Q, states, rl):
			# localize variables for speed
			oR = offsetRange
			Pb = PriSbin
			lA2P = A2P
			linvA2P = invA2P
			
			res = {}
			# avoid the offset position
			for state in states:
				d = 0 # number of mismatches in the alignment
				p = 1 # read/position specific likelihood
				if state == 'Y':# or state == ('Y','Y'):
					for i in oR[rl]:
						p *= (S[i] == R[i] and lA2P[Q[i]]) or linvA2P[Q[i]]
						d += int(S[i] != R[i])
				elif state == 'N':# or state == ('N','N'):
					for i in oR[rl]:
						p *= (S[i] == R[i] and linvA2P[Q[i]]) or lA2P[Q[i]]
						d += int(S[i] == R[i])
				
				try: res[state] = Pb[rl][d] + opLambda*p
				except IndexError: print 'error', state, rl, d, Pb[rl]
			return res
		
		
	elif Lambda == 1.0:
		# binomial only model when lambda = 1
		print ' - Using binomial only read likelihood model for P(ri | s, G)'
		
		def statePriS(S, R, Q, states, rl):
			# localize variables for speed
			oR = offsetRange
			Pb = PriSbin
			
			res = {}
			for state in states:
				res[state] = (state=='Y' and Pb[rl][sum(map(lambda i: int(S[i] != R[i]), oR[rl]))]) or (state=='N' and Pb[rl][sum(map(lambda i: int(S[i] == R[i]), oR[rl]))])
				
				# d = 0 # number of mismatches in the alignment
				# if state == 'Y':# or state == ('Y','Y'):
				# 	d = sum(map(lambda i: int(S[i] != R[i]), oR[rl]))
				# 	# for i in oR[rl]: d += int(S[i] != R[i])
				# elif state == 'N':# or state == ('N','N'):
				# 	for i in oR[rl]: d += int(S[i] == R[i])
				# res[state] = PriSbin[rl][d]
			return res
	elif Lambda < 0:
		# default to a constant model of read alignment probability
		print ' - Using constant read likelihood model for P(ri | s, G)'
		
		def statePriS(S, R, Q, states, rl):
			# localize variables for speed
			oR = offsetRange
			Pb = PriSbin
			
			res = {}
			val = 0.99
			for state in states:
				res[state] = (state=='Y' and val) or (state=='N' and val)
			return res
		
	# DONE CHOOSING METHOD
	
	# initialize factorial precomputations
	if USECCODE: crap = snipercore.initprefact()
	
	# # switches between C code and Python code
	# compLik = lambda S,R,Q,O,RL,L,E: USECCODE and snipercore.ntPriS(S,R,Q,O,RL,L,E) or ntPriS(S,R,Q,O,RL)
	# 
	# # DONE CHOOSING METHOD
	# if USECCODE: crap = snipercore.initprefact()
	
	
	# PriSj = pr. that Sj was the template out of all SGri templates
	# --------------------------------------------------------------
	maximumcoverage = 100001 # this will crash if any locus has greater read depth
	PriSj = [nastr]+map(lambda n: 1/float(n), range(1,maximumcoverage)) # uniform, precomputed
	
	# =========================================================================
	
	
	# Set up bounds for subset of genome to process
	# =========================================================================
	GRNA = {}
	if len(bounds.keys()):
		flg = 0
		print ' - We have genotyping bounds...'
		
		null = [nastr, 'NA']
		for tmpk in GR.keys():
			master = GR[tmpk]['seq'][:]
			
			SEEN = False
			for k in [tmpk, tmpk.lower(), tmpk.upper()]:
				if k in bounds and not SEEN:
					SEEN = True
					lastpos = 0
					# may be many positions of interest
					for low,high,flag in sorted(bounds[k]):
						if flag == '+':
							# remove positions not specified in a range
							
							if low not in null and high not in null:
								# we have a range
								tmps = master
								tmp = ''.join(map(lambda x: 'N', tmps[0:low]))
								tmp += tmps[low:high+1]
								tmp += ''.join(map(lambda x: 'N', tmps[high+1:len(tmps)]))
								master = tmp
							elif low not in null:
								# we have a single position of interest
								# since we've sorted these positions, 
							
								# erase between lastpos and low
								tmp = master[0:lastpos+1]
								tmp += ''.join(map(lambda x: 'N', master[lastpos+1:low]))
								tmp += master[low]
								tmp += master[low+1:len(master)]
								master = tmp
								lastpos = low
							# else just leave the entire chromosome intact
						elif flag == '-':
							# set ranges to NA
							if low not in null and high not in null:
								tmps = master
								tmp = tmps[0:low]
								tmp += ''.join(map(lambda x: 'N', tmps[low:high+1]))
								tmp += tmps[high+1:len(tmps)]
								master = tmp
							
					# now I can finish erasing the end
					if lastpos:
						master = master[0:lastpos] + ''.join(map(lambda x: 'N', master[lastpos:len(master)]))
					
					# write back
					GRNA[tmpk] = {'id':tmpk, 'seq':master}
					flg += len(master)
	else: GRNA = GR
	
	# Initialize variables for likelihood computation
	# =========================================================================
	totalLoci = 0     # count number of reference positions passed
	numHypotheses = 0 # count up number of hypotheses we make
	
	outfh = open(outfile, 'w') # mother file
	
	# blank template list for likelihood calculations below (copy from this)
	LLg0 = [0 for i in genotypes]
	
	# now proceed to genotype each reference genomic position
	
	chromosomesOfInterest = GRNA.keys()
	chromosomesOfInterest.sort()
	nchroms = len(chromosomesOfInterest)
	
	# split chromosomes into LIPieces loading groups
	chromgroups = []
	copydir = chromosomesOfInterest[:]
	if LIPieces > nchroms: LIPieces = max(1,nchroms)
	itemspergroup = math.ceil(nchroms / float(LIPieces))
	
	while len(copydir):
		currgrp = []
		while len(currgrp) < itemspergroup and len(copydir):
			currgrp += [copydir.pop(0)]
		chromgroups += [currgrp]
	
	globRu = {'PE':{}, 'SE':{}} # unique pileup
	globRd = {} # degenerate pileup
	
	# Open pointer to sorted unique files
	# reopen for each chromosome because sam sorted order might differ from scan order
	readhandle = {}
	if len(readfiles['SEuni']) and 'SEuni' in readmaph: readhandle['SE'] = readmaph['SEuni']
	if len(readfiles['PEuni']) and 'PEuni' in readmaph: readhandle['PE'] = readmaph['PEuni']
	RHK = readhandle.keys()
	
	# Load the index for this file
	mapindex = {}
	for pese in RHK:
		mapindex[pese] = dict(readTable(readfiles[pese+'uni']+'.index', header=0, rownames=0)[0])
		for k in mapindex[pese].keys(): mapindex[pese][k] = int(mapindex[pese][k])
	
	
	# Open pointer to sorted deg index
	dreadhandle = {}
	if PESECOMBO and len(readfiles['SEdeg']) and 'SEdeg' in readmaph:
		dreadhandle['SE'] = open(btmapdir+sampleID+'.degenerate.map.index')
	else:
		if len(readfiles['SEdeg']) and 'SEdeg' in readmaph: 
			dreadhandle['SE'] = open(btmapdir+sampleID+'.SE.degenerate.map.index')
		if len(readfiles['PEdeg']) and 'PEdeg' in readmaph: 
			dreadhandle['PE'] = open(btmapdir+sampleID+'.PE.degenerate.map.index')
	dRHK = dreadhandle.keys()
	# if unique mapping, then set this to null
	if maxnumalign == 1: dRHK = []
	
	USEDEGINDEX = False; degmapindex = {}
	if os.access(readfiles['SEdeg']+'.metaindex', os.F_OK) or os.access(readfiles['PEdeg']+'.metaindex', os.F_OK):
		USEDEGINDEX = True
		# Load the deg index for this file
		pipeit(' - Loading degenerate map index...')
		for pese in dRHK:
			if not os.access(readfiles[pese+'deg']+'.metaindex', os.F_OK): continue
			degmapindex[pese] = dict(readTable(readfiles[pese+'deg']+'.metaindex', header=0, rownames=0)[0])
			for k in degmapindex[pese].keys(): degmapindex[pese][k] = int(degmapindex[pese][k])
		pipeit('Done.',1)
	
	lipct = 0
	for chromosomesOfInterest in chromgroups:
		lipct += 1
		
		# load the relevant reads for this chromosome group
		# -------------------------------------------------
		PESEs = ['SE', 'PE']
		if not USESINGLETONS: PESEs = ['PE']
		
		# # keep unique PCR duplciates loaded in memory
		# # unfortunately reads don't have chrom info so can't just load portion of interest
		# pipeit(' - Loading unique PCR duplicate reads...')
		# filterids = loadPCRduplicates(libraryErrorFiles['unique'])
		# pipeit('Done.',1)
		
		# for pese in PESEs:
		# 	pipeit(' - Indexing unique %s reads | loading group %s/%s... '%(pese, lipct, LIPieces))
		# 	globRu[pese] = loadUniqueReads(readfiles[pese+'uni'], chromosomesOfInterest, filterids, subsampleFold)
		# 	pipeit('Done.',1)
		
		filterids = {}
		
		globRd = None
		if not USEDEGINDEX:
			# index degenerate reads
			badIDs = {}
			if len(libraryErrorFiles['degenerate']):
				pipeit(' - Loading degenerate PCR duplicate reads...')
				badIDs = loadPCRduplicates(libraryErrorFiles['degenerate'])
				pipeit('Done.',1)
			pipeit(' - Indexing degenerate map | loading group %s/%s...'%(lipct, LIPieces))
			globRd = loadDegReads(readfiles, chromosomesOfInterest, badIDs, subsampleFold)
			pipeit('Done.',1)
			del badIDs # delete the bad reads to save memory
		
		
		# keep unique PCR duplciates loaded in memory
		# unfortunately reads don't have chrom info so can't just load portion of interest
		if len(libraryErrorFiles['unique']):
			pipeit(' - Loading unique PCR duplicate reads...')
			filterids = loadPCRduplicates(libraryErrorFiles['unique'])
			pipeit('Done.',1)
		
		for chrom in chromosomesOfInterest:
			# reset variables per chromosome
			GRseq = GR[chrom]['seq'] # want actual sequence for proper derror calculations
			chromlen = len(GRseq)
			
			# print status update every 10% per chromosome/scaffold
			printcutoff = chromlen > 10000 and int(chromlen/10.) or 10000
			
			# clear read pileup and marked list every new chromosome
			SGr.clear() # stores read pileup - this is a global variable
			prevGi = -maxreadlen # remember the last position Gi we visited
			
			REACHEDENDOFFILE = {'SE':False, 'PE':False}
			dREACHEDENDOFFILE = {'SE':False, 'PE':False}
			
			# set the file pointer at first row for this chromosome
			for pese in RHK:
				try: readhandle[pese].seek(mapindex[pese][chrom])
				except KeyError: REACHEDENDOFFILE[pese] = True
			
			for pese in dRHK:
				try: dreadhandle[pese].seek(degmapindex[pese][chrom])
				except KeyError: dREACHEDENDOFFILE[pese] = True
			
			# if we have a deg index file, we don't need to worry about LIPieces, just load what we need per chromosome
			# using the meta index, load the relevant index for the chromosome group
			if USEDEGINDEX and not DEGREADSTREAM:
				pipeit(' - Loading degenerate index for %s...'%(chrom))
				globRd = {}
				for pese in dRHK:
					globRd[pese] = loadDegIndex(dreadhandle[pese], [chrom], bounds=bounds, maxaln=maxnumalign)
				pipeit('done.',1)
			
			
			# save the previous line read
			curLine = {'SE':[0,0,0,0], 'PE':[0,0,0,0]}
			dcurLine = {'SE':[0,0,0,0], 'PE':[0,0,0,0]}
			
			pipeit(' - Iterating nucleotide positions...',1)
			
			# ITERATE OVER REFERENCE GENOMIC POSITIONS
			# ----------------------------------------
			for Gi in xrange(chromlen):
				if GRNA[chrom]['seq'][Gi] == 'N': continue
				
				ref = GRseq[Gi]
				totalLoci += 1
				
				Gi % printcutoff == 0 and \
					pipeit('   %s| %s %s/%s = %s (Hypos=%s | %.2f'%(sampleID, chrom, Gi,chromlen, ref, numHypotheses, 100*(totalLoci/float(flg)))+'% complete)',1)
				
				# GET READ PILEUP FOR THIS POSITION
				# ---------------------------------------------------------------
				# add entries from reads within some range of positions
				
				# potentially improve this for larger inserts - just avoid the middle
				# for pos in xrange(max(0, Gi-fragRad+1, prevGi+1), min(Gi+fragRad+1, chromlen)):
				# print Gi, len(xrange(max(0, Gi-fragRad+1), min(Gi+fragRad+1, chromlen))), len(SGr)
				
				# for pos in xrange(max(0, Gi-fragRad+1), min(Gi+fragRad+1, chromlen)):
				
				# slight code duplication this way but fewer if statements
				if DEGREADSTREAM:
					for pos in xrange(max(0, Gi-fragRad+1, prevGi+1), min(Gi+fragRad+1, chromlen)):
						while len(SGr) > maxdepth: del SGr[random.choice(SGr.keys())]
						
						# having group-indexed the reads by starting position we can do this
						# each key of SGr is the read id, returning a list of pos,seq,phred triples
						
						# FETCH CORRESPONDING UNIQUE READS
						# --------------------------------
						for pese in RHK:
							# having group-indexed the reads by starting position we can do this
							# each key of SGr is the read id, returning a list of pos,seq,phred triples
							
							# first check previous line
							hlo = curLine[pese][0]
							if curLine[pese][2]==chrom and int(curLine[pese][3])-1 == pos and not cfPCRduplicates(filterids, hlo):
								if subsampleFold and random.random() > subsampleFold: continue
								strand = int(curLine[pese][1])&(0x10) and '-' or '+'
								#                         pos         seq    qual     readlength, strand
								SGr[(hlo,pese)] = [(int(curLine[pese][3])-1, curLine[pese][9], curLine[pese][10], len(curLine[pese][9]), strand)]
							
							if not REACHEDENDOFFILE[pese]:
								try:
									while curLine[pese][2] < chrom or (curLine[pese][2]==chrom and int(curLine[pese][3])-1 <= pos):
										try:
											curLine[pese] = readhandle[pese].next()[:-1].split('\t')
											# check chrom, position, and filter file
											hlo = curLine[pese][0]
											if curLine[pese][2]==chrom and int(curLine[pese][3])-1 == pos and not cfPCRduplicates(filterids,hlo):
												if subsampleFold and random.random() > subsampleFold: continue
												strand = int(curLine[pese][1])&(0x10) and '-' or '+'
												SGr[(hlo,pese)] = [(int(curLine[pese][3])-1, curLine[pese][9], curLine[pese][10], len(curLine[pese][9]), strand)]
										except IndexError: 
											# reload buffer
											retval = reloadBuffer(URB[pese], readhandle[pese], bufferlimit)
											if retval == -1: break # end of file
								except StopIteration: REACHEDENDOFFILE[pese] = True
						
						# FETCH CORRESPONDING DEGENERATE READS
						# ------------------------------------
						# stream-read the deg index files
						# after we fetch a line, then we need to fetch the sequence and quality information from the actual SAM file
						for pese in dRHK:
							# first check previous line
							if dcurLine[pese][0]==chrom and int(dcurLine[pese][1]) == pos:
								# parse this line: chrom,pos,alns
								# glob = [map(lambda x: map(int,x.split(':')), dcurLine[pese][2].split(','))]
								
								startBytes = int(dcurLine[pese][2])
								nlines = int(dcurLine[pese][3])
								glob = [[[startBytes,nlines>maxnumalign and maxnumalign or nlines]]]
								
								# now fetch information from SAM file
								for lineblock in map(lambda bytelst: fetchLines(readmaph['%sdeg'%(pese)],bytelst), glob):
									altaln = []; ehitid = None
									for line in lineblock:
										if subsampleFold and random.random() > subsampleFold: continue
										
										lst = line[:-1].split('\t')
										ohitid = lst[0]; ochrom = lst[2]; opos=int(lst[3])-1; oseq=lst[9]; oqual=lst[10]
										
										if opos == pos:
											ehitid = ohitid
											estrand = (int(lst[1])&(0x10) and '-' or '+')
											SGr[(ehitid,pese)] = [(opos, oseq, oqual, len(oseq), estrand)]
										else: altaln += [(opos, ochrom, oseq, oqual)]
									try: SGr[(ehitid,pese)] += altaln
									except KeyError: pass # left over - there is no alignment that overlaps this position
							
							if not dREACHEDENDOFFILE[pese]:
								try:
									while dcurLine[pese][0] < chrom or (dcurLine[pese][0]==chrom and int(dcurLine[pese][1]) <= pos):
										try:
											dcurLine[pese] = DRB[pese].pop(0)[:-1].split('\t')
											# dcurLine[pese] = dreadhandle[pese].next()[:-1].split('\t')
											
											if dcurLine[pese][0]==chrom and int(dcurLine[pese][1]) == pos:
												# parse this line: chrom,pos,alns
												# glob = [map(lambda x: map(int,x.split(':')), dcurLine[pese][2].split(','))]
												
												startBytes = int(dcurLine[pese][2])
												nlines = int(dcurLine[pese][3])
												glob = [[[startBytes,nlines>maxnumalign and maxnumalign or nlines]]]
												
												# now fetch information from SAM file
												for lineblock in map(lambda bytelst: fetchLines(readmaph['%sdeg'%(pese)],bytelst), glob):
													altaln = []; ehitid = None
													for line in lineblock:
														if subsampleFold and random.random() > subsampleFold: continue
														
														lst = line[:-1].split('\t')
														ohitid = lst[0]; ochrom = lst[2]; opos=int(lst[3])-1; oseq=lst[9]; oqual=lst[10]
														
														if opos == pos:
															ehitid = ohitid
															estrand = (int(lst[1])&(0x10) and '-' or '+')
															SGr[(ehitid,pese)] = [(opos, oseq, oqual, len(oseq), estrand)]
														else: altaln += [(opos, ochrom, oseq, oqual)]
													try: SGr[(ehitid,pese)] += altaln
													except KeyError: pass # left over - there is no alignment that overlaps this position
										except IndexError: 
											# reload buffer
											retval = reloadBuffer(DRB[pese], dreadhandle[pese], bufferlimit)
											if retval == -1: break # end of file
								except StopIteration: dREACHEDENDOFFILE[pese] = True
					
					
					
				else:
					for pos in xrange(max(0, Gi-fragRad+1, prevGi+1), min(Gi+fragRad+1, chromlen)):
						while len(SGr) > maxdepth: del SGr[random.choice(SGr.keys())]
						
						# having group-indexed the reads by starting position we can do this
						# each key of SGr is the read id, returning a list of pos,seq,phred triples
						
						# FETCH CORRESPONDING UNIQUE READS
						# --------------------------------
						for pese in RHK:
							# having group-indexed the reads by starting position we can do this
							# each key of SGr is the read id, returning a list of pos,seq,phred triples
							
							# first check previous line
							hlo = curLine[pese][0]
							if curLine[pese][2]==chrom and int(curLine[pese][3])-1 == pos and not cfPCRduplicates(filterids, hlo):
								if subsampleFold and random.random() > subsampleFold: continue
								strand = int(curLine[pese][1])&(0x10) and '-' or '+'
								#                         pos         seq    qual     readlength, strand
								SGr[(hlo,pese)] = [(int(curLine[pese][3])-1, curLine[pese][9], curLine[pese][10], len(curLine[pese][9]), strand)]
							
							if not REACHEDENDOFFILE[pese]:
								try:
									while curLine[pese][2] < chrom or (curLine[pese][2]==chrom and int(curLine[pese][3])-1 <= pos):
										try:
											curLine[pese] = readhandle[pese].next()[:-1].split('\t')
											# check chrom, position, and filter file
											hlo = curLine[pese][0]
											if curLine[pese][2]==chrom and int(curLine[pese][3])-1 == pos and not cfPCRduplicates(filterids,hlo):
												if subsampleFold and random.random() > subsampleFold: continue
												strand = int(curLine[pese][1])&(0x10) and '-' or '+'
												SGr[(hlo,pese)] = [(int(curLine[pese][3])-1, curLine[pese][9], curLine[pese][10], len(curLine[pese][9]), strand)]
										except IndexError: 
											# reload buffer
											retval = reloadBuffer(URB[pese], readhandle[pese], bufferlimit)
											if retval == -1: break # end of file
								except StopIteration: REACHEDENDOFFILE[pese] = True
						
						# FETCH CORRESPONDING DEGENERATE READS
						# ------------------------------------
						# non streaming version
						for pese in PESEs:
							if pos not in globRd[pese][chrom] or not readmaph['%sdeg'%(pese)]: continue
							
							for lineblock in map(lambda bytelst: fetchLines(readmaph['%sdeg'%(pese)],bytelst), globRd[pese][chrom][pos]):
								altaln = []; ehitid = None; estrand = None
								for line in lineblock:
									if subsampleFold and random.random() > subsampleFold: continue
									lst = line[:-1].split('\t')
									ohitid = lst[0]; ochrom = lst[2]; opos=int(lst[3])-1; oseq=lst[9]; oqual=lst[10]
									if opos == pos:
										ehitid = ohitid; estrand = (int(lst[1])&(0x10) and '-' or '+')
										SGr[(ehitid,pese)] = [(opos, oseq, oqual, len(oseq), estrand)]
									else: altaln += [(opos, ochrom, oseq, oqual)]
								try: SGr[(ehitid,pese)] += altaln
								except KeyError: pass # left over - there is no alignment that overlaps this position
						
				# DONE OBTAINING READ PILEUP
				prevGi = Gi+fragRad # update previous position we visited
				
				# figure out why only portion of data set is getting processed...
				# if len(SGr) > 0: numHypotheses += 1
				# continue
				
				####################################################################
				
				# ASSESS THE GENOTYPE AT THIS POSITION
				# ------------------------------------
				# Compute product likelihood of read pileup across genotypes
				LL = 0 #LLg0[:]                             # read-product likelihoods
				multiplicity = 0 # number of reads that are degenerate (have >1 alignment; some k < n)
				totalmultiplicity = [] # total number of degenerate alignments (d-n)
				# mmpr = [] # number of alignment mismatches per read
				marked = [] # list of readids marked to be deleted
				
				# what is final read depth at position Gi in resequenced genome?
				# this depends on the diameter - so number of reads used to compute
				# a probability at this position
				depth = 0
				
				for readidpese,alignments in SGr.iteritems():
					(pos,seq,qual,curRL,dirdir) = alignments[0]
					
					# ignore mis-directed reads
					if dirdir != (pos <= Gi and '+' or '-'): continue
					# ignore lame-duck reads
					if pos <= Gi-fragRad: marked += [readidpese]; continue
					
					depth += 1
					altAlignments = alignments[1:]
					
					# if len(altAlignments):
					# 	print >> sys.stderr, 'first', pos, seq, curRL
					# 	print >> sys.stderr, 'remain', altAlignments
					
					# what is the size of SG(l), the k-bounded string set 
					# for the resequenced genome SG at position l
					# that is, total number of degenerate alignments + alignment at this position
					# This is d+1 = |alignments|
					
					# However, it is important to note that PE reads no longer carry /1 and /2,
					# so we need to distinguish different ends by quality
					# each mate will have one of 2 qualities (reverse)
					# this is actually an assumption in sniper - we break up the degenerate paired end reads
					if len(altAlignments) > 0 and readidpese[1] == 'PE':
						# that is, degenerate paired end reads
						# we only want the group of alignments carrying the alignment of interest
						groups = {}
						for tup in altAlignments:
							cq = tup[-1]; rcq = cq[::-1]
							if cq in groups: groups[cq] += [tup]
							elif rcq in groups: groups[rcq] += [tup]
							else: groups[cq] = [tup]
						
						# now which group do we want here?
						if qual in groups: altAlignments = groups[qual]
						else: 
							try: altAlignments = groups[qual[::-1]]
							except KeyError: altAlignments = [] # no deg alignments of interest for this read
					
					# final update to number of templates
					numTemplates = len(altAlignments)+1 # = d+1, where d is multiplicity for this read
					
					# now have d alignments, which are all degenerate
					multiplicity += (numTemplates>1 and 1) or 0 # are there any degenerate alignments?
					totalmultiplicity += [numTemplates-1] # add d = len(alignments)
					
					#    likelihood = (1/d+1)P(ri|Gl=x) + (1/d+1)Pdeg1 + ... + (1/d+1)Pdegd
					#               = (1/d+1)[P(ri|Gl=x) + Pdeg1 + ... + Pdegd]
					#               = (1/d+1)[P(ri|Gl=x) + sum_d(Pdeg)]
					#               ~ (1/d+1)P(ri|Gl=x) + (1/d+1)d x Pdeg
					#    so we are just computing a weighted average
					
					# now the 1/d+1 term is the probability of one genomic segment 
					# serving as the template for sequencing: P(s | G)
					# here this is uniformly distributed across reads
					# below is an implementation where P(s | G) ~ P(s | ri, G)
					# where we solve in a Bayesian fashion:
					# P(s | ri, G) = P(ri | s, G)P(s,G) / Sum P(ri | s, G)P(s,G)
					# with uniform prior this reduces to
					# P(s | ri, G) = P(ri | s, G) / Sum P(ri | s, G)
					# which is a weighting based on relative strength of the read aligning to 
					# that genomic region
					# basically this says that for degenerate reads, give the most weight 
					# to the location with fewer mismatches or higher quality
					
					# if ploidy = 2, then really have 2d + 2 alignments for this read
					Pcond = ploidy == 1 and PriSj[numTemplates] or PriSj[2*numTemplates]
					# Note, we could formulate Pcond as a read-specific probability...
					
					# compute the likelihood for this read
					
					# check if previously computed
					Prl = None
					Pdeg = None
					if readidpese in PRELIK:
						Prl = PRELIK[readidpese][0] # unique
						Pdeg = PRELIK[readidpese][1] # degenerate
					else:
						# Prlo = statePriS(seq, bestRefSeq(GRseq[pos:pos+curRL], seq), qual, states, curRL)
						# Prlo = 

						# C code
						# reference sequence may have SNPs, so choose the string 
						# with lowest edit distance to observed read
						
						# test with new code on haploid ref
						# test with old code on hap and dip ref
						# is blank ref getting used? 
						# is there weird bug? should throw an index error if ref and obs don't match length
						
						# new | DO NOT USE
						# Prl = Lambda < 0 and {'Y':0.99} or snipercore.statePriS(seq, bestRefSeq(GRseq[pos:pos+curRL], seq), qual, curRL, Lambda, expError)
						
						# SNIPERCORE CODE | USE THIS
						Prl = snipercore.statePriS(seq, bestRefSeq(GRseq[pos:pos+curRL], seq), qual, curRL, Lambda, expError)
						# print >> sys.stderr, 'Prl', Prl, curRL
						
						# old | DO NOT USE
						# Prl = snipercore.statePriS(seq, GRseq[pos:pos+curRL], qual, curRL, Lambda, expError)
						
						# if Prl != Prlo:
						# 	print 'Compare'
						# 	print 'orig', Prlo
						# 	print 'new ', Prl
						
						# sum up the degenerate probability mass: info = (opos, ochrom, seq, oqual)
						# -------------------------------------------------------------
						# do not consider reference variant positions
						
						# PURE PYTHON CODE
						# Pdeg = sum( map(lambda info: statePriS(info[2],GR[info[1]]['seq'][info[0]:info[0]+readlen],info[3])[Ystate], alignments[1:]) )*Pcond
						# Pdeg = sum( map(lambda info: statePriS(info[2],GR[info[1]]['seq'][info[0]:info[0]+curRL],info[3], ['Y'], curRL)[Ystate], altAlignments) )*Pcond
						
						# new | DO NOT USE
						# Pdeg = sum( map(lambda info: Lambda < 0 and 0.99 or snipercore.statePriS(info[2], bestRefSeq(GR[info[1]]['seq'][info[0]:info[0]+len(info[2])], info[2]), info[3], curRL, Lambda, expError)[Ystate], altAlignments) )*Pcond
						
						# THIS IS FINAL SNIPERCORE CODE!!!!!
						Pdeg = sum( map(lambda info: snipercore.statePriS(info[2], bestRefSeq(GR[info[1]]['seq'][info[0]:info[0]+len(info[2])], info[2]), info[3], curRL, Lambda, expError)[Ystate], altAlignments) )*Pcond
						
						# old | DO NOT USE
						# Pdeg = sum( map(lambda info: snipercore.statePriS(info[2],GR[info[1]]['seq'][info[0]:info[0]+curRL], info[3], curRL, Lambda, expError)[Ystate], altAlignments) )*Pcond
						# except KeyError: continue
						
						PRELIK[readidpese] = (Prl,Pdeg)
						
					# LL = ploidy == 1 and map(lambda g: LL[g] + log(Pcond*Prl[states[g]] + Pdeg), rangelenstates) or \
				 	                 # map(lambda g: LL[g] + log(Pcond*(Prl[states[g][0]] + Prl[states[g][1]]) + 2*Pdeg), rangelenstates)
					
					try:
						LL = ploidy == 1 and LL + log(Pcond*Prl[Ystate] + Pdeg) or \
					 	                 LL + log(Pcond*(Prl[Ystate[0]] + Prl[Ystate[1]]) + 2*Pdeg)
					except ValueError:
						# not sure how to resolve this, so just skip the 
						continue
						
						
						# print 'ERROR:', ploidy, 'pcond', Pcond, 'Prl', Prl, 'Pdeg', Pdeg
						# print 'alt', altAlignments
						# print 'Ps', map(lambda info: snipercore.statePriS(info[2], bestRefSeq(GR[info[1]]['seq'][info[0]:info[0]+len(info[2])], info[2]), info[3], curRL, Lambda, expError)[Ystate], altAlignments)
						# print 'info[2]', altAlignments[0][2]
						# print 'bestref', bestRefSeq(GR[altAlignments[0][1]]['seq'][altAlignments[0][0]:altAlignments[0][0]+len(altAlignments[0][2])], altAlignments[0][2])
						# print 'info[3]', altAlignments[0][3]
						# print 'curRL', curRL
						# print 'parms', Lambda, expError
						# sys.exit('EXIT WITH INF ERRROR')
						
						
				# done loop over reads
				map(delete,marked) # delete marked reads from SGr
				
				# don't print out info if position has no depth
				if depth == 0: continue
				
				# convert mismatches to weights / frequencies that sum to 1
				# smpr = float(sum(mmpr)); flpr = float(len(mmpr))
				# mmpr = map(lambda x: smpr > 0 and x/smpr or 1/flpr, mmpr)
				
				# Multiplicity ratio: d-n / n (or Ndeg / N), ranging from 0 to inf
				# aka the average number of additional alignments per read / locus
				# Sum(d_read) / Sum(reads) per locus
				udalpha = round(harmonicMean(totalmultiplicity), 3)
				
				# Record statistics
				# -----------------
				info = [chrom, str(Gi), ref, str(LL), str(depth), str(multiplicity), str(sum(totalmultiplicity)), str(udalpha)]
				print >> outfh, '\t'.join(info) # if printing all base calls
				numHypotheses += 1 # we are evaluating a hypothesis
			# DONE LOOP OVER POSITIONS
		# DONE LOOP OVER CHROMOSOMES
	# DONE GENOME
	outfh.close()
	
	pipeit('   Done processing %s (n=%s hypotheses, m=%s total loci)'%(sampleID, numHypotheses, totalLoci),1)
	
	for k in readmaph.keys(): readmaph[k].close()
	
	return numHypotheses



# END METHODS
# -----------------------------------------------------------------
# -----------------------------------------------------------------

# =========================================================
# PROGRAM INITIALIZATION
# =========================================================

# CLI VARIABLES
# -------------------------------------------------------------------
outdir = slash(os.getcwd())

btmapdir = ''; SETMAPDIR = False
findir = ''
mapfile = ''
zkey = 'Summary' # default folder name for processed samples
bfaref = ''
MAPCUT = False

PID = None # unique process id for this call

################## GLIMMRLIKELIHOOD PARAMETERS ###########################
##########################################################################

maxnumalign = None

STRATEGY = 'all'                # broad mapping strategy
                                # one of: all, unique, best, bestno [default: all]
                                # unique: genotype using uniqely aligned reads only (--bestm -m 1)
                                # best: genotype using best-guess aligned reads only (--bestk -m 1)
                                # if UNIQ = BEST = False, align using --bestk/m -m d
                                # best-noguess: exclude best reads if multiply best alignments exist

MAPWITHQUALITY = False          # quality-blind or quality-aware mapping


# UNIQ = False                    # genotype using uniqely aligned reads only (--bestm -m 1)
# BEST = False                    # genotype using best-guess aligned reads only (--bestk -m 1)
#                                 # if UNIQ = BEST = False, align using --bestk/m -m d
# BESTNO = False                  # best-noguess; exclude best reads if multiply best alignments exist

PESECOMBO = False               # save PE and SE alignments to single map file

# PART III: SNP CALLING
# ---------------------
SNPCALLING = False              # whether or not to perform genotyping
REDOSNPCALLING = False          # whether to redo SNP calling on a sampleID which has already been typed
REDOFINALSNPCALLING = False     # correct for multiple testing to return final significant snps
USESINGLETONS = True            # genotype using both PE and SE reads
boundsfile = ''                 # file containing chrom\tlow\thigh rows for restricting genotyping
singlesampleid = None           # name of single sampleID to genotype
LOWSID = None; HIGHSID = None   # if user wants to bound which samples are genotyped in a map file

LIPieces = 5                    # Break loading of unique reads into this number of pieces: min(nchroms, LIPieces)
                                # greatly reduces memory requirements (but requires rescanning the read set this many times)

maxreadlen = 150                # hardwired length of all reads
Lambda = 0.67                   # weighting factor; larger values up-weight global binomial pr.
MAXALPHA = None                 # don't call any SNPs that have alpha > 10
MINALPHA = 0.0

# defined by a windowsize parameter - 
fragDiameter = 150              # mean of the bioanalyzer trace
fragUpperBound = 25             # extra ~1SD from the bioanalyzer size distribution


# Estimate of average per-nucleotide sequencing error for global binomial likelihood
# use --ge flag to estimate
expPhred = 35.                 # from alignable unique PE 19-anc reads

ploidy = 1 # consider haploid

ALLELESPECIFIC = False

REDUNDANCYFILTER = True       # mask PCR duplicates
badReadDir = None            # defined below

DEGREADSTREAM = True         # stream deg reads

# if one channel has drastically excess coverage may want to subsample randomly to avoid excessive bias
subsampleCh = None; subsampleFold = None

maxdepth = 100                # maximum reads to evaluate a nucleotide locus (if there are more only the first n will be used)

# OTHER PARAMETERS
# ----------------
REDOALL = False               # redo everything
DRYRUN = False                # dryrun the pipeline
VERBOSE = False
PYTHONPATH = False

##########################################################################
##########################################################################


# CLI ARGS
# -------------------------------------------------------------------
help = '\nHELP for glimmrLikelihood.py\n'
nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False
# print >> sys.stderr, 'GL given', ' '.join(args)
while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('out|^o$', arg): outdir = slash(val)
	elif re.match('map|^m$', arg): mapfile = val
	elif re.match('key|^z$', arg): 
		zkey = val
		if zkey[len(zkey)-1] == '/': zkey = zkey[0:len(zkey)-1]
	elif re.match('^readmapdir$|^btmapdir$|^bmd$|^i$', arg): SETMAPDIR = True; btmapdir = val
	
	elif arg == 'combo': PESECOMBO = True; ai-=1
	
	elif re.match('^ref$|^b$', arg): bfaref = val
	elif re.match('^mcap|^d$$', arg): maxnumalign = int(val)
	elif re.match('^mapcut$|^mc$', arg): MAPCUT = True; ai-=1
	elif re.match('usesingletons|pese|sepe|^g$', arg): USESINGLETONS = True; ai-=1
	elif re.match('peonly|^pe$', arg): USESINGLETONS = False; ai-=1
	
	elif arg == 'as' or arg == 'alleles': ALLELESPECIFIC = True; ai-=1
	
	# mapping strategy selection
	elif re.match('^all$|^total$', arg): STRATEGY = 'all'; ai-=1
	elif re.match('^uniq$|^unique$', arg): 
		# only return an alignment if it is unique in the reference genome
		# given specified maximum mismatches
		STRATEGY = 'unique'; alignmentPolicy = 'bestm'; maxnumalign = 1; ai-=1
	elif re.match('^best$|^bestguess$', arg): 
		# return an alignment if it is the one with the fewest mismatches
		# given specified maximum mismatches; quasi-randomly break ties
		STRATEGY = 'best'; alignmentPolicy = 'bestk'; maxnumalign = 1; ai-=1
	elif re.match('^bestno$|^bestnoguess$', arg):
		# return an alignment if it is the one with the fewest mismatches
		# but discard read if multiple equally best alignments exist
		STRATEGY = 'bestno'; alignmentPolicy = 'bestk'; maxnumalign = 2; ai-=1
		# technically we map allowing 2 alignments and then SNP call with
		# the unique maps only
	elif re.match('^mapqual$', arg): MAPWITHQUALITY = True; ai-=1
	
	elif arg == 'stream': DEGREADSTREAM = True; ai -= 1
	elif arg == 'nostream': DEGREADSTREAM = False; ai -= 1
	elif arg == 'buffer': bufferlimit = int(val)
	
	elif re.match('^lik$', arg): SNPCALLING = True; ai-=1
	elif re.match('^relik$', arg): SNPCALLING = True; REDOSNPCALLING = True; ai-=1
	# elif re.match('^redocorrection$|^rc$', arg): SNPCALLING = True; REDOFINALSNPCALLING = True; ai-=1
	# elif re.match('^redocorrectionpost$|^rcf$', arg): SNPCALLING = True; REDOFINALSNPCALLING = True; ai-=1
	
	elif re.match('^loadinpieces$|^lip$', arg): LIPieces = int(val)
	
	elif re.match('^maxalpha$|^mxa$', arg): MAXALPHA = float(val)
	elif re.match('^minalpha$|^mna$', arg): MINALPHA = float(val)
	elif re.match('^minquality$|^Q$', arg): PvalueRange = [float(val)]
	
	elif re.match('^norepfilter$|^nrf$|^nopcr$', arg): REDUNDANCYFILTER = False; ai -=1
	elif arg == 'pcrdir': badReadDir = val
	
	elif arg == 'subsample': subsampleFold = float(val)
	
	elif re.match('^restrict|^rest$', arg): boundsfile = val
	elif re.match('^sampleid|^id$', arg): singlesampleid = val
	elif re.match('^subset|^sub$', arg): LOWSID,HIGHSID = map(int,val.split(','))
	
	elif re.match('readlen|^rl$', arg): maxreadlen = int(val)
	elif re.match('^diameter$', arg): fragDiameter = int(val)
	elif re.match('^delta$', arg): fragUpperBound = int(val)
	
	# elif re.match('usecal|^uc$', arg):  USEQUALITYCALIBRATION = True; ai-=1
	elif re.match('lambda|^l$', arg): Lambda = float(val)
	elif re.match('expphred|^ep$', arg): expPhred = float(val)
	
	elif arg == 'pid': PID = val
	elif arg == 'pythonpath': PYTHONPATH = val
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)
# -------------------------------------------------------------------

if PYTHONPATH: 
	sys.path += [PYTHONPATH]
	# print 'Adding to path', PYTHONPATH


if not maxnumalign: sys.exit("Must specify --mcap INT in glimmrHet.py")

pipeit('\n'+' '.join(args)+'\n',1)

# contingency of flags
# --------------------
if REDOALL:
	RUNBOWTIE = True
	RUNMAQ = True
	SNPCALLING = True

# if not USEQUALITYCALIBRATION: Qcalibrations = None


# declaration of directories
# --------------------------
zdir = outdir+zkey+'/'

if not SETMAPDIR: btmapdir = zdir+'read_maps_all/'
unmapdir = zdir+'unmapped_reads/'
pmadir = zdir+'alignment_mismatches/'

if STRATEGY=='unique':
	if not SETMAPDIR: btmapdir = zdir+'read_maps_uni'
	unmapdir = zdir+'unmapped_reads_uni'
	pmadir = zdir+'alignment_mismatches_uni'
if STRATEGY=='best':
	if not SETMAPDIR: btmapdir = zdir+'read_maps_best'
	unmapdir = zdir+'unmapped_reads_best'
	pmadir = zdir+'alignment_mismatches_best'
if STRATEGY=='bestno':
	if not SETMAPDIR: btmapdir = zdir+'read_maps_bestno'
	unmapdir = zdir+'unmapped_reads_bestno'
	pmadir = zdir+'alignment_mismatches_bestno'

if MAPWITHQUALITY: 
	btmapdir += 'q'
	unmapdir += 'q'
	pmadir += 'q'

zinsdir = zdir+'insertsize/'

tmpstr = 'ALL'
if STRATEGY=='unique': tmpstr = 'UNI'
elif STRATEGY=='best': tmpstr = 'BEST'
elif STRATEGY=='bestno': tmpstr = 'BESTNO'
if MAPWITHQUALITY: tmpstr += 'Q'
resdir = zdir+'results-L%s,E%s,%s,%s/'%(Lambda, expPhred, USESINGLETONS and 'PESE' or 'PE', tmpstr)
tmpdir = resdir+'tmp/'

# file suffix depending on map type
ubamap = '-%s'%(STRATEGY)
if MAPWITHQUALITY: ubamap += 'q'

completedir = zdir+'completed/'
zinsdir = zdir+'insertsize/'

if not badReadDir: badReadDir = zdir+'libraryErrors/'

# create new directories
# ----------------------
createdir(zdir)
createdir(outdir);   createdir(btmapdir)
createdir(zinsdir)
createdir(unmapdir); createdir(pmadir)

createdir(resdir)

createdir(badReadDir)#; createdir(degReadDir)


# total fragDiameter with fragUpperBound
fragDiameter += fragUpperBound
fragRad = int(fragDiameter/2.) # half of a fragment

# parse map file and group replicates
# ----------------------------------------
# parse map file and group replicates
# ----------------------------------------
fh = open(mapfile)
runs = []; mapd = []; head = []; postcombos = []
INPOSTCOMBOS = False
ct = 0
for line in fh:
	# store header
	if line[0:2] == '#!': head = line[2:-1].split('\t')
	# enter posterior part
	if line[0:2] == '#@': INPOSTCOMBOS = True; continue
	
	if line[0] != '#' and not INPOSTCOMBOS:
		tmp = line[:-1].split('\t')
		if len(tmp) > 1:
			runs += [tmp.pop(0)]
			mapd += [tmp]
		
	elif line[0] != '#' and INPOSTCOMBOS:
		tmp = line[:-1].split('\t')
		# this parses the enrichment compared to input
		if len(tmp) == 2 and tmp[0] == 'Perr':
			Perr = float(tmp[1])
		# this parses the cross antibody relative efficiency
		elif len(tmp) and tmp[0] != '':
			postcombos += [tmp]
	
	ct += 1
	
# works if you dont have the post combos section
# mapd,runs,head = readTable(mapfile,keepKey=True, header=False)

seqRuns = {} # keyed by sampleID
for i in range(len(runs)):
	run,lane,sampleID,alias,refGenome,paired,path = runs[i], mapd[i][0], mapd[i][1], mapd[i][2], mapd[i][3], int(mapd[i][4]), mapd[i][5]
	path = path.strip()
	if len(path) and path[len(path)-1] != '/': path += '/'
	try: seqRuns[sampleID][(run,lane)] = [refGenome, paired, path, nastr, alias]
	except KeyError: seqRuns[sampleID] = {(run,lane):[refGenome, paired, path, nastr, alias]}

sortseqruns = seqRuns.keys()
sortseqruns.sort()

keeps = sortseqruns[LOWSID:HIGHSID]
for sid in sortseqruns:
	if sid not in keeps:
		del seqRuns[sid]

if singlesampleid and singlesampleid in seqRuns:
	bak = seqRuns[singlesampleid]
	seqRuns = {singlesampleid:bak}

sortseqruns = seqRuns.keys()
sortseqruns.sort()


if LOWSID or HIGHSID:
	print ' - Processing samples %s'%(sortseqruns[LOWSID:HIGHSID])

# load bounds file
boundslabel = ''
seqbounds = {}
if boundsfile and os.access(boundsfile, os.F_OK):
	d,r,c = readTable(boundsfile, header=False)
	for i in range(len(d)):
		tmp = []
		try: tmp += [intna(d[i][0])]
		except IndexError: tmp += [nastr]
		except ValueError: tmp += [nastr]
		try: tmp += [intna(d[i][1])]
		except IndexError: tmp += [nastr]
		except ValueError: tmp += [nastr]
		try: tmp += [d[i][2]] # flag - + or -
		except IndexError: tmp += ['+']
		keyname = r[i].strip()
		try: seqbounds[keyname] += [tmp]
		except KeyError: seqbounds[keyname] = [tmp]
	
	parts = boundsfile.split('/')
	crap, boundslabel, crap = re.split('(.+)\..+$', parts[len(parts)-1])
	boundslabel = '.'+boundslabel
	
elif boundsfile:
	# this may also be a string of a particular bound: chr01, 500, 1000
	splt = boundsfile.split(',')
	if len(splt) == 1:
		loc = splt[0]
		seqbounds[loc] = [('NA', 'NA', '+')]
		# boundslabel = '.%s,%s,%s,%s'%(loc,'NA','NA','+')
		boundslabel = '.%s'%(loc)
	elif len(splt) == 3:
		loc, low, high = map(lambda x: x.strip(' '), splt)
		seqbounds[loc] = [(intna(low), intna(high), '+')]
		boundslabel = '.%s,%s,%s,%s'%(loc,low,high,'+')
	elif len(splt) == 4:
		loc, low, high,flag = map(lambda x: x.strip(' '), splt)
		seqbounds[loc] = [(intna(low), intna(high), flag)]
		boundslabel = '.%s,%s,%s,%s'%(loc,low,high,flag)
	else:
		sys.exit('Restriction flag not specified properly: '%(boundsfile))
	
# DONE INITIALIZATION



tmpdir = resdir+'tmp/'
createdir(tmpdir)

sampleIDs = sorted(seqRuns.keys())

# fix sampleIDs to include allele specific analysis if specified
# haploidSampleIDs = sampleIDs
if ALLELESPECIFIC:
	haploidSampleIDs = []
	for sampleID in sampleIDs:
		haploidSampleIDs += [sampleID+'.h1']
		haploidSampleIDs += [sampleID+'.h2']
	sampleIDs = haploidSampleIDs

if singlesampleid:
	if singlesampleid in sampleIDs:
		sampleIDs = [singlesampleid]
	else:
		sys.exit('ID specified by --sampleid "%s" not found in map.'%(singlesampleid))

# Only load reference if re-estimating...
NEEDTOLOADREFERENCE = False
for sampleID in sampleIDs:
	motherfile = resdir+sampleID+'.mother'   # genotype at all informative positions
	motherfile += boundslabel+'.txt'
	if not os.access(motherfile, os.F_OK) or REDOSNPCALLING: NEEDTOLOADREFERENCE = True
GR = None
flg = None
if NEEDTOLOADREFERENCE:
	sys.stdout.write(' - Loading reference genome... ')
	sys.stdout.flush()
	GR = readFasta(bfaref, TOUPPER=True)
	flg = sum(map(lambda x: len(x['seq']), GR.values()))
	sys.stdout.write('N=%s loci, M=%s chromosomes\n'%(flg,len(GR)))

theresdir = resdir
if DRYRUN: theresdir = tmpdir


for sampleID in sampleIDs:
	# set up SNP calling for PE vs PE+SE / unique reads vs unique + deg reads
	pref = btmapdir+sampleID
	rf = {'PEuni':'', 'PEdeg':'', 'SEuni':'', 'SEdeg':''}
	
	if PESECOMBO:
		rf['SEuni'] = pref+'.unique.map'
		if STRATEGY!='unique' and STRATEGY!='bestno': 
			rf['SEdeg'] = pref+'.degenerate.map'
	else:
		rf = {'PEuni':pref+'.PE.unique.map', 'PEdeg':'', 'SEuni':'', 'SEdeg':''}
	
		if USESINGLETONS:
			rf['SEuni'] = pref+'.SE.unique.map'
			if STRATEGY!='unique' and STRATEGY!='bestno': 
				rf['PEdeg'] = pref+'.PE.degenerate.map'
				rf['SEdeg'] = pref+'.SE.degenerate.map'
		else:
			if STRATEGY!='unique' and STRATEGY!='bestno':
				rf['PEdeg'] = pref+'.PE.degenerate.map'
	
	
	# pref = btmapdir+sampleID
	# rf = {'PEuni':pref+'.PE.unique.map', 'PEdeg':'', 'SEuni':'', 'SEdeg':''}
	# 
	# if USESINGLETONS:
	# 	rf['SEuni'] = pref+'.SE.unique.map'
	# 	if not UNIQ: 
	# 		rf['PEdeg'] = pref+'.PE.degenerate.map'
	# 		rf['SEdeg'] = pref+'.SE.degenerate.map'
	# else:
	# 	if not UNIQ:
	# 		rf['PEdeg'] = pref+'.PE.degenerate.map'
	
	print ' - Processing sample %s...'%(sampleID)
	
	# 2 relevant file names
	# if using seq bounds file, distribute genotyping across processors
	# assumes sun grid engine
	motherfile = theresdir+sampleID+'.mother'   # genotype at all informative positions
	# nominalfile = theresdir+sampleID+'.nominal' # nominally significant snp variants
	
	motherfile += boundslabel+'.txt'
	# nominalfile += boundslabel+'.txt'
	
	# retrieve 
	PCRreplicatedReads = {'unique':'', 'degenerate':''}
	if REDUNDANCYFILTER:
		for ud in ['unique', 'degenerate']:
			afile = badReadDir+'%s_%s_redundantIDs.txt'%(sampleID,ud)
			PCRreplicatedReads[ud] = afile
			if not os.access(afile, os.F_OK): 
				pipeit('Warning: redundancy filter ON but file not found: %s'%(afile),1)
				
	# perform genotyping if this hasn't been done yet
	numHypos = 0 # number of informative positions; for multiple test correction
	
	# process likelihoods for this file if the file doesn't exist, if user specifies a redo or ? 
	# os.stat(motherfile).st_size == 0 # if file has no data
	if not os.access(motherfile, os.F_OK) or REDOSNPCALLING:
		numHypos = glimmrLikelihood(readfiles=rf, GR=GR, outfile=motherfile, expPhred=expPhred, maxreadlen=maxreadlen, Lambda=Lambda, tmpdir=tmpdir, bounds=seqbounds, verbose=VERBOSE, ploidy=1, libraryErrorFiles=PCRreplicatedReads, LIPieces=LIPieces)
	
	
# finished processing, so write to an outfile so parent process knows to launch another
createdir(completedir)
fin = open(completedir+'likelihood_%s.done'%(PID), 'w'); fin.close()

pipeit('Results for %s saved to %s'%(PID, resdir), 1)