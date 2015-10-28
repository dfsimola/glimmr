#! /usr/bin/env python

import os, re, sys, math, time, stat, time, random, subprocess, decimal
from operator import add,mul
from string import *

GLIMMR_VERSION = 1.7

version = """\nglimmrAccessories.py, version %s
 
 Requires Unix/Mac OS X/CYGWIN with Python 2.5+, bowtie v1 or v2 
 (preferably 64-bit for memory allocation > 4GB)
 
 Daniel F. Simola, PhD (simola@upenn.edu)
 Laboratory of Shelley L. Berger, PhD
 University of Pennsylvania
 October 2015
 
 Copyright (c) 2015, Daniel F. Simola and Shelley L. Berger, University of
 Pennsylvania. You may not use this file except in compliance with the terms of our
 License, GNU GENERAL PUBLIC LICENSE v3.0, a copy of which is located in the LICENSE file.
"""%(GLIMMR_VERSION)

"""%s"""%(version)

multichromtoken = '.'         # separator for --rest <chr1.chr2.....chrn> separation
nastr = '"NA"'                # represents missing value
nan = nastr
delim = '\t'
infy = decimal.Decimal('inf') # 

GLOB_FONT_FACE = 'Helvetica'
GLOB_FONT_SIZE = 18
GLOB_FORMAT = 'png'
GLOB_INF = 1e308

 # --------------------------------------------------------------------

makeDictionary = lambda mylist,value: dict(map(lambda x: (x,value), mylist))
def makeDict(mylist,value=None):
	return makeDictionary(mylist,value)

makeDictFromPairs = lambda mylist: dict(map(lambda x: (x[0],x[1]), mylist))


getLabel = lambda astr: '.'.join(astr.split('/')[-1].split('.')[:-1])
getPath = lambda astr: '/'.join(astr.split('/')[:-1])+'/'
getSuffix = lambda astr: astr.split('/')[-1].split('.')[-1]

vslice = lambda dat, j: map(lambda row: row[j], dat)

def valueWithinRange(interval, val):
	list(interval).sort()
	if val < interval[0] or val > interval[1]: return False
	else: return True

def intervalWithinRange(baseint, queryint):
	"""Returns 0,1,2,3 according to these rules:
	0: queryint is not within base int
	1: queryint is within base int
	2: queryint spans the beginning of baseint
	3: queryint spans the end of baseint
	baseint/queryint need not already be sorted
	ie, [1000, 2000] and [2000,1000] will return same value
	"""
	baseint = sorted(baseint); queryint = sorted(queryint)
	
	# if the end of the element comes before the start of the base
	if (queryint[1] < baseint[0] or baseint[1] < queryint[0]):
		return 0
	# if element lies inside the base region
	elif (queryint[0] >= baseint[0] and queryint[1] <= baseint[1]):
		return 1
	# if we span a region
	elif (queryint[0] < baseint[0] and queryint[1] >= baseint[0]):
		return 2
	else:
		return 3

def intervalWithin(baseint, queryint):
	"""Deprecated"""
	return intervalWithinRange(baseint, queryint)


def overlap(baseint, queryint):
	"Returns number of shared bases if overlapping, otherwise 0"
	
	baseint = sorted(baseint); queryint = sorted(queryint)
	
	# if the end of the element comes before the start of the base
	if (queryint[1] < baseint[0] or baseint[1] < queryint[0]):
		return 0
	# if element lies inside the base region
	elif (queryint[0] >= baseint[0] and queryint[1] <= baseint[1]):
		return abs(queryint[0]-queryint[1])
	# if base region lies inside the query region
	elif (baseint[0] >= queryint[0] and baseint[1] <= queryint[1]):
		return abs(baseint[0]-baseint[1])
	# if we span a region
	elif (queryint[0] < baseint[0] and queryint[1] >= baseint[0]):
		return abs(baseint[0] - queryint[1])
	else:
		return abs(baseint[1] - queryint[0])



def maxna(lst):
	d = filterstr(nastr)(lst)
	if len(d)>0: return max(d)
	return nastr

def minna(lst):
	d = filterstr(nastr)(lst)
	if len(d) > 0: return min(d)
	return nastr


def isdir(path):
	try: return stat.S_ISDIR(os.stat(path)[stat.ST_MODE])
	except OSError: return False


def mergeSamFiles(lst,outfile,header=False):
	for x in lst:
		assert os.access(x, os.F_OK), 'Cannot access %s'%(x)
	
	f1 = lst.pop(0)
	# assume first file has desired sam header
	if header==True:
		os.system('cp "%s" "%s"'%(f1,outfile))
	else:
		fho = open(outfile,'a')
		fh = open(f1, 'r')
		for line in fh:
			if line[0] != '@': print >> fho, line,
		fh.close()
		fho.close()
		
	# read subsequent lst files and concatenate non-header lines to outfile
	fho = open(outfile,'a')
	for l in lst:
		fh = open(l, 'r')
		for line in fh:
			if line[0] != '@': print >> fho, line,
		fh.close()
	fho.close()
	

def getFilename(astr):
	return astr.split('/')[-1]


def isnan(x):
    return str(x) == str(1e400*0)


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

null2na = lambda lst: map(lambda x: len(x) and x or nastr, lst)

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
						# fix oddities in item
						item = item.replace('+','\+')
						
						# print 'TEST', item,'vs',f
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
						# else: print 'no match', item,'vs',f,re.match('.*%s.*'%(item), f)
				# else:
					# print 'ODDITY', f
					
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

# def readList(fname, header=False):
# 	crap,thelist,crap = readTable(fname, header=header)
# 	return flatten(thelist)
def readList(fname, dtype=None, delim='\n'):
	fh = open(fname)
	lst = None
	if dtype: lst = map(lambda x: dtype(x[:-1]), fh)
	else: lst = map(lambda x: x[:-1], fh)
	fh.close()
	if delim != '\n': return lst[0].split(delim)
	return lst

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


def unique(lst):
	"""Returns a subset of lst containing the unique elements."""
	d = {}
	for i in lst:
		if i not in d: d[i] = 0
		else: d[i] += 1
	e = d.keys(); e.sort()
	return e


def unzip(thread):
	"""Reverse of the built-in method zip."""
	try:
		if len(thread):
			return [[thread[i][j] for i in range(len(thread))] for j in range(len(thread[0]))]
		else: return []
	except ValueError: return []


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

def argmin(v):
	"Return argmin of a list"
	mi = 0
	mv = v[0]
	for i in range(1,len(v)):
		if v[i] < mv:
			mv = v[i]
			mi = i
	return mi


def slash(astr):
	if not len(astr): return '/'
	elif astr[-1] !='/': astr+='/'
	elif astr[len(astr)-2:len(astr)] == '//': astr = astr[:-1]
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
sumna = lambda lst: sum(filterstr(nastr)(map(floatna, lst)))

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

stdev = sd

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

def quantiles(v=[], na=True, nastr=nastr, Q=[.05,.25,.5,.75,.95]):
	v2 = v
	if na: v2 = filterstr(nastr)(v)
	return dict( map(lambda q: (q, percentile(v2,q,na=False)), Q) )

getQuants = lambda dat: ','.join([str(q)+':'+str(round(percentile(dat, q),4)) for q in [0.05, 0.25, 0.5, 0.75, 0.95]])

# list operations
head = lambda x: x[0]
tail = lambda x: x[len(x)-1]
car = lambda x: x[0]
cdr = lambda x: x[len(x)-1]

# arithmetic operations with missing values
addna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) + float(v)
subna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) - float(v)
multna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) * float(v)
divna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) / float(v)
roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))
sumna = lambda lst: sum(map(float, filterstr(nastr)(lst)))
sqrtna = lambda x: x==nastr and nastr or math.sqrt(x)
roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))

addna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(m) + float(v)
subna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(v) - float(m)
multna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(v) * float(m)
divna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(v) / float(m)

r1 = lambda x: roundna(1)(x)
r2 = lambda x: roundna(2)(x)
r3 = lambda x: roundna(3)(x)
r4 = lambda x: roundna(4)(x)
r5 = lambda x: roundna(5)(x)

def makeRange(xmin,xmax,interval=None,n=None):
	if n == None: n = int(abs(xmax-xmin)/float(interval))
	else: 
		interval = (xmax-xmin)/float(n)
		interval += interval/float(n-1)
	return [xmin+interval*i for i in range(n)]


def histogram(dat, nbins=100, bins=[], categorical=False, yfreq=False, logscale=False, logbase=10, shift=0, dtype=None):
	import math
	if logscale: dat = map(lambda x: math.log(x,logbase), dat)
	
	if len(bins) > 0: 
		# print 'problem?'
		nbins = len(bins)
		givenbins = True
	else: givenbins = False
	
	# need to create a histogram
	filt = filterstr(nastr)(dat)
	if dtype: filt = map(dtype, filterstr(nastr)(dat))
	dat = None
	if not len(filt): return [[0,0]]
	
	counts = [0 for i in range(nbins)]
	if not givenbins and not categorical: bins = makeRange(min(filt), max(filt), n=nbins)
	elif not givenbins: bins = unique(filt)
	
	# minor bin correction for plotting
	if shift != 0: bins = map(lambda x: x+shift, bins)
	
	# this is like O(N2) too slow
	if not categorical:
		for j in range(len(filt)):
			k = 0
			# print j, len(filt)
			while k < nbins and float(filt[j]) > bins[k]: 
				k +=1
			if k == nbins: k -= 1
			counts[k] += 1
	else:
		bindct = {}
		for b in bins: bindct[b] = 0
		bk = bindct.keys()
		for j in range(len(filt)): bindct[filt[j]] += 1
		# convert back to a histogram
		bins = []; counts = []
		for key in bindct.keys(): bins += [key]; counts += [bindct[key]]
	
	tmp = []
	intnbins = int(nbins)
	# counts or frequency?
	if yfreq:
		tot = float(sum(counts))
		if logscale: tmp = [[logbase**bins[k],counts[k]/tot] for k in range(intnbins)]
		else: tmp = [[bins[k],counts[k]/tot] for k in range(intnbins)]
	if logscale: tmp = [[logbase**bins[k],counts[k]] for k in range(intnbins)]
	else: tmp = zip(bins,counts)
	
	if categorical: tmp = sorted(tmp)
	
	return tmp


def reHist(hist,bins, freq=False):
	# aggregates an existing histogram in terms of new bins
	
	newhist = dict( map(lambda x: (x,0), bins) )
	nbins = len(bins)
	# which bins corresponds to this key?
	binmap = {}
	
	for k in transpose(hist)[0]:
		newbini = 0
		while newbini < nbins and bins[newbini] < k:
			newbini += 1
		if newbini >= nbins: newbini = nbins-1
		binmap[k] = bins[newbini]
	
	for k,v in hist:
		newhist[binmap[k]] += v
	
	ret = sorted(newhist.items())
	
	#
	if freq:
		tot = float(sum(transpose(ret)[1]))
		# if logscale: tmp = [[logbase**bins[k],counts[k]/tot] for k in range(intnbins)]
		# else: tmp = [[bins[k],counts[k]/tot] for k in range(intnbins)]
		ret = [[bins[k],ret[k][1]/tot] for k in range(len(bins))]
	
	return ret



log = lambda x: math.log(x, 10)
lg = lambda x: math.log(x, 2)
ln = lambda x: math.log(x, math.e)

# dictionary keyed by ascii character returning phred q-score
# ascii 33 to 126 - full range
fastqchars = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '\{', '\|', '\}', '~']
A2Q = dict(zip(fastqchars, range(len(fastqchars))))

# string conversion to integer (potentially a long)
encodeNt = lambda x: x=='A'and'1' or x=='T'and'2' or x=='C'and'3' or x=='G'and'4' or x=='N'and'5'
encodeSeq = lambda seq: int(''.join(map(encodeNt, seq)))
encodeSeqStr = lambda seq: ''.join(map(encodeNt, seq))
decodeNt = lambda x: x=='1'and'A' or x=='2'and'T' or x=='3'and'C' or x=='4'and'G' or x=='5'and'N'
decodeSeq = lambda val: ''.join(map(decodeNt, str(val)))

encodeQ = dict(zip(fastqchars, map(str,range(10,len(fastqchars)+10)))) # every character is 2 integer digits
decodeQ = dict(zip(map(str,range(10,len(fastqchars)+10)), fastqchars))
encodeQual = lambda qual: int(''.join(map(lambda x: encodeQ[x], qual)))
pairUP = lambda qual: map(lambda x: qual[x-1]+qual[x] ,range(1,len(qual),2))
decodeQual = lambda qual: ''.join(map( lambda x: decodeQ[x], pairUP(str(qual)) ))

isNT = lambda x: (x=='A' or x=='T' or x=='C' or x=='G' and True) or False

def rc(seq=''):
	"Translates a nucleotide sequence (including degenerate symbols) into its reverse complement. Always returns uppercase."
	
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
		
	rcseq = map(srev, seq)
	rcseq.reverse()
	rcseq = ''.join(rcseq)
	return rcseq

#
def threadListPairs(A,B,inner=True,lena=1):
	pairs = []
	genes = []
	Adct = {}
	for i in range(len(A)): 
		# print 'storing', A[i][1], 'into key', A[i][0]
		Adct[A[i][0]] = A[i][1]
	
	for i in range(len(B)):
		name = B[i][0]
		if name in Adct:
			pairs += [[Adct[name], B[i][1]]]
			# print 'adding', pairs[-1]
			genes += [name]
		elif not inner:
			if lena > 1: pairs += [([nastr for x in xrange(lena)], B[i][1])]
			else: pairs += [(nastr, B[i][1])]
			genes += [name]
	# print 'returning', pairs
	return genes, pairs
#
def threadSetPairs(lst,inner=True,flat=True):
	if len(lst) ==1: 
		n,d = unzip(lst[0])
		d = [[x] for x in d]
		return n,d
	if not len(lst): return [],[]
	
	# if list contents are themselves lists or tuples, make sure they are listed
	if type(lst[0][1][1]) == type(()):
		# print 'we have tuples', lst[0][1]
		# print lst[0]
		lst[0] = [[lst[0][i][0], list(lst[0][i][1])] for i in xrange(len(lst[0]))]
		lst[1] = [[lst[1][i][0], list(lst[1][i][1])] for i in xrange(len(lst[1]))]
	
	ans = threadListPairs(lst[0], lst[1], inner)
	# print 'answer', ans
	
	lena = 2
	for i in range(2,len(lst)):
		ans = threadListPairs(zip(*ans), lst[i], inner, lena=lena)
		lena += 1
	# final flattening of data
	if flat:
		# print 'flattening', ans[1]
		# return [ans[0],[flatten(row) for row in ans[1]]]
		return [ans[0],map(flatten,ans[1])]
	else:
		return [ans[0],ans[1]]
#
def joinSet(lst,inner=True,flat=True):
	return threadSetPairs(lst,inner,flat=flat)






# smoothing methods
def kernelsmoothlist(A, w=5, kernel='normal', weights=None, naconvert=None):
	"""Apply a smoothing kernel along the rows of a list A, using window 
	 of length w. A Gaussian distribution is laid over the window, s.t.:
	 Ai is a sum of the n adjacent Aj values weighted by the Gaussian."""
	
	
	if naconvert != None:
		A = map(lambda x: x != nastr and x or naconvert, A)
	
	# if no smoothing return original data vector
	if not w: return A
	
	# make sure n is odd, nb division is integer division
	if (w/2)*2 == w: w -= 1
		
	weight = range(w) # weight vector of the kernel
	# populate W using Binomial distribution
	# idea is that we have a window of length w-1, and we can observe
	# any event from 0 to w-1, where P(event) = 1/2
	# this gives us a discrete approximation to a gaussian distribution for
	# w events
	
	
	if kernel == 'normal':
		if w > 25:
			weight = map(lambda i: binomial(w-1,i,.5,approx=1), range(w))
			# we are approximating so renormalize
			z = float(sum(weight))
			weight = map(divna1(z),weight)
			# print 'weight vector', sum(weight)
		else:
			weight = map(lambda i: binomial(w-1,i,.5,approx=0), range(w))
	elif kernel == 'uniform':
		weight = map(lambda x: 1/float(w), range(w))
	
	# if weights:
	# 	if len(weights) != len(weight): 
	# 		print 'Error provided weights (n=%s) do not match bin size w=%s'%(len(weights),w)
	# 		return nastr
	# 	
	# 	print 'old weights:', weight
	# 	weight = map(ut.multna, weight, weights)
	# 	print 'new weights:', weight
		
	
	#print >> sys.stderr, "sum of probs", sum(weight)
	
	B = A[:] # copy A, to store transformed data
	n = len(A)
	for i in range(n):
		# build kernel
		half = int(math.floor(w/2.0)) # half of the kernel window
		lbuf = half - i # our left buffer for the row
		kernel = range(w)
		
		# kw = [1 for zz in range(w)]
		
		# local kernel
		lw = weight[:]
		if weights:
			mx = maxna(weights)
			if mx == nastr: return A
			# print 'max weight', 2*mx
			for k in kernel:
				# update the kernel based on user input
				val = nastr
				if lbuf > 0:
					val = weights[i]
				elif k - (n-i+1) < 0: # rbuf < n
					val = weights[i+k-half]
				else:
					val = weights[n-1]
				
				if val in [nastr,0]: 
					# print 'MAX', mx
					val = mx
				lw[k] *= 1/val
				
		# normalize back to 1
		lw = map(lambda x: x/float(sum(lw)) , lw)
		# print 'LW', lw, A
		for k in kernel:
			# if we are near the first or last columns of a row,
			# pad the kernel with the terminal values. i.e the first 
			# window at Ai0 (n=3) will be [Ai0, Ai0, Ai0, Ai1, Ai2]
			if lbuf > 0:
				# kernel is past left end of row
				kernel[k] = lw[k]*A[i]
				lbuf -= 1
			elif k - (n-i+1) < 0: # rbuf < n
				# kernel is wholly within the row
				# print 'here', k, i+k-half, weight[k], A[i+k-half], len(A)
				kernel[k] = lw[k]*A[i+k-half]
			else:
				# kernel will go past right end of row
				kernel[k] = lw[k]*A[n-1]
		# now store new Bij value
		B[i] = sum(kernel)
	return B


def kernelsmooth(A, w=5, n=1, kernel='normal', weights=None, naconvert=None): 
	if len(A) < w: w = 0
	A = kernelsmoothlist(A,w,kernel,weights,naconvert)
	for i in range(n-1):
		A = kernelsmoothlist(A,w,kernel,weights,naconvert)
	return A


################
# smoothing methods
def kernelsmoothlistweights(A, w=5, kernel='normal', weights=None, naconvert=None):
	"""Apply a smoothing kernel along the rows of a list A, using window 
	 of length w. A Gaussian distribution is laid over the window, s.t.:
	 Ai is a sum of the n adjacent Aj values weighted by the Gaussian."""
	
	
	if naconvert != None:
		A = map(lambda x: x != nastr and x or naconvert, A)
	
	# if no smoothing return original data vector
	if not w: return A
	
	# make sure n is odd, nb division is integer division
	if (w/2)*2 == w: w -= 1
		
	weight = range(w) # weight vector of the kernel
	# populate W using Binomial distribution
	# idea is that we have a window of length w-1, and we can observe
	# any event from 0 to w-1, where P(event) = 1/2
	# this gives us a discrete approximation to a gaussian distribution for
	# w events
	
	
	if kernel == 'normal':
		if w > 25:
			weight = map(lambda i: binomial(w-1,i,.5,approx=1), range(w))
			# we are approximating so renormalize
			z = float(sum(weight))
			weight = map(divna1(z),weight)
			# print 'weight vector', sum(weight)
		else:
			weight = map(lambda i: binomial(w-1,i,.5,approx=0), range(w))
	elif kernel == 'uniform':
		weight = map(lambda x: 1/float(w), range(w))
	
	# if weights:
	# 	if len(weights) != len(weight): 
	# 		print 'Error provided weights (n=%s) do not match bin size w=%s'%(len(weights),w)
	# 		return nastr
	# 	
	# 	print 'old weights:', weight
	# 	weight = map(ut.multna, weight, weights)
	# 	print 'new weights:', weight
		
	
	#print >> sys.stderr, "sum of probs", sum(weight)
	
	B = A[:] # copy A, to store transformed data
	n = len(A)
	for i in range(n):
		# build kernel
		half = int(math.floor(w/2.0)) # half of the kernel window
		lbuf = half - i # our left buffer for the row
		kernel = range(w)
		
		# kw = [1 for zz in range(w)]
		
		# local kernel
		lw = weight[:]
		if weights:
			mx = maxna(weights)
			if mx == nastr: return A
			# print 'max weight', 2*mx
			for k in kernel:
				# update the kernel based on user input
				val = nastr
				if lbuf > 0:
					val = weights[i]
				elif k - (n-i+1) < 0: # rbuf < n
					val = weights[i+k-half]
				else:
					val = weights[n-1]
				
				if val in [nastr,0]: 
					# print 'MAX', mx
					val = mx
				lw[k] *= 1/val
				
		# normalize back to 1
		lw = map(lambda x: x/float(sum(lw)) , lw)
		
		for k in kernel:
			# if we are near the first or last columns of a row,
			# pad the kernel with the terminal values. i.e the first 
			# window at Ai0 (n=3) will be [Ai0, Ai0, Ai0, Ai1, Ai2]
			if lbuf > 0:
				# kernel is past left end of row
				kernel[k] = lw[k]*A[i]
				lbuf -= 1
			elif k - (n-i+1) < 0: # rbuf < n
				# kernel is wholly within the row
				# print 'here', k, i+k-half, weight[k], A[i+k-half], len(A)
				kernel[k] = lw[k]*A[i+k-half]
			else:
				# kernel will go past right end of row
				kernel[k] = lw[k]*A[n-1]
		# now store new Bij value
		B[i] = sum(kernel)
	return B

def kernelsmoothweights(A, w=5, n=1, kernel='normal', weights=None, naconvert=None): 
	if len(A) < w: w = 0
	A = kernelsmoothlist(A,w,kernel,weights,naconvert)
	for i in range(n-1):
		A = kernelsmoothlist(A,w,kernel,weights,naconvert)
	return A




