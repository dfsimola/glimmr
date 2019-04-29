#! /usr/bin/env python

import os, re, sys, math, time, stat, time, random, subprocess, decimal, glob
from operator import add,mul
from string import *

GLIMMR_VERSION = 1.8

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

# parts of a filename
getLabel = lambda astr: '.'.join(astr.split('/')[-1].split('.')[:-1])
getFilename = lambda astr: astr.split('/')[-1]
getDirectory = lambda astr: astr.split('/')[-2]
getPath = lambda astr: '/'.join(astr.split('/')[:-1])+'/'
getSuffix = lambda astr: astr.split('/')[-1].split('.')[-1]


vslice = lambda dat, j: map(lambda row: row[j], dat)

def slice(dat, rows=[], cols=[]):
	if rows == [] and cols == []: return dat
	
	if rows == [] and len(cols):
		return map(lambda row: [row[colj] for colj in cols], dat)
	elif cols == [] and len(rows):
		return [dat[r] for r in range(len(dat)) if r in rows]
	elif len(rows) and len(cols):
		return map(lambda row: [row[colj] for colj in cols], [dat[r] for r in range(len(dat)) if r in rows])

def multislice(dat,lst):
	return map(lambda row: [row[colj] for colj in lst], dat)


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

# SET OPERATIONS
# --------------
def intersection(a=[], b=[]):
	"""Returns a list containing those items common to both a and b.
	"""
	# return filter(lambda item: item in a, b) # O(n^2) time
	
	# faster version for big items
	ad = {}
	for item in a: ad[item] = 1
	return filter(lambda item: item in ad, b) # O(n) time

def setIntersection(a=[], b=[]):
	return intersection(a,b)

def setUnion(a=[],b=[],sorted=True):
	return unique(a + b, sorted=sorted)

def setDifference(a=[], b=[]):
	"""Setwise a - b: returns a list containing those items unique to a.
	"""
	# return filter(lambda item: item not in b, a)
	
	# faster version for big items
	bd = {}
	for item in b: bd[item] = 1
	return filter(lambda item: item not in bd, a)


def setGroup(lst,op):
	"""Perform the union/intersection/difference of all lists within lst."""
	if len(lst) == 0: return lst
	elif len(lst) == 1: return lst[0]
	base = op(lst[0],lst[1])
	for i in range(2,len(lst)): base = op(base,lst[i])
	return base

def groupUnion(lst): return setGroup(lst,setUnion)
def groupDifference(lst): return setGroup(lst,setDifference)
def groupIntersection(lst): return setGroup(lst,setIntersection)



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

def retypeArray(X,type):
	"""Typecast every element in X to type."""
	X2 = X[:]
	for i in range(len(X)):
		X2[i] = map(type, X[i])
	return X2

def castemat(m, cast):
	return retypeArray(m,cast)


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

def slash(astr):
	if not len(astr): return '/'
	elif astr[-1] !='/': astr+='/'
	elif astr[len(astr)-2:len(astr)] == '//': astr = astr[:-1]
	return astr

# READ METHODS
# ------------
def getFiles(indirfile, include=[], exclude=[], literal=False):
	
	if isdir(indirfile) or os.access(indirfile,os.F_OK):
		if type(include)!=type(list()): include = [include]
		if type(exclude)!=type(list()): exclude = [exclude]
		return getContents(indirfile, typ='files', include=include, exclude=exclude, literal=literal)
	
	return glob.glob(indirfile)

def getDirectories(indirfile, include=[], exclude=[], literal=False):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	files = getContents(indirfile, typ='dirs', include=include, exclude=exclude, literal=literal)
	return map(slash, files)

def getContents(indirfile, typ='all', include=[], exclude=[], literal=False):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	keepfiles = []
	if len(include) == 0: include = ['']
	wildcards = ['-', '+', ':', '*', '.']
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
						# escape wild cards
						if literal:
							for esc in wildcards: item = item.replace(esc, '\\'+esc)
						
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
					for esc in wildcards:
						item = item.replace(esc, '\\'+esc)
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

def getContentsDeprecated(indirfile, typ='all', include=[], exclude=[]):
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
	lst = map(lambda x: x.strip(delim), lst)
	# if delim != '\n': return lst[0].split(delim)
	return lst

# def printList(lst, title="", delim="\n", pipe=sys.stdout, file='', sort=False, newline='\n', ioMethod='w'):
# 	if file: pipe = open(file,ioMethod)
# 	if sort: lst.sort()
# 	if title: 
# 		lst = [title]+lst#print >> pipe, title
# 	for i in range(len(lst)):
# 		if i == len(lst)-1:
# 			print >> pipe, '"%s"'%(str(lst[i])),
# 		else:
# 			print >> pipe, str(lst[i])+delim,
# 	if newline: print >> pipe, newline,
# 	if file: pipe.close()

def printList(lst, title="", delim="\n", pipe=sys.stdout, file='', sort=False, newline='\n', ioMethod='w'):
	if file: pipe = open(file,ioMethod)
	if sort: lst.sort()
	if title: 
		lst = [title]+lst#print >> pipe, title
	for i in range(len(lst)):
		if i == len(lst)-1:
			if newline: print(str(lst[i]), file=pipe)
			else: 
				print(str(lst[i]), file=pipe)
				# print 'test "%s"'%(lst[i])
		else:
			print(str(lst[i])+delim, file=pipe, end='')
	# if newline: print >> pipe, newline,
	if file: pipe.close()

# def readTable(filename, header=True, rownames=True, delim="\t", comments=True, keepKey=False):
# 	fh = open(filename)
# 	table = fh.readlines()
# 	fh.close()
# 	if comments: table = trimComments(table)
	
# 	tabHeader = []; rowIDs = []; data = []
# 	if header: 
# 		tabHeader = table.pop(0).strip('\n').split(delim)
# 		if not keepKey: tabHeader.pop(0)
# 		# else: tabHeader = table.pop(0).strip('\n').split(delim)
# 	for line in table:
# 		sline = line.strip('\n').split(delim)
# 		if rownames:
# 			rowIDs.append(sline.pop(0))
# 		data += [sline]
# 	return data, rowIDs, tabHeader

# changed keepKey default to True - watch for bugs!
def readTable(filename, header=True, rownames=True, delim="\t", newline='\n', comments=True, keepKey=False, fill=nastr, nrows=-1, dtype=None):
	fh = open(filename)
	table = []; ct = 0
	for row in fh:
		if nrows > 0 and ct > nrows: break
		table += [row]
		ct += 1
	# table = fh.readlines()
	fh.close()
	if not len(table): return [],[],[]
	if comments: table = trimComments(table)
	
	tabHeader = []; rowIDs = []; data = []
	if header: 
		tabHeader = table.pop(0).strip(newline).split(delim)
		if not rownames: keepKey = True
		if not keepKey: tabHeader.pop(0)
		# else: tabHeader = table.pop(0).strip('\n').split(delim)
	for line in table:
		sline = line.strip(newline).split(delim)
		if fill != None:
			sline = list(map(lambda x: x=='' and fill or x, sline))
		# print 'sline', sline
		if rownames:
			rowIDs.append(sline.pop(0))
		data += [sline]
	
	if dtype: data = castemat(data,dtype)
	
	data = [list(row) for row in data]

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
		if VERBOSE: print('- Loading %s'%(f))
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
			print(theid, ": This key has already been used to store another fasta entry!")
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
	
	if verbose: print(shorts)
	if verbose: print(taxa, alilen)
	if verbose: print(map(len, seqs))
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
	e = sorted(d.keys())#; e.sort()
	return e


def unzip(thread):
	"""Reverse of the built-in method zip."""
	t2 = list(thread)
	try:
		if len(t2):
			return [[t2[i][j] for i in range(len(t2))] for j in range(len(t2[0]))]
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


def floatnastr(x,na=nastr):
	y = None
	try: y = floatna(x)
	except ValueError: y = str(x)
	return y

def floatnacrap(x,na=nastr):
	if x == 'nan': return na
	y = na
	try: y = float(x)
	except ValueError: return na
	return y

floatorstr = floatnacrap
floator0 = lambda x: floatorstr(x,0)


def factorize(lst):
	# hash strings to integers
	dct = dict([(lst[i],i) for i in range(len(lst))])
	# return integers in same order as input
	return [dct[lst_i] for lst_i in lst]

def sign(x, returnBool=True):
	# return True or False 
	return (returnBool and floatna(x) > 0 and True or False) or (floatna(x) > 0 and '+' or '-')

# list manipulation
head = lambda x: x[0]
tail = lambda x: x[len(x)-1]
car = lambda x: x[0]
cdr = lambda x: x[len(x)-1]
append = lambda astr: lambda lst: map(lambda x: str(x)+str(astr), lst)





def factit(x):
	"""Iterative implementation of factorial function."""
	if x <= 1: return 1
	a = [0 for i in range(x+1)]
	a[0] = 1; a[1] = 1
	for i in range(2,x+1): a[i] = i*a[i-1]
	return a[x]

nchoosek = lambda n,k: factit(n)/(factit(k)*factit(n-k)) # binomial coefficient
prod = lambda lst: reduce(mul, lst) # product function

def binomial(n,x,p, approx=False):
	"""Returns the binomial probability of x successes in n trials"""
	
	if approx:
		# Normal(np, np(1-p))
		import maths
		return maths.GaussianPDF(x, n*p, n*p*(1-p))
	
	c = factit(n)/float(factit(x)*factit(n-x))
	return c*(p**x)*((1-p)**(n-x))


# def avg(lst):
# 	try: return reduce(add, lst)/float(len(lst))
# 	except TypeError: return nastr

def avg(v, w=[], na=True, ret=nastr, nastr=nastr, cast=floatna):
	"""When using a weight vector w, this assumes that NAs
	in list v have weight 0.0.
	"""
	if len(v) == 0: return ret
	if len(w) > 0 and len(w) != len(v): return ret
	v2 = v[:]
	if na: v2 = filter(eq(nastr), v)
	
	v2 = list(map(cast, filter(eq(None), v2)))
	if len(v2) == 0: return ret
	elif len(unique(v2))==1 and v2[0] == nastr: return ret
	if len(w):
		return sum([v2[i]*float(w[i]) for i in range(len(v2))])
	else:
		try:
			return sum(v2)/float(len(v2))
		except TypeError:
			raise TypeError('utilities.py | avg(): type error on input list: %s...'%(v2[0:10]))


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

def stdev(v, na=True, ret=nastr):
	var = variance(v, na=na, ret=nastr)
	return (var != nastr and math.sqrt(var)) or nastr

def sd(lst):
	try: return math.sqrt(variance(lst))
	except TypeError: return nastr

def sd(v, na=True, ret=nastr):
	return stdev(v,na,ret)

def MAD(v,na=True,ret=nastr):
	"""Median absolute deviation"""
	v2 = v
	if na: v2 = filterstr(nastr)(v)
	gmed = median(v2)
	return median([abs(vi - gmed) for vi in v2])

def cv(v, na=True, ret=nastr, biased=False):
	v2 = na and filterstr(nastr)(v) or v
	V = multna(100,divna(stdev(v2,na=False,ret=ret), avg(v2,na=False,ret=ret)))
	return V==nastr and nastr or (biased and V) or (1 + (1/float(4*len(v2)))) * V

def stderr(v, na=True, ret=nastr, cast=float):
	"""Compute standard error of the mean"""
	v2 = v[:]
	if na: v2 = filterstr(nastr)(map(cast, filterstr(nastr)(v)))
	n = len(v2)
	return divna(stdev(v2,na=False,ret=ret),math.sqrt(n))

powna = lambda b: lambda e: ((b==nastr or e==nastr) and nastr) or float(b)**float(e)

def zfunc(X,center=avg,scale=stdev):
	m = center(X)
	sd = scale(X)
	return map(lambda xi: divna(subna(xi,m),sd), X)


def percentile(v=[], k=0.5):
	"""
	Return the value of array that corresponds to the kth percentile of the data.
	"""
	if len(v) <= 0: return 0.0
	temp = sorted(list(map(float, v)))
	# temp.sort()
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
	return dict( map(lambda q: (q, percentile(v2,q)), Q) )

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

arcsin = lambda x: (x==nastr and x) or math.asin(x)

glg = lambda l: lambda x: (x==nastr and x) or math.log((float(x)+math.sqrt(float(x)**2+l))/2.0, 2)
slg = lambda c: lambda x: (x==nastr and x) or math.log(addna(x,c), 2)
slog = lambda c: lambda x: (x==nastr and x) or math.log(addna(x,c), 10)
lg = lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),2)
log10 = lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),10)
log = log10
ln = lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),math.e)
log2 = lg
logb = lambda b: lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),b)

absna = lambda x: (x==nastr and x) or abs(x)

lgm = lambda x, y: math.log(float(x),2)+math.log(float(y),2)
lg2 = lambda n,d: math.log(float(n),2)-math.log(float(d),2)


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

def complement(seq=''):
	srev = lambda c:\
		(c=='A'and'T') or (c=='T'and'A') or (c=='C'and'G') or (c=='G'and'C') or (c=='U'and'A') or\
		(c=='a'and't') or (c=='t'and'a') or (c=='c'and'g') or (c=='g'and'c') or (c=='u'and'a') or\
		(c=='-'and'-') or (c=='M'and'K') or (c=='K'and'M') or (c=='R'and'Y') or (c=='Y'and'R') or\
		(c=='-'and'-') or (c=='m'and'k') or (c=='k'and'm') or (c=='r'and'y') or (c=='y'and'r') or\
		(c=='S'and'S') or (c=='V'and'B') or (c=='B'and'V') or (c=='W'and'W') or \
		(c=='s'and's') or (c=='v'and'b') or (c=='b'and'v') or (c=='w'and'w') or \
		(c=='H'and'D') or (c=='D'and'H') or 'N' or \
		(c=='h'and'd') or (c=='d'and'h') or 'n'
		
	return ''.join(map(srev, seq))

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
def threadListPairs(A,B,inner=True,lena=1,lenadct=1):
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
			if lena > 1: pairs += [([[nastr for y in xrange(lenadct)] for x in xrange(lena)], B[i][1])]
			else: pairs += [([nastr for y in xrange(lenadct)], B[i][1])]
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
	
	# print 'HEY', len(lst[0][1][1])
	ans = threadListPairs(lst[0], lst[1], inner, lena=1, lenadct=len(lst[0][1][1]))
	# print 'answer', ans
	
	lena = 2
	for i in range(2,len(lst)):
		ans = threadListPairs(zip(*ans), lst[i], inner, lena=lena, lenadct=len(ans[0][1]))
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
	"""input is a list of lists, e.g. [A, B, C], where A = [(name1, val1), (name2, val2), ...]"""
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




#
filterstr = lambda s: lambda v: filter(eq(s), v)

def filterStrFromLists(meta, thestr=nastr):
	# this has potential bug when given a completely na list
	n = len(meta); m = len(meta[0])
	threaded = [ [meta[i][j] for i in range(n)] for j in range(m)]
	# filter by string
	nathreaded = []
	for i in range(m):
		toss = False
		for j in range(n):
			if threaded[i][j] == thestr: 
				toss = True
				break
		if not toss:
			nathreaded.append(threaded[i])
	# now unthread
	epi = []
	m = len(nathreaded)
	if m == 0: return [[],[]]
	
	n = len(nathreaded[0])
	return [[nathreaded[i][j] for i in range(m)] for j in range(n)]
	#return nathreaded

def rank(x,na=True,ties='split'):
	if na: x = filterstr(nastr)(x)
	xr = [[floatna(x[i]),i,0] for i in range(len(x))] # i is original ordering
	xr.sort() # sort by value
	
	# assign ranks by labeling with sequential integer
	xr = [[xr[i][1],xr[i][0],i] for i in range(len(xr))] # i is ranking
	
	# print 'sorted by value', xr
	
	if ties == 'split':
		# split ranks on ties
		rank = 1; rankcount = 1
		i = 1
		xr[0][2] = rank
		while i < len(xr):
			# print 'current:', xr[i],
			if xr[i][1] == xr[i-1][1]:
				# rank += 1
				rankcount += 1
				# print '1: rankcount', rankcount, 'x',xr[i][1], 'y', xr[i-1][1]
			elif rankcount > 1: 
				# set all previous indices to rank/rankcount
				
				for j in range(i-rankcount,i):
					xr[j][2] = rank + 1/float(rankcount)
				
				# print '2: rankcount', rankcount, 'rat',rank/float(rankcount)
				
				rank += 2
				# reset ranking
				rankcount = 1
				
			else: 
				rank += 1
				xr[i][2] = rank
				# print '3: rankcount', rankcount
				rankcount = 1
				
			i += 1
	# for i in range(len(xr)):
		# print 'NOW:', xr[i]
	
	xr.sort() # revert to original ordering
	return [xr[i][2] for i in range(len(xr))] # these are the ranks in original order

def variancem(v, mu, na=True, nastr=nastr, ret=nastr):
	s = 0.0; mu = float(mu)
	v2 = v[:]
	if na: v2 = map(float, filterstr(nastr)(v))
	if len(v2) == 1 or len(v2) == 0: return ret
	for i in v2: s += (i - mu)**2
	return s/float((len(v2)-1))

def variance(v, w=[], na=True, ret=nastr):
	return variancem(v, avg(v, w=w, na=na, ret=na), na=na, ret=ret)

def stdev(v, na=True, ret=nastr):
	var = variance(v, na=na, ret=nastr)
	return (var != nastr and math.sqrt(var)) or nastr

def covariance(x,y=[], na=True, nastr=nastr, ret=nastr):
	if not len(y):
		# x is threaded as pairs
		x2 = map(floatna, [x[i][0] for i in range(len(x))])
		y2 = map(floatna, [x[i][1] for i in range(len(x))])
		x = x2
		y = y2
	
	assert(len(x)==len(y)),\
	  'input vectors must match in length: x='+str(len(x))+' y='+str(len(y))
	
	if na: x,y = filterStrFromLists([x,y])
	if not len(x) or not len(y): return ret
	n = len(x); mx = avg(x); my = avg(y)
	return 1/float(n-1)*sum([(x[i]-mx)*(y[i]-my) for i in range(n)])

def correlation(x, y=[], na=True, ret=nastr, method='pearson', cast=floatna):
	if not len(y):
		# x is threaded as pairs
		x2 = []; y2 = []
		for i in xrange(len(x)):
			x2 += [x[i][0]]
			y2 += [x[i][1]]
		
		# x2 = [x[i][0] for i in range(len(x))]
		# y2 = [x[i][1] for i in range(len(x))]
		x = x2; y = y2

	def castWrapper(x):
		try: return cast(x)
		except ValueError: return ret
		
	if cast != None:
		x = map(castWrapper, x)
		y = map(castWrapper, y)
	
	if na: x,y = filterStrFromLists([x,y])
	if not len(x): return ret

	try:
		if method.lower() == 'spearman':
			x = rank(x,ties='split')
			y = rank(y,ties='split')
	except IndexError:
		print('RANKING PROBLEM in utilities.py::correlation(method="spearman")')
		# print 'X', x
		# print 'Y', y
		return ret

	
	
	try:
		R = divna(covariance(x,y,na=False), floatna(multna(stdev(x,na=False), stdev(y,na=False))))
		return R
	except ZeroDivisionError: return ret
	except IndexError: return ret


def corWithPvalue(x, y=[], na=True, method='pearson', pmethod='Z', tails='2', p='true', ret=nastr, invertpp=False):
	"""Returns a dictionary with correlation 'r', 'p-value', 'n', and 't'."""
	
	if not len(y):
		# x is threaded as pairs
		x2 = [x[i][0] for i in range(len(x))]
		y2 = [x[i][1] for i in range(len(x))]
		x = x2; y = y2
	if not len(x): return {'r':nastr, 't':nastr, 'P':nastr, 'n':0}
	x,y = filterStrFromLists([x,y])
	
	if not len(x) or not len(y): return {'r':nastr, 't':nastr, 'P':nastr, 'n':0}
	
	r = correlation(x,y,na=na,method=method)
		
	if pmethod == 'Z':
		# or Fisher transform
		n=nastr
		try:
			import maths
			n = len(x)
			Fr = 0.5*ln((1.+r)/(1.-r))
			Z = math.sqrt((n-3)/1.06)*Fr
			# print 'test', Fr, Z
			
			P = maths.GaussianCDF(Z,0,1)
			P = min(P,1-P)
			
			# P = None
			# if tails == '2':
			# 	P = maths.GaussianCDF(Z,0,1)
			# 	P = 2*min(P,1-P)
			# elif tails == '1' or tails == '-1':
			# 	P = maths.GaussianCDF(Z,0,1)
			# 	P = min(P,1-P)
			# elif tails.lower() == 'cdf':
			# 	P = maths.GaussianCDF(Z,0,1)
			# print 'post P', P
			
			if invertpp and r > 0: P = 1-P
			
			return {'r':r, 't':Fr, 'z':Z, 'P':P, 'n':n}
		
		except Exception(e):
			# print 'Exception', e
			return {'r':r, 't':nastr, 'z':nastr, 'P':nastr, 'n':n}
	elif pmethod == 'T':
		# Biometry p575-6
		# Standard t-test of product-moment correlation coefficient
		try:
			n = len(x)                          # sample size
			tstat = nastr
			pval = nastr
			if n >= 500:
				tr = r / math.sqrt((1-r**2)/(n-2))  # t-value
				pval = tdistPDF(tr,n-2,tails=tails) # p-value
				tstat = tr
			else:
				# from Biometry p575
				z = 0.5 * ln((1+r)/(1-r))
				zast = z - (3*z + r) / (4*(n-1))
				tz = zast * math.sqrt(n-1)
				# compare to t_alpha[inf]
				pval = tdistPDF(tz,1e3,tails=tails)
				tstat = tz
			if invertpp and r > 0: pval = 1-pval
			return {'r':r, 't':tstat, 'P':pval, 'n':n}
			# if p == 'true': return {'r':r, 't':tr, 'P':pval, 'n':n}
			# else: return {'r':r, 't':tz, 'P':pvalz, 'n':n}
		except ZeroDivisionError: return {'r':r, 't':nastr, 'P':nastr, 'n':n}
		except TypeError:         return {'r':r, 't':nastr, 'P':nastr, 'n':n}
		except ValueError:        return {'r':r, 't':nastr, 'P':nastr, 'n':n}

