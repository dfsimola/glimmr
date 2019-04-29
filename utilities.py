#! /usr/bin/env python

# utilities.py
# Daniel F. Simola
# Junhyong Kim Lab @ University of Pennsylvania
# July 2004
# Thu, 10 Apr 2008, 15:46:11

"""
"""

import sys, math, time, re, os, random
from operator import add, mul
from string import *


def isdir(f):
	return stat.S_ISDIR(os.stat(f)[stat.ST_MODE])


nastr = '"NA"'
nan = nastr
delim = '\t'
GLOB_FONT_FACE = 'Helvetica'
GLOB_FONT_SIZE = 18
GLOB_FORMAT = 'png'
GLOB_INF = 1e308


def overflow(x,ret=nastr):
	try:
		return x
	except OverflowError:
		return ret


def slash(astr):
	if not len(astr): return '/'
	elif astr[-1] !='/': astr+='/'
	elif astr[len(astr)-2:len(astr)] == '//': astr = astr[:-1]
	return astr

def deslash(astr):
	if astr[len(astr)-1] =='/': return astr[:-1]
	return astr


def absPath(apath):
	import os
	if apath[0] != '/': return slash(os.getcwd())+apath
	return apath


def splitString(x):
	return [x[i] for i in range(len(x))]

def maxlenval(tup):
	mlv = 0
	for i in tup:
		x = len(str(i))
		if x>mlv:mlv = x
	return mlv


enQuote = lambda q: lambda s: str(q)+str(s)+str(q)

def caseChangePos(seq, side='left', exclude=[]):
	"""Returns the first position where there is a change in the case of a string"""
	pos = 0
	ucase = seq[0].isupper()
	
	i = 1
	while i < len(seq) and (seq[i].isupper() == ucase or seq[i] in exclude):
		i += 1
	
	# we broke, so case changed at i
	if side == 'left' and i < len(seq): return seq[i:len(seq)]
	elif i < len(seq): return seq[0:i]
	return seq


def fixspacing(item,nsp):
	v = str(item)
	if len(v) == nsp: return v
	delt = nsp-len(v)
	if delt < 0:
		v = v[0:nsp]
	else:
		for i in range(delt): v+=' '
	return v


def polynomial(B, x):
	M = range(len(B)) # specify exponents of polynomial
	vv = map(lambda m, b: b*pow(x,m), M, B)
	return sum(vv) # final predicted response value

estimated = lambda B: lambda x: polynomial(B, x)
inv_estimated = lambda B: lambda x: 1/polynomial(B,x)
inv = lambda p: lambda x: 1/float(x**p)

#tosmall = lambda y: lambda x: (x=='NA' and 'NA') or (float(x)==0.0 and y) or x
null2na = lambda lst: map(lambda x: len(x) and x or nastr, lst)


subna = lambda m, v: ((v==nastr or m==nastr) and nastr) or float(m) - float(v)
addna = lambda m, v: ((v==nastr or m==nastr) and nastr) or float(m) + float(v)
divna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) / float(v)
multna = lambda m, v: ((v==nastr or m==nastr) and nastr) or float(m) * float(v)
powna   = lambda b: lambda e: ((b==nastr or e==nastr) and nastr) or float(b)**float(e)
pownaeb = lambda e: lambda b: ((b==nastr or e==nastr) and nastr) or float(b)**float(e)
modna = lambda m, v: ((v==nastr or m==nastr) and nastr) or float(m) % float(v)

def pownaof(b,e):
	try:
		return powna(b)(e)
	except OverflowError: return nastr

lenna = lambda lst: len(filterstr(nastr)(lst))

# maxna = lambda lst: max(filterstr(nastr)(lst))
def maxna(lst):
	d = filterstr(nastr)(lst)
	if len(d)>0: return max(d)
	return nastr

def minna(lst):
	d = filterstr(nastr)(lst)
	if len(d) > 0: return min(d)
	return nastr

# sumna = lambda lst: sum(filterstr(nastr)(map(floatna, lst)))

sumna = lambda lst: sum(map(float, filterstr(nastr)(lst)))

sqrtna = lambda x: x==nastr and nastr or math.sqrt(x)
roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))
# round = lambda e: lambda x: roundna(e)(x)

addna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(m) + float(v)
subna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(v) - float(m)
multna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(v) * float(m)
divna1 = lambda m: lambda v: ((v==nastr or m==nastr) and nastr) or float(v) / float(m)

# aliases
# add = lambda m, v: addna(m,v)
# sub = lambda m, v: subna(m,v)
# mult = lambda m, v: multna(m,v)
# div = lambda m, v: divna(m,v)

r1 = lambda x: roundna(1)(x)
r2 = lambda x: roundna(2)(x)
r3 = lambda x: roundna(3)(x)
r4 = lambda x: roundna(4)(x)
r5 = lambda x: roundna(5)(x)

arcsin = lambda x: (x==nastr and x) or 2*math.asin(math.sqrt(x))
logit = lambda x: (x==nastr and x) or log(x/(1-x))
# arcsin = lambda x: math.asin(math.sqrt(x))


glg = lambda l: lambda x: (x==nastr and x) or math.log((float(x)+math.sqrt(float(x)**2+l))/2.0, 2)
slg = lambda c: lambda x: (x==nastr and x) or math.log(addna(x,c), 2)
slog = lambda c: lambda x: (x==nastr and x) or math.log(addna(x,c), 10)
lg = lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),2)
log10 = lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),10)
log = log10
ln = lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),math.e)
log2 = lg
logb = lambda b: lambda x: (x==nastr and x) or (float(x)<=0.0 and nastr) or math.log(float(x),b)

absna = lambda x: (x==nastr and x) or abs(float(x))

lgm = lambda x, y: math.log(float(x),2)+math.log(float(y),2)
lg2 = lambda n,d: math.log(float(n),2)-math.log(float(d),2)
def lg1(n):
	try: return math.log(float(n),2)
	except OverflowError: return 0
	except ZeroDivisionError: return 0


logFunc = lambda f: lambda x: protectedLog(f,x)
def protectedLog(f, x):
	"""Log function wrapper that catches overflow errors assoc. 
	with large values"""
	try: return f(x)
	except OverflowError: return 0.0
	except ZeroDivisionError: return 0.0


postCat = lambda s: lambda t: str(t)+str(s)
preCat = lambda s: lambda t: str(s)+str(t)


Q2P = lambda e: (e == 0 and 1-1e-16) or 10**(-e/10.) # avoid 1
P2Q = lambda p: -10*log(p)


def writecol(data, v, j):
	if len(data) != len(v): return 0
	for i in range(len(data)): data[i][j] = v[i]
	return 1

def addcol(data, v, j):
	if len(data) != len(v): return 0
	for i in range(len(data)): data[i].insert(j,v[i])
	return 1


# type casting for missing values
castna = lambda typ: lambda x:  x == nastr or typ(x)
floatna = lambda x: x == nastr and nastr or x == '' or float(x)
intna = lambda x: (x == nastr or x == '') and nastr or int(x)

def floatnastr(x,na=nastr):
	y = None
	try: y = floatna(x)
	except ValueError: y = str(x)
	return y

def floatnacrap(x,na=nastr):
	y = None
	try: y = float(x)
	except ValueError: y = na
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

# time manipulation
# =================
def toseconds(astr, tomin=False):
	components = astr.split(' ')
	tsec = 0
	for c in components:
		if 'h' in c:
			tsec += float(c.split('h')[0]) * 3600. # hour to seconds
		elif 'm' in c:
			# print 'test', float(c.split('m')[0])
			tsec += float(c.split('m')[0]) * 60. # minutes to seconds
		elif 's' in c:
			tsec += float(c.split('s')[0]) # seconds
	if tomin: return int(tsec/60)
	return int(tsec)



# matrix manipulation
roundmatrix = lambda digits: lambda mat: [[roundna(digits)(x) for x in row] for row in mat]
casteMatrix = lambda typ: lambda mat: [map(typ,row) for row in mat]

# parts of a filename
getLabel = lambda astr: '.'.join(astr.split('/')[-1].split('.')[:-1])
getFilename = lambda astr: astr.split('/')[-1]
getDirectory = lambda astr: astr.split('/')[-2]
getPath = lambda astr: '/'.join(astr.split('/')[:-1])+'/'
getSuffix = lambda astr: astr.split('/')[-1].split('.')[-1]

def unzip(thread):
	"""Reverse of the built-in method zip."""
	try:
		if len(thread):
			return [[thread[i][j] for i in range(len(thread))] for j in range(len(thread[0]))]
		else: return []
	except ValueError: return []


def reFilterFromList(lst,pattern):
	ret = []
	for f in lst:
		if not re.match('^.*%s.*$'%(pattern), f): ret += [f]
	return ret

def findIndices(astr,lst,rc=False):
	p1 = [m.start() for m in re.finditer(astr, lst)]
	if rc:
		import seqDefs as sd
		rastr = sd.rc(astr)
		p2 = [m.start() for m in re.finditer(rastr, lst)]
		pu = p1+p2
		pu.sort()
		return pu
	return p1

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
			if lena > 1:
				if type(B[i][1]) == type([]):
					pairs += [ list(Adct[name])+B[i][1] ]
				else:
					pairs += [ list(Adct[name])+[B[i][1]] ]
			else:
				pairs += [[ Adct[name], B[i][1] ]]
			
			# if type(Adct[name] == type(list)) and lena > 1:
			# 	if type(B[i][1]) == type(list):
			# 		pairs += [ list(Adct[name])+B[i][1] ]
			# 	else:
			# 		pairs += [ list(Adct[name])+[B[i][1]] ]
			# else:
			# 	pairs += [[ list(Adct[name]), B[i][1] ]]

							
			genes += [name]
		elif not inner:
			if lena > 1: pairs += [([nastr for x in xrange(lena)], B[i][1])]
			else: pairs += [(nastr, B[i][1])]
			genes += [name]
	# print 'returning', pairs
	return genes, pairs


def join(A,B,inner=True):
	return threadListPairs(A,B,inner)


def threadSetPairs(lst,inner=True,flat=True):
	if len(lst) == 1: 
		n,d = unzip(lst[0])
		if type(d) == type([]): d = [x for x in d]
		else: d = [[x] for x in d]
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


def joinSet(lst,inner=True,flat=True):
	return threadSetPairs(lst,inner,flat=flat)


def cubicSpline(v, T, Tnew, s=0):
	"""Cubic spline interpolates discrete time curve, returning CS function"""
	
	pairs = [[v[i],T[i]] for i in range(len(v))]
	np = pairedfilterstr(nastr)(pairs)
	y = [float(np[i][0]) for i in range(len(np))]
	t = [float(np[i][1]) for i in range(len(np))]
	
	# y = kernelsmoothlist(y, w=kw) # smooth expression data
	tck = splrep(t,y, k=3, s=0, task=0, per=0)
	
	# reevaluate at high resolution
	dt = max(T)/float(len(Tnew))
	myrange = [i*dt for i in range(len(Tnew))]
	
	smoox = splev(myrange, tck)
	# print myrange
	
	# smoox = np.transpose(splev(Tnew, tck))
	
	# trimx = smoox[release:len(smoox)] # ignore before alpha release
	# trimtimes = newT[release:len(newT)]
	
	return smoox

def integrate(v, dt=1, shift=False):
	"""Integrate expression profile with respect to mean expression"""
	base = 0.0
	if shift: base = abs(min(v))
	vc = [v[i] + base for i in range(len(v))]
	return sum([vc[i]*dt for i in range(len(v)-1)])


def runRscript(scpt, outdir=''):
	import os
	fh = open(outdir+'Rscript.R', 'w')
	print >> fh, scpt
	fh.close()
	
	os.system("R --vanilla < %s > %s" % (outdir+'Rscript.R', outdir+'Rlog.out.R'))
	
	return 1



def strictHistogram(dat, nbins=100, bins=[], categorical=False, yfreq=False, logscale=False, logbase=10, shift=0, dtype=floatna):
	import math
	
	if categorical: dtype = None
	elif logscale: dat = map(lambda x: math.log(x,logbase), dat)
	
	if len(bins) > 0: 
		nbins = len(bins)
		givenbins = True
	else: givenbins = False
	# print 'we have nbins', nbins, len(bins), bins
	
	# need to create a histogram
	filt = filterstr(nastr)(dat)
	if dtype: filt = map(dtype, filt)
	dat = None
	if not len(filt): return [[0,0]]
	
	# mf = minna(filt); xf = maxna(filt)
	# delta = (xf-mf)/float(nbins)
	
	counts = [0 for i in range(nbins)]
	if not givenbins and not categorical: bins = makeRange(min(filt), max(filt), n=nbins)
	elif not givenbins: bins = unique(filt)
	
	# counts = []
	# 	for j in range(int(nbins)): 
	# 		if not givenbins: bins += [mf + delta*j]
	# 		counts.append(0)
	
	# bins = [mf + delta*j for j in range(int(nbins))]
	# counts = [0 for j in range(int(nbins))]
	# print delta, len(bins), len(counts), 'bins', min(bins), max(bins)
	
	# minor bin correction for plotting
	if shift != 0: bins = map(lambda x: x+shift, bins)
	
	# this is like O(N2) too slow
	if not categorical:
		for j in range(len(filt)):
			k = nbins-1
			while k > 0 and float(filt[j]) < bins[k]: k -= 1
			
			# while k < nbins and float(filt[j]) >= bins[k]: k += 1
			# if k == nbins: k -= 1
			
			counts[k] += 1
	else:
		bindct = {}
		for b in bins: bindct[b] = 0
		bk = bindct.keys()
		for j in range(len(filt)): bindct[filt[j]] += 1
		# convert back to a histogram
		bins = []; counts = []
		for key in bindct.keys(): bins += [key]; counts += [bindct[key]]
		
		
		
	# print 'second bins/counts'
	# print zip(bins,counts)
	
	# count correction for last bin count = 0
	# for j in range(nbins):
		# if counts[j] == 0: counts[j] = 1
	
	tmp = []
	intnbins = int(nbins)
	# counts or frequency?
	if yfreq:
		# print 'HIST YFREQ'
		tot = float(sum(counts))
		if logscale: tmp = [[logbase**bins[k],counts[k]/tot] for k in range(intnbins)]
		else: tmp = [[bins[k],counts[k]/tot] for k in range(intnbins)]
	if logscale: tmp = [[logbase**bins[k],counts[k]] for k in range(intnbins)]
	else: tmp = zip(bins,counts)
	
	if categorical: tmp = sorted(tmp)
	
	return tmp


def histogram(dat, nbins=100, bins=[], categorical=False, strict=True, yfreq=False, logscale=False, logbase=10, shift=0, dtype=floatna):
	
	if strict:
		dat2 = strictHistogram(dat, nbins=nbins, bins=bins, yfreq=False, logscale=False, shift=shift, logbase=logbase, dtype=dtype, categorical=categorical)
	else:
		if categorical: dtype = None
		dat = filterstr(nastr)(dat)
		if dtype != None: dat = map(dtype, dat)
		
		dat2 = {}
		if bins != []:
			dat2 = dict( map(lambda x: (x,0), bins) )
			lenbins = len(bins)
			for i in dat:
				try: dat2[i] += 1
				except KeyError: 
					idx = 0
					while idx < lenbins and i>bins[idx]:idx += 1
					if idx == lenbins: idx -= 1
					dat2[bins[idx]] += 1
		else:
			dat2 = {}
			for i in dat:
				try: dat2[i] += 1
				except KeyError: dat2[i] = 1
	
			# fill in gaps / fix bin number
			if bins != []: therange = bins
			elif nbins > 0: therange = range(int(min(dat2.keys())),nbins)
		
			for i in therange: 
				if i not in dat2: dat2[i] = 0
		
			# pileup any bins that are beyond what user wants
			mxtr = therange[len(therange)-1]
			for i in setDifference(dat2.keys(), therange):
				dat2[mxtr] += dat2[i]
				del dat2[i]
		
		# convert dat2 to list of pairs
		dat2 = dat2.items()
	
	if yfreq:
		denom = float(sum(vslice(dat2,1)))
		dat2 = [(bin[0],divna(bin[1],denom)) for bin in dat2]
	
	return dat2



def reHist(hist,bins, freq=False):
	# aggregates an existing histogram in terms of new bins
	
	newhist = dict( map(lambda x: (x,0), bins) )
	nbins = len(bins)
	# which bins corresponds to this key?
	binmap = {}
	for k in vslice(hist,0):
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
		tot = float(sum(vslice(ret,1)))
		# if logscale: tmp = [[logbase**bins[k],counts[k]/tot] for k in range(intnbins)]
		# else: tmp = [[bins[k],counts[k]/tot] for k in range(intnbins)]
		ret = [[bins[k],ret[k][1]/tot] for k in range(len(bins))]
	
	return ret





def groupPairedData(dat):
	
	dct = {}
	for x,y in dat:
		if x in dct: dct[x] += [y]
		else: dct[x] = [y]
	labels = dct.keys(); labels.sort()
	return labels,[dct[x] for x in labels]

def counts(lst):
	counts = {}
	for item in lst:
		try: counts[item] += 1
		except KeyError: counts[item] = 1
	return counts

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


def Fstat(mat, verbose=False, groupByRow=True, ret='f', cutoff=0):
	"""return the F statistic for a 2D matrix. replicates of same group are 2nd dim
	and of length len(mat[i]), and number of groups is 1st dim length len(mat).
	This can be reversed using groupByRow
	"""
	
	# make float
	for i in range(len(mat)):
		mat[i] = map(floatna, mat[i])
		#for j in range(len(mat[i])):
		#	mat[i][j] = floatna(mat[i][j])
	
	if not groupByRow:
		mat = transpose(mat)
	
	try:
		# filter out groups with single replicate
		tup = [[avg(mat[i]),variance(mat[i]),mat[i]] for i in range(len(mat))]
		natup = pairedfilterstr(nastr)(tup) # this will check 1st two entries of tup
		a = len(natup)
		if a < 2: raise TypeError
		
		newmat = [tup[i][2] for i in range(a)]
		
		ns = [len(newmat[i]) for i in range(a)]
		n = sum(ns)
		t = sum([sum(newmat[i]) for i in range(a)])
		gmean = t/float(sum(ns))
		
		# print 'groups', a, 'reps', n/float(a)
		
		avgs = [natup[i][0] for i in range(a)]
		va = [natup[i][1] for i in range(a)]
		
		WGSS = sum([va[i] for i in range(a)])
		BGSS = sum([ns[i]*(avgs[i] - gmean)**2 for i in range(a)])
		
		try:
			F = BGSS/float(WGSS) * (n-a)/float((a-1))
			WGMS = WGSS/float(n-a)
			BGMS = BGSS/float(a-1)
			
			if verbose:
				print '### Fstat ###'
				print 'F', F, 'BGMS', BGMS, 'WGMS', WGMS
				print 'a',a,'ns',ns,'n',n,'t',t
				print 'BGSS',BGSS,'WGSS',WGSS,'n-a',n-a,'a-1',a-1
				print 'averages\n', [avg(newmat[i]) for i in range(a)]
				print 'variance\n', [variance(newmat[i]) for i in range(a)]
		
		except ZeroDivisionError:
			#print 'error', 'BGSS',BGSS,'WGSS',WGSS,'n-a',n-a,'a-1',a-1
			return nastr
		if F >= cutoff:
			if ret == 'f': return F
			elif ret == 'b': return BGMS
			elif ret == 'w': return WGMS
		else:
			return nastr
	except TypeError:
		raise TypeError, "F statistic is NaN."


def trimstr(X, astr):
	"""Remove rows from 2D array X in which an element in the row is astr."""
	X2 = []
	for i in range(len(X)):
		addrow = True
		for j in range(len(X[i])):
			if X[i][j] == astr:
				addrow = False
		if addrow:
			X2.append(X[i])
		#else:
		#	print 'pop row', i, X[i]
	return X2


def hms(fpd):
	if fpd < 60:
		return fpd
	elif fpd < 60**2:
		return "%s:%s" % (int(fpd/60), fpd-int(fpd/60)*60)
	else:
		h = int(fpd/60**2)
		fpd -= h*60**2
		m = int(fpd/60)
		fpd -= m*60
		s = fpd
		return "%s:%s:%s" % (h, m, s)



def deepLookup(name2id, name, kalias):
	"""Normally you would call name2id[name][id] to get the id associated
	with the name, however you might want to do a full search for aliases...
	i.e. sometimes name2id will list name under the heading key instead 
	of name. This def searches the entire dictionary for the value of 
	kname for name that is listed under kalias."""
	
	# Go through all values under key heading
	for tabname in name2id.keys():
		# parse out compound alias fields
		aliases = name2id[tabname][kalias].strip('\"').split('|')
		#print "aliases:", aliases
		if name in aliases:
			#print "aliases for", name, ":", aliases
			# now do a lookup for the kid
			return tabname
	return nastr # if nothing found


def unique(lst, sorted=True):
	"""Returns a subset of lst containing the unique elements."""
	
	if not sorted:
		d = []
		for i in lst:
			if i not in d: d += [i]
		return d
	else:
		d = {}
		for i in lst:
			if i not in d: 
				d[i] = 0
			else: d[i] += 1
		e = d.keys(); e.sort()
		return e

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

def linearize(mat):
	"""similar to flatten, except faster for 2D lists."""
	a = []
	for i in range(len(mat)): a += mat[i]
	return a


def subset(pairs,filtlabels):
	labs = {}
	for i in range(len(filtlabels)):
		labs[filtlabels[i]] = 1
	return filter(lambda x: x[0] in labs, pairs)


#
def zipPairs(*lst):
	if len(lst) <= 1: return lst
	x,y = ut.threadListPairs(lst[0],lst[1])
	for i in range(2,len(lst)):
		x,y = ut.threadListPairs(zip(x,y),lst[i])
	return x,map(lambda x: ut.flatten(x), y)





# now we want to replace d[:,c] with v
def replaceVector(d,c,v):
	def cp(i,x): d[i][c] = x
	
	map(cp, range(len(d)), v)
	return 1

# colvec returns a column vector from a dictionary d, column j
# e.g. colvec(dct)(key)
colvec = lambda d: lambda j: lambda i: float(d[i][j])
select = lambda d: lambda i: d[i]

#select = lambda d: lambda i: d[i]
vcol = lambda j,dat: map(col(j), dat)#map(colvec(d)(j), d.keys())
# vslice = lambda dat, j: map(col(j), dat)
vslice = lambda dat, j: map(lambda row: row[j], dat)
fcol = lambda j: lambda row: float(row[j])
col = lambda j: lambda row: row[j]

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
colslice = multislice


def listRange(lst):
	flst = filterstr(nastr)(lst)
	if not len(flst): return 0
	return float(max(flst) - min(flst))


makeDictionary = lambda mylist,value: dict(map(lambda x: (x,value), mylist))
makeDictFromList = makeDictionary

def makeDict(mylist,value=None):
	return makeDictionary(mylist,value)

makeDictFromPairs = lambda mylist: dict(map(lambda x: (x[0],x[1]), mylist))


def makeRange(xmin,xmax,interval=None,n=None):
	if n == None: n = int(abs(xmax-xmin)/float(interval))
	else: 
		interval = (xmax-xmin)/float(n)
		interval += interval/float(n-1)
	return [xmin+interval*i for i in range(n)]


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






# def spline(x,y,order=4,knots=3, fit='sse'):
# 	"""Fit a smooth spline function of order to the (x,y) data based on knots 
# 	knot points. Default is cubic spline (intercept, 1st, 2nd, 3rd derivatives).
# 	Continuity is preserved up to order-1 derivatives.
# 	"""
# 	
# 	# divide x-range into knots evenly spaced intervals
# 	minx = min(x); maxx = max(x)
# 	delta = (maxx - minx)/float(knots-2)
# 	intervals = [minx]
# 	cknot = minx
# 	for i in range(knots-2):
# 		intervals.append(cknot+delta)
# 	intervals.append(maxx)
# 	
# 	# fit a polynomial to the y data in each subinterval
# 	
# 	
# 

def pivotVectors(R, cols=[]):
	"""Returns a dictionary associating values in cols with their pivot cols.
	if cols not given, key is the index of R from which was seen a pivot col."""
	
	S = array(R)
	P = dict()
	# for each column, except the augmented one
	idx = -1	# the row of the current column's leading 1
	for i in range(len(S[0])-1):
		# following the rules for a rre matrix, we want the columns
		# whose sum is 1, and successive columns must have the leading
		# one in lower rows
		if sum(S[:,i]) == 1:
			cidx = S[:,i].tolist().index(1)
			if cidx > idx:
				if cols == []:
					P[i] = S[:,i]
				else:
					P[cols[i]] = S[:,i]
				idx = cidx
	return P


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






def bedOverlap(set1,set2,minbases, returnProportions=False):
	"""Given two lists of triples (chrom,start,stop), return those that overlap at least minbases"""
	
	# first make a dictionary keyed by chrom
	d1 = {}
	for i,j,k in set1:
		try: d1[i] += [(float(j),float(k))]
		except KeyError: d1[i] = [(float(j),float(k))]
	
	tab = []
	P = [] # overlap proportions per pair of coordinates
	for i,j,k in set2:
		ol = 0
		cands = []
		if i in d1:
			for query in d1[i]:
				ol = overlap((float(j),float(k)),query)
				
				if returnProportions:
					P += [(ol/abs(float(j)-float(k)) + ol/abs(float(query[1])-float(query[0])))/2.]
				
				if ol >= minbases: cands += [(ol,query)]
		cands.sort()
		
		if len(cands):
			keep = cands[-1]
			# find the largest overlapping peak
			# ol = max([overlap((j,k),query) for query in d1[i]])
			entry = [keep[0], i, keep[1][0], keep[1][1], j, k]
			tab += [entry]
	
	if returnProportions: return tab,P
	return tab


def transposeOLD(mat):
	"""Given a 2d array of arrays, will invert the dimensions."""
	if not len(mat): return mat
	newm = []
	for c in range(len(mat[0])):
		newrow = []
		for r in range(len(mat)):
			newrow.append(mat[r][c])
		newm.append(newrow)
	return newm

def transpose(mat):
	"Fast matrix transposition using list comprehensiion. The list map converts tuples to lists."
	return map(list, zip(*mat))


###########################################################################
# BASIC MATHS
###########################################################################
add = lambda x,y: x+y
sub = lambda x,y: x-y
mult = lambda x,y: x*y
div = lambda x,y: x/y

############
# STATISTICS
############
eq = lambda s: lambda i: str(s)!=str(i)

prod = lambda lst: reduce(mul, lst)

def prodna(lst, na=True, nastr=nastr):
	if na: 
		lst = filterstr(nastr)(lst)
		if len(lst): return reduce(mul, lst)
		return nastr
	return reduce(mul, lst)

def sumna(lst, na=True, nastr=nastr):
	if na:
		tmp = filterstr(nastr)(lst)
		return len(tmp)>0 and reduce(add, tmp) or nastr
	return reduce(add, lst)

def avg(v, w=[], na=True, ret=nastr, nastr=nastr, cast=floatna):
	"""When using a weight vector w, this assumes that NAs
	in list v have weight 0.0.
	"""
	if len(v) == 0: return ret
	if len(w) > 0 and len(w) != len(v): return ret
	v2 = v[:]
	if na: v2 = filter(eq(nastr), v)
	
	v2 = map(cast, filter(eq(None), v2))
	if len(v2) == 0: return ret
	elif len(unique(v2))==1 and v2[0] == nastr: return ret
	if len(w):
		return sum([v2[i]*float(w[i]) for i in range(len(v2))])
	else:
		try:
			return sum(v2)/float(len(v2))
		except TypeError:
			raise TypeError, 'utilities.py | avg(): type error on input list: %s...'%(v2[0:10])


def wavg(v1, v2, w=0.5, na=True):
	"""Weighted average of 2 vectors"""
	assert len(v1)==len(v2),\
	 "wavg: input lists must have same length: %s vs %s"%(len(v1), len(v2))
	assert len(v1)>0,\
	 "wavg: empty list."
	
	combo = v1[:]
	pairwavg = lambda w: lambda i1, i2: \
			   (i1 != nastr and i2 != nastr and w*float(i1)+(1.0-w)*float(i2)) or \
			   (i1 != nastr and i2 == nastr and float(i1)) or \
			   (i2 != nastr and i1 == nastr and float(i2)) or \
			   (i1 == nastr and i2 == nastr and nastr)
	miniavg = pairwavg(w)
	# def faltr(x): 
	# 	if x == False: return 0.0
	# 	else: return x
	# 
	# return map(faltr, map(miniavg, v1, v2))
	return map(lambda x: x==True and x or 0., map(miniavg, v1, v2))



def arithmeticMean(lst, ret=nastr, na=True, nastr=nastr, cast=floatna):
	return avg(lst,ret=ret,na=na,nastr=nastr,cast=cast)

def harmonicMean(lst, ret=nastr, na=True, nastr=nastr, cast=floatna, minval=1e-20):
	if len(lst) == 0: return ret
	v = lst
	if na: v = map(cast, filter(eq(nastr), v))
	else: v = map(cast, v)
	
	if len(v) == 0: return ret
	
	# prevent 0 division
	v = map(lambda x: x==0 and minval or x, v)
	return len(v) / sum(map(lambda x: 1/x, v))

def geometricMean(lst, ret=nastr, na=True, nastr=nastr, cast=floatna):
	if len(lst) == 0: return ret
	v = lst
	if na: v = filterstr(nastr)(v)
	v = map(cast, filter(eq(None), v))
	if len(v) == 0: return ret
	return prod(v)**(1/float(len(v)))



def percentile(v=[], k=0.5, na=True, nastr=nastr, exact=False):
	"""
	Return the value of array that corresponds to the kth percentile of the data.
	"""
	if exact: return truePercentile(v, k, na, nastr)
	
	v2 = v[:]
	if na: v2 = filterstr(nastr)(v)
	if len(v2) <= 0: return 0.0
	temp = map(float, v2)
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

# call through percentile function
def truePercentile(v=[], k=0.5, na=True, nastr=nastr):
	"""
	Return the value of array that corresponds to the kth percentile of the data.
	"""
	v2 = v[:]
	if na: v2 = filterstr(nastr)(v)
	if len(v2) <= 0: return 0.0
	temp = map(float, v2)
	temp.sort()
	n = len(temp)
	idx = float(k)*(n+1)
	lidx = math.floor(idx)
	hidx = math.ceil(idx)
	if idx < 1.0:
		return float(temp[0])
	elif idx > n-1:
		return float(temp[n-1])
	else:
		return float(temp[int(idx)])


#
def argPercentile(v=[], k=0.5, na=True, nastr=nastr):
	"""
	Return the value of array that corresponds to the kth percentile of the data.
	"""
	v2 = v[:]
	if na: v2 = filterstr(nastr)(v)
	if len(v2) <= 0: return 0.0
	temp = map(float, v2)
	temp.sort()
	n = len(temp)
	idx = float(k)*(n+1)
	lidx = math.floor(idx)
	hidx = math.ceil(idx)
	if idx < 1.0:
		return 0
	elif idx > n-1:
		return n-1
	else:
		return int(idx)



def quantiles(v=[], na=True, nastr=nastr, Q=[.05,.25,.5,.75,.95]):
	v2 = v
	if na: v2 = filterstr(nastr)(v)
	return dict( map(lambda q: (q, percentile(v2,q,na=False)), Q) )


def pvalue(lst, true, nastr=nastr):
	"""The p-value is the number of times a random value <= true value
	divided by the number of trials
	"""
	nullvals = filterstr(nastr)(lst)
	if len(nullvals) == 0: return nastr
	nullvals.sort()
	
	count = 0
	CONT = True
	while CONT and count < len(nullvals):
		if nullvals[count] >= true: CONT = False
		else: count += 1
		
	
	return count/float(len(nullvals))


def pvalue_BAD(v, d, na=True, nastr=nastr, exact=True, limit=.0001, maxiter=1000, delta=.01):
	"""Returns the pvalue of d in empirical distribution v.
	THIS DOESN'T WORK VERY WELL.
	"""
	
	if d == nastr: return nastr
	diff = 1e18
	k = 0.0
	bestk = k
	count = 0
	while k < 1. and diff > limit and count < maxiter:
		pct = percentile(v, k, exact=exact)
		diff = absna(d - pct)
		bestk = k
		k += delta
		# print k, bestk, diff, count, pct
		count += 1
	
	
	
	
	# ppct = 0.
	# delta = nastr
	# k = .5
	# count = 0
	# while (delta > limit and count < 1000):
	# 	pct = percentile(v,k, exact=True)
	# 	if d <= pct:
	# 		k += limit
	# 	else:
	# 		k -= limit
	# 	delta = pct - ppct
	# 	count += 1
		
	return bestk


def median(v, na=True, nastr=nastr):
	return percentile(v, 0.5, na, nastr=nastr)


def mode(v, na=True, nastr=nastr):
	v2 = v[:]
	if na: v2 = filterstr(nastr)(v)
	if len(v2) <= 0: return 0.0
	
	# mode is most frequent value
	entries = {}
	for i in v2:
		try: entries[i] += 1
		except KeyError: entries[i] = 1
	keys = v2[:]; keys.sort()
	values = [entries[k] for k in keys]
	pairs = zip(values,keys); pairs.sort(reverse=True)
	return pairs[0][1]

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

def MAD(v,na=True,ret=nastr):
	"""Median absolute deviation"""
	v2 = v
	if na: v2 = filterstr(nastr)(v)
	gmed = median(v2)
	return median([abs(vi - gmed) for vi in v2])

def sd(v, na=True, ret=nastr):
	return stdev(v,na,ret)

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

def zfunc(X,center=avg,scale=stdev):
	m = center(X)
	sd = scale(X)
	return map(lambda xi: divna(subna(xi,m),sd), X)

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

def cov(x,y=[], na=True, nastr=nastr, ret=nastr):
	return covariance(x,y, na, nastr, ret)

def correlation(x,y=[],na=True,ret=nastr, method='pearson', cast=floatna):
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
		if method == 'spearman':
			x = rank(x,ties='split')
			y = rank(y,ties='split')
	except IndexError:
		print 'RANKING PROBLEM in utilities.py::correlation(method="spearman")'
		# print 'X', x
		# print 'Y', y
		return ret

	try:
		R = divna(covariance(x,y,na=False), floatna(multna(stdev(x,na=False), stdev(y,na=False))))
		return R
	except ZeroDivisionError: return ret
	except IndexError: return ret

def cor(x,y=[],na=True,ret=nastr, method='pearson', cast=floatna):
	return correlation(x,y,na,ret, method=method, cast=cast)
def abscor(x,y=[],na=True,ret=nastr, method='pearson', cast=floatna):
	return absna(correlation(x,y,na,ret, method=method, cast=cast))
def negcor(x,y=[],na=True,ret=nastr, method='pearson', cast=floatna):
	return subna(1,correlation(x,y,na,ret, method=method, cast=cast))

def upperTri(mat, includeIden=True):
	mat2 = [[nastr for j in range(len(mat))] for i in range(len(mat))]
	if not includeIden:
		for i in range(len(mat)):
			for j in range(i): mat2[i][j] = mat[i][j]
	else:
		for i in range(len(mat)):
			for j in range(i+1): mat2[i][j] = mat[i][j]
	return mat2


def offDiag(mat):
	return upperTri(mat, includeIden=False)

def flatOffDiag(mat):
	res = []
	for i in range(len(mat)):
		for j in range(i):
			res += [mat[i][j]]
	return res


def matrixCorrelation(f1,f2,type='file',nperms=1000,cormethod='spearman',mantelfile='/Users/simola/bin/mantel.test.R',outdir='',rows=[],header=[]):
	"""Given 2 data matrices, compute each of their correlation matrices and then compare the
	correlation matrices using a Mantel test.
	"""
	import sio as io
	if type == 'data':
		# save them to file
		tmp1 = outdir+'data1.txt'; tmp2 = outdir+'data2.txt'
		io.printMatrix(f1,rows,header,file=tmp1)
		io.printMatrix(f2,rows,header,file=tmp2)
		f1 = tmp1; f2 = tmp2
	
	# load the mantel test code
	fh = open(mantelfile)
	scpt = ''.join(fh.readlines())
	fh.close()
	
	# using the mantel matrix correlation coefficient
	scpt += """
	X = read.delim("%s")
	Y = read.delim("%s")
	X = X[,2:length(X)]
	Y = Y[,2:length(Y)]
	
	cmx = cor(t(X), method="%s")
	cmy = cor(t(Y), method="%s")
	
	mant = mantel.fcn(cmx, cmy, nperm=%d)
	R = mant$R
	p = mant$UpperP
	
	write(c(R,p), file="%s")
	""" % (f1, f2, cormethod, cormethod, nperms, outdir+'rcorout.txt')
	fh = open(outdir+'rcorscript.R', 'w')
	print >> fh, scpt
	fh.close()
	
	os.system('R --vanilla < %s > %s' % (outdir+'rcorscript.R', outdir+'rcorout.out.R'))
	
	d,r,c = io.readMatrix(outdir+'rcorout.txt', header=False)
	# print 'R/pval:', r
	vals = map(float, r[0].split(' '))
	R, P = vals[0], vals[1]
	return {'R':R, 'P':P}


def rescaleVector(v, scale=1, scalefunc=avg):
	"""Take a length |v| vector and collapse into a length w vector"""
	
	return [scalefunc(v[curr:curr+scale]) for curr in makeRange(0, len(v), scale)]

rescaleVectorLambda = lambda scale: lambda v: rescaleVector(v,scale)


def jackknifeCorrelation(x,y=[],rem=0,na=True,ret=nastr,method='pearson', func=min):
	if not len(y):
		# x is threaded as pairs
		x2 = [x[i][0] for i in range(len(x))]
		y2 = [x[i][1] for i in range(len(x))]
		x = x2; y = y2
	if not len(x): return ret
	
	x = map(floatna, x)
	y = map(floatna, y)
	
	# correlate all single observations
	pairs = zip(x,y)
	
	jkpairs = []
	lst = range(len(pairs))
	if rem > 0:
		lst = random.sample(range(len(pairs)), rem)
	for i in lst:
		# remove ith observation
		# print 'control len orig pairs', len(pairs)
		tmp = pairs[:]
		tmp.pop(i)
		# print '  len tmp', len(tmp)
		jkpairs += [ tmp ]
	
	# if rem > 0:
	# 	# choose rem randomly
	# 	jkpairs = random.sample(jkpairs,rem) 
	
	# print 'computing'
	docor = lambda pair: correlation(pair,na=na,ret=ret,method=method)
	jkcors = map(docor, jkpairs)
	# print 'jkcors', jkcors
	# print ' with min', min(jkcors)
	return func(jkcors)
	
	


def fisherTransform(p,n):
	import maths
	# n = len(x)
	# print p, n, (1.+p)/(1.-p), ln((1.+p)/(1.-p))
	Fr = 0.5*ln((1.+p)/(1.-p))
	return math.sqrt((n-3)/1.06)*Fr

def fisherMethod(plst,minp=1e-10):
	"""Fisher's method for combining p-values from independent tests bearing on the same H0"""
	
	#x2 =  -2 sum (ln pi)
	plstF = map(lambda x: x!=0 and x or minp, plst)
	plstF = filter(lambda x: x!=nastr, plstF)
	try:
		X2 = multna(-2,sum(map(ln, plstF)))
	except TypeError:
		print 'hey', plst
		sys.exit()
	df = 2*len(plstF)
	P = X2!=nastr and ChiSquarePDF(X2, dof=df) or nastr
	return {'X2':X2, 'P':P, 'df':df}
	


def corWithPvalue(x,y=[],na=True,method='pearson', pmethod='Z', tails='2', p='true', scale=None, smooth=None, ret=nastr, invert=False, invertpp=False, jacknife=False, jackpct=0.1, jackiter=100):
	"""Returns a dictionary with correlation 'r', 'p-value', 'n', and 't'."""
	
	if not len(y):
		# x is threaded as pairs
		x2 = [x[i][0] for i in range(len(x))]
		y2 = [x[i][1] for i in range(len(x))]
		x = x2; y = y2
	if not len(x): return {'r':nastr, 't':nastr, 'P':nastr, 'n':0}
	x,y = filterStrFromLists([x,y])
	
	if not len(x) or not len(y): return {'r':nastr, 't':nastr, 'P':nastr, 'n':0}
	
	if scale:
		x = rescaleVector(x,scale)
		y = rescaleVector(y,scale)
	
	if smooth:
		x = kernelsmooth(x,w=smooth)
		y = kernelsmooth(y,w=smooth)
	

	r = nastr
	if jacknife:
		pairs = zip(x,y)
		lst = range(len(pairs))
		allR = []
		
		# simple jacknife computes correlation removing all single values and takes min R
		# for i in range(len(pairs)):
		# 	# remove ith observation
		# 	tmp = pairs[:]
		# 	tmp.pop(i)
		# 	allR += [correlation(tmp,na=na,method=method)]
		
		# this removes X% of data points sampled Y times and returns minimum
		nr = int(jackpct * len(pairs))
		import random
		for i in range(jackiter):
			allR += [correlation(random.sample(pairs,len(pairs)-nr),na=na,method=method.lower())]

		r = minna(allR)
	else:
		r = correlation(x,y,na=na,method=method.lower())
	
	if invert: r = subna(1,r)
	# if collapsep and r > 0: r = collapsep
	
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
		
		except Exception, e:
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


def percentileForValues(X, vals, na=False):
	V = vals
	if type(vals) != type([]) or type(vals) != type(()): V = list(vals)
	if na: X = filterstr(nastr)(X)
	if not len(V): return [nastr]
	R = rank(X+V, na=False)
	D = X+filterstr(nastr)(V)
	return [V[i]==nastr and nastr or D.index(V[i])/float(len(D)) for i in range(len(V))]



def rankPairs(x,y=[], na=True, ties='split'):
	if not len(y):
		# x is threaded as pairs
		x2 = [x[i][0] for i in range(len(x))]
		y2 = [x[i][1] for i in range(len(x))]
		x = x2
		y = y2
	
	if na:
		pairs = [[x[i],y[i]] for i in range(len(x))]
		pairs = pairedfilterstr(nastr)(pairs)
		x = [pairs[i][0] for i in range(len(pairs))]
		y = [pairs[i][1] for i in range(len(pairs))]
	
	x = rank(x,na=na,ties=ties)
	y = rank(y,na=na,ties=ties)
	
	return [[x[i],y[i]] for i in range(len(x))]

def rankCorrelation(x,y=[],na=True,ret=nastr):
	return correlation(rankPairs(x,y), method='pearson',na=na)



# mannWhitneyTest = rankSumTest
# wilcoxonTest = signedRankTest


# confirm this is equivalent to rank(split=True)
def rankList(X, integer=False, sort=False, reverse=False):
	"""Returns a list of ranks for input list. Breaks ties."""
	
	X2 = X
	if sort:
		X2 = X[:]
		X2.sort(reverse=reverse)
	
	R = []
	c = 1 # starting rank
	i = 0
	
	while i < len(X2):
		# print 'outside', i, len(X)
		equiv = []
		j = i
		while j < len(X2) and X2[j] == X2[i]:
			# print '  inside', j, len(X)
			equiv += [c]
			c += 1
			j += 1
			
		# assign equal ranks
		# R += [avg(equiv) for e in equiv]
		if not integer:
			R += [avg(equiv) for e in equiv]
		else:
			R += range(equiv[0],equiv[0]+len(equiv))
		i = len(R)
	
	return R


def rankSumTest(X,Y=[], tails=2, verbose=False, minsamplesize=20, minsingle=3):
	# mann-whitney U test for independent samples
	import maths
	
	blanket = {'n':len(X),'m':len(X), 'W+':nastr, 'W-':nastr, 'W':nastr, 'U':nastr, 'P':nastr, 'z':nastr}
	# blanket = {'n':len(X),'m':len(Y), 'n1':n1, 'n2':n2, 'U1':U1, 'U2':U2, 'U':U, 'z':nastr, 'P':nastr}
	
	if not len(Y): Y = vslice(X,1)
	
	# 1. rank all observations
	# compute differences
	Z = [(xi,'X') for xi in X] + [(yi,'Y') for yi in Y]
	Z = pairedfilterstr(nastr)(Z)
	Z.sort()
	# print 'Z', Z
	# print 'OK, ranking'
	R = rankList(vslice(Z,0))
	# print 'done ranking'
	# count up ties
	tiegroups = []
	cur = []
	for i in range(len(R)):
		if i == 0: cur = [R[i]]
		else:
			# print 'test', i, R[i], cur
			if R[i] == cur[0]: cur += [R[i]]
			else:
				if len(cur)>1: tiegroups += [len(cur)]
				cur = [R[i]]
	tiegroups += [len(cur)]
	
	if verbose:
		print 'Z', len(Z), Z
		print 'R', len(R), R
		print 'ties', tiegroups
	
	
	# 2. extract ranks from sample 1
	R1 = [R[i] for i in range(len(Z)) if Z[i][1]=='X']
	n1 = len(R1)
	R1 = sum(R1)
	U1 = R1 - n1*(n1+1)/2.
	
	R2 = [R[i] for i in range(len(Z)) if Z[i][1]=='Y']
	n2 = len(R2)
	R2 = sum(R2)
	U2 = R2 - n2*(n2+1)/2.
	
	# print 'R1',R1, 'R2', R2
	# print 'U1', U1, 'U2', U2
	
	U = minna([U1,U2])
	
	blanket = {'n':len(X),'m':len(Y), 'U1':U1, 'U2':U2, 'U':U, 'z':nastr, 'P':nastr}
	
	# if min(n1,n2) < 3 or n1*n2 < 15: # should be ~20
	# minsingle should be 3
	# minsamplesize should be ~20
	if min(n1,n2) < minsingle or n1*n2 < minsamplesize:
		# print >> sys.stderr, 'Not yet implemented'
		return blanket
	else:
		# Gaussian approximation
		
		mu = n1*n2/2.
		N = n1+n2
		sd = math.sqrt(n1*n2*(N+1)/12.)
		# print 'uncsd', sd
		if len(tiegroups) > 0:
			# special correction
			# print 'test', sum([(tiegroups[j]**3 - tiegroups[j])/12. for j in range(len(tiegroups))])
			# print 'test2', ( (N**3 - N)/12. - sum([(tiegroups[j]**3 - tiegroups[j])/12. for j in range(len(tiegroups))]) )
			sd = math.sqrt(n1*n2/float(N*(N-1)) * ( (N**3 - N)/12. - sum([(tiegroups[j]**3 - tiegroups[j])/12. for j in range(len(tiegroups))]) ))
		
		
		# Compute Z-score
		Z = divna(U-mu,sd)
		P = nastr
		if Z != nastr:
			P = None
			if tails == 2:
				P = maths.GaussianCDF(Z,0,1)
				P = 2*min(P,1-P)
			elif tails == 1 or tails == -1:
				P = maths.GaussianCDF(Z,0,1)
				P = min(P,1-P)
			elif tails.lower() == 'cdf':
				P = maths.GaussianCDF(Z,0,1)
		
			if P == 0: P = 1e-16
			elif P == 1: P = 1-1e-16
		
		
		blanket['z'] = Z
		blanket['P'] = P
	
	
	if verbose:
		print 'U1=%s, U2=%s, U=%s, Z=%s, P=%s'%(U1,U2,U,Z,P)
	
	
	return blanket


Utest = rankSumTest
mannWhitneyTest = rankSumTest

def kruskalWallisTest():
	"""Need to implement!"""
	return 0

# need to include 1 sample testing
def signedRankTest(X, Y=[], na=True, nastr=nastr, lenx=None, leny=None, tails=2, verbose=False):
	"""Wilcoxon two-sample (paired) signed rank W test. Will eventually include 1-sample testing."""
	import maths
	
	Z = None
	
	if not lenx: lenx = len(X)
	if not leny: leny = len(Y)
	
	blanket = {'n':lenx,'m':leny, 'W+':nastr, 'W-':nastr, 'W':nastr, 'P':nastr, 'z':nastr}
	
	# compute differences
	if len(Y): 
		# assert len(X)==len(Y), 'signedRankTest(): X and Y lengths must be equal.'
		if len(X) != len(Y): return blanket
		Z = map(subna,X,Y)
	else:
		Z = map(subna,*X)
		Y = vslice(X,1)
	
	if na:
		Z = filterstr(nastr)(Z)
		lenx = leny = len(Z)
		blanket = {'n':lenx,'m':leny, 'W+':nastr, 'W-':nastr, 'W':nastr, 'P':nastr, 'z':nastr}
		

	# trim 0 differences
	Zm = filter(lambda x: x!=0, Z)
	# compute sample size
	n = len(Z)
	m = len(Zm)
	
	blanket = {'n':n,'m':m, 'W+':nastr, 'W-':nastr, 'W':nastr, 'P':nastr, 'z':nastr}
	if len(Zm)==0: return blanket
	
	# compute absolute values
	Za = [abs(zi) for zi in Zm]
	# compute sign
	Zp = map(lambda x: x>0 and 1. or 0., Zm)
	Zn = map(lambda x: x<0 and 1. or 0., Zm)
	
	Z = zip(Za,Zm,Zp,Zn)
	# sort
	Z.sort() # sort from smallest to largest using |differences|
	Za = vslice(Z,0)
	Zm = vslice(Z,1)
	Zp = vslice(Z,2)
	Zn = vslice(Z,3)
	# assign ranks (assumes sorted)
	Zr = rankList(Z)
	
	# compute wilcoxon statistic
	Wp = [Zp[i]*Zr[i] for i in xrange(len(Zp))]
	Wn = [Zn[i]*Zr[i] for i in xrange(len(Zp))]
	
	sWp = sum(Wp)
	sWn = sum(Wn)
	
	W = min(sWp,sWn)
	# Wx = max(sWp,sWn)
	
	if n < 5: # should be 10
		print >> sys.stderr, 'Not yet implemented'
		return blanket
	else:
		# Gaussian approximation to the binomial
		# see Biometry Sokal and Rholf
		
		# Compute Z-score
		Z = (W - 0.5 - m*(m+1)/4.) / math.sqrt(m*(m+0.5)*(m+1)/12.) 
		
		P = None
		if tails == 2:
			P = maths.GaussianCDF(Z,0,1)
			P = 2*min(P,1-P)
		elif tails == 1 or tails == -1:
			P = maths.GaussianCDF(Z,0,1)
			P = min(P,1-P)
		elif tails.lower() == 'cdf':
			P = maths.GaussianCDF(Z,0,1)
		
		if verbose: 
			import sio
			sio.printFormattedTable(zip(Zm,Za,Zp,Zn,Wp,Wn),['Zm','Za','Zp','Zn','Wp','Wn'])
			print 's=%s, z=%s, P=%s'%(s,Z,P)
		
	return {'n':n, 'm':m, 'W':W, 'W+':sWp,'W-':sWn, 'P':P, 'z':Z}


wilcoxonTest = signedRankTest

def ZdistPDF(Z, tails=2):
	if Z == nastr: return nastr
	import maths
	P = None
	if tails == 2:
		P = maths.GaussianCDF(Z,0,1)
		P = 2*min(P,1-P)
	elif tails == 1 or tails == -1:
		P = maths.GaussianCDF(Z,0,1)
		P = min(P,1-P)
	elif tails.lower() == 'cdf':
		P = maths.GaussianCDF(Z,0,1)
	return P



def getR2(x,xhat,adj=False,coefs=0,na=True, nastr=nastr):
	"""R2 coefficient of determination between two vectors"""
	
	assert(len(x)==len(xhat)), \
	 "Vectors must have same length."
	if adj==True:
		assert(coefs > 0),\
	 	 "Number of coefficients estimated must be provided for adjusted R2."
	y = x[:]
	yhat = xhat[:]
	if na:
		mat = [[x[i], xhat[i]] for i in range(len(x))]
		mat = pairedfilterstr(nastr)(mat)
		# recover Y, still in original input order
		y = map(float, [mat[i][0] for i in range(len(mat))])
		# recover X, still in original input order
		yhat = map(float, [mat[i][1] for i in range(len(mat))])
	else:
		y = map(float, y)
		yhat = map(float, yhat)
	
	n = len(y)
	avgy = avg(y)
	sse = float(sum(map(euclid2d, y, yhat)))
	sst = float(sum(map(euclid2d, y, [avgy]*n)))
	
	if adj: 
		return 1 - (sse*(n-1)/float(sst*(n-coefs-1)))
	elif sst > 0:
		return 1 - (sse/sst)
	else:
		return 0


nchoosek = lambda n,k: factorial(n)/float(factorial(k)*factorial(n-k))
# nchoosek = lambda n,k: factit(n)/float(factit(k)*factit(n-k))
binomialCoefficient = nchoosek
bincoef = nchoosek

def binomial(n,x,p, approx=False):
	"""Returns the binomial probability of x successes in n trials"""
	
	if approx:
		# Normal(np, np(1-p))
		import maths
		return maths.GaussianPDF(x, n*p, n*p*(1-p))
	
	c = factorial(n)/float(factorial(x)*factorial(n-x))
	return c*(p**x)*((1-p)**(n-x))


def binomialPDF(n,k,p,approx=False, lower=True):
	if approx:
		import maths
		return maths.GaussianCDF(k, n*p, n*p*(1-p))
	tot = 0
	if lower:
		for i in range(0,k+1):
			c = factorial(n)/(factorial(i)*factorial(n-i))
			tot += c * p**i * (1.0-p)**(n-i)
		return tot
	else:
		for i in range(k,n+1):
			c = factorial(n)/(factorial(i)*factorial(n-i))
			tot += c * p**i * (1.0-p)**(n-i)
		return tot



def multinomial(ks,ps):
	return factorial(sum(ks))/prod(map(factorial, ks)) * prod(map(lambda x,e: x**e, ps, ks))


def factorialO(x):
	if x <= 1: return 1
	else: return x*factorial(x-1)

def factit(x):
	"""Iterative implementation of factorial function."""
	
	if x <= 1: return 1
	
	a = [0 for i in range(x+1)]
	a[0] = 1; a[1] = 1
	for i in range(2,x+1): a[i] = i*a[i-1]
	return a[x]

fact = lambda x: factit(x)
# def fact(x): return factit(x)

factorial = lambda x: factit(x)

# def factorial(n):
# 	import maths
# 	return maths.factit(n)
# 

poisson = lambda l,k: l**k * math.e**(-l)/fact(k)

def poissonPDF(l,k,tol=1e-10):
	k_ = k+1
	p = poisson(l,k)
	last = p
	while abs(p-last) < tol:
		last = poisson(l,k_)
		p += last
		k_ += 1
	return p


# def binomialCoefficient(n,k):
# 	return factorial(n)/(factorial(k)*factorial(n-k))
# 

binomialCoefficient = lambda n,k: fact(n)/(fact(k)*fact(n-k))
nchoosek = lambda n,k: fact(n)/(fact(k)*fact(n-k))

def slowFibonacci(x):
	if x == 0: return 0
	elif x == 1: return 1
	else: return fibonacci(x-1) + fibonacci(x-2)

def fibonacci(x):
	fibs = [0 for i in range(x+1)]
	if x > 0: fibs[1] = 1
	for i in range(2,x+1): fibs[i] = fibs[i-1] + fibs[i-2]
	return fibs[x]

def dictAvg(dctlist, key):
	"""take the average of a particular key for a list of dictionaries, given the
	list of dictionaries and a key"""
	values = []
	# extract value from each dict using key
	for d in dctlist:
		values.append(d[key])
	# average the values
	return avg(values)

def dictSum(dctlist, key):
	"""take the sum of a particular key for a list of dictionaries, given the
	list of dictionaries and a key"""
	values = []
	# extract value from each dict using key
	for d in dctlist:
		values.append(d[key])
	# sum the values
	return sum(values)

def dictMin(dctlist, key):
	values = []
	# extract value from each dict using key
	for d in dctlist:
		values.append(d[key])
	# sum the values
	return min(values)

def dictMax(dctlist, key):
	values = []
	# extract value from each dict using key
	for d in dctlist:
		values.append(d[key])
	# sum the values
	return max(values)

def listAvg(lst, idx):
	""""""
	values = []
	# extract value from each list using idx
	for i in range(len(lst)):
		values.append(lst[idx])
	# average the values
	return avg(values)


def makeCombinations(nentries, nthresholds, start_idx=0, iterator=False):
	import xpermutations as xp
	A = range(start_idx, nentries+start_idx)
	 
	if not iterator:
		B = []
		for i in xp.xuniquecombinations(A, nthresholds):
			B.append(i)
		return B
	else:
		return xp.xuniquecombinations(A, nthresholds)

def combination(lst,nitems):
	import xpermutations as xp
	return list(xp.xuniquecombinations(lst,nitems))


def permutations(L):
    # needs xp.permutations?
    if len(L) <= 1:
        yield L
    else:
        a = [L.pop(0)]
        for p in permutations(L):
            for i in range(len(p)+1):
                yield p[:i] + a + p[i:]

def randomizeRow(lst, withReplacement=True):
	"""Randomizes a list with/without replacement"""
	
	import random
	
	new = []
	usedIndices = []
	
	i = len(lst)
	while i > 0:
		idx = random.randrange(len(lst))
		if not withReplacement:
			while idx in usedIndices: idx = random.randrange(len(lst))
		new += [ lst[idx] ]
		i -= 1
		usedIndices += [idx]
	return new

def permute(lst):
	"""Permutes a list"""
	return randomizeRow(lst, withReplacement=False)


# DISTANCE MEASURES AND GEOMETRY
# ------------------------------
euclid2d = lambda x,y: (x-y)*(x-y)
def dot(x,y):
	assert(len(x)==len(y)),\
	  'vectors must have same length.'
	return sum([float(x[i])*float(y[i]) for i in range(len(x))])

def norm(x, type='euclidean'):
	"""Return magnitude of a vector"""
	if type == 'euclidean':
		return sum( map(lambda y: y**2, x) )**.5
		# return math.sqrt(sum( [x[i]**2 for i in range(len(x))] ))
	else: 
		print 'no other norm types implemented yet.'
		return 0.0

def distance(x, y=None,v=None,CI=None, method='euclidean', na=True):
	# return distance between two vectors x and y
	if not y:
		y = vslice(x,1)
		x = vslice(x,0)
	
	assert(len(x)==len(y)),\
	  'input vectors must match in length.'
	if na: x,y = filterStrFromLists([x,y])
	x = map(float,x)
	y = map(float,y)
	
	if method == 'euclidean':
		return math.sqrt(sum([(x[i]-y[i])**2 for i in range(len(x))]))
	elif method == 'mahalanobis':
		import numpy as np
		if CI == None:
		# if type(CI) != np.core.defmatrix.matrix:
			sys.exit("must provide a covariance matrix")
		CI = np.matrix(CI) # need to be a matrix
		md = np.matrix(x).T - np.matrix(y).T
		# print 'dim of diff', md.shape#, md
		# print 'shapes:', CI.shape, md.shape
		return np.sqrt(md.T*CI*md)[0,0]
	elif method == 'mahalanobis2':
		if type(v) != float: v = variance(x)
		return math.sqrt(sum([(x[i] - y[i])**2/v for i in range(len(x))]))

	elif method == 'hamming':
		return hamdist(x,y)
mahalanobis1 = lambda V: lambda x,y: distance(x,y,CI=V,method='mahalanobis')

def d(x, y=None,v=None,CI=None, method='euclidean', na=True):
	return distance(x,y,v,CI,method,na)

def eucdist(s1,s2):
	assert len(s1) == len(s2)
	return sum( (ch1 - ch2)**2 for ch1,ch2 in zip(s1,s2) )**.5

def eucdistna(s1,s2):
	# tackles missing data removal
	assert len(s1) == len(s2)
	pairs = pairedfilterstr(nastr)(zip(s1,s2))
	return sum( (ch1 - ch2)**2 for ch1,ch2 in pairs )**.5


def eucdistlist(X):
	return [eucdist(X[i],X[j]) for i,j in makeCombinations(len(X),2)]


def eucdistavg(X):
	E = [eucdist(X[i],X[j]) for i,j in makeCombinations(len(X),2)]
	return sum(E)/float(len(E))


def difflist(X):
	# report the list of differences between all pairs of vectors in X
	# where the vector difference is defined as the maximum absolute deviation between paired entries
	return [max(map(lambda x,y: abs(x-y), X[i],X[j])) for i,j in makeCombinations(len(X),2)]

diffpairs = lambda X: [subna(X[i],X[j]) for i,j in makeCombinations(len(X),2)]

# CV with low sample correction
CVc = lambda X: multna1(addna(1,divna(1,(4*lenna(X))))) (divna(stdev(X), avg(X)))
	
# 
def groupDiff(X,Y, op=subna, agglom=avg):
	# compute distance between all pairs then take unweighted average
	D = []
	for i in X:
		for j in Y:
			try:
				D += [op(i,j)]
			except TypeError:
				D += [op([i,j])]
	return agglom(D)
	
upgma = lambda x,y: groupDiff(x,y, op=subna, agglom=avg)


def mindist(X,Y):
	"""Return minimum distance between two pairs"""
	
	return min(map(abs,[X[0]-Y[0], X[1]-Y[0],X[0]-Y[1], X[1]-Y[1]]))
	

def hamdist(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def cosangle(x, y, coord='radians'):
	try:
		a = dot(x,y)/(norm(x)*norm(y))
		if coord == 'radians': return a
		else: return math.degrees(a)
	except ZeroDivisionError: return nastr

def cosanglena(x, y, coord='radians'):
	x,y = unzip(pairedfilterstr(nastr)(zip(x,y)))
	a = dot(x,y)/(norm(x)*norm(y))
	if coord == 'radians': return a
	else: return math.degrees(a)


rmse = lambda x,y: math.sqrt( sum(map(lambda a,b: (a-b)*(a-b), x,y))/float(len(x)) )

def angle(x, y, coord='radians'):
	a = nastr
	try:
		a = math.acos(dot(x,y)/(norm(x)*norm(y)))
		if coord == 'radians': return a
		else: return math.degrees(a)
	except ZeroDivisionError: return a

def anglena(x, y, coord='radians'):
	a = nastr
	x,y = unzip(pairedfilterstr(nastr)(zip(x,y)))
	try:
		a = math.acos(dot(x,y)/(norm(x)*norm(y)))
		if coord == 'radians': return a
		else: return math.degrees(a)
	except ZeroDivisionError: return a

def angulardistance(x,y):
	theta = anglena(x,y,coord='radians')
	return 1 - theta/math.pi

# def absangulardistance(x,y):
# 	theta = anglena(x,y,coord='radians')
# 	return 1 - 2*theta/math.pi

def reflectedAngle(x, y, coord='degrees'):
	ang = angle(x,y,coord=coord)
	if ang > 90.0 and ang != nastr: ang = abs(180.0 - ang)
	return ang

def sign(x):
	if x < 0: return -1
	elif x == 0: return 0
	else: return 1


def anglePairs(A, B, coord='degrees', reflect=True):
	"""
	Returns a list of angles generated as all pairs of columns 
	of matrices A and B, forming a background distribution.
	"""
	At = transpose(A)
	Bt = transpose(B)
	
	# get index pairs
	nx = len(At); ny = len(Bt)
	# print 'nx', nx, 'ny', ny, 'total', nx+ny
	# pairs = makeCombinations(nx,2)
	
	angdist = []
	if reflect:
		for i in range(nx):
			for j in range(ny):
		
		# for i,j in pairs:
				# print 'i,j', i,j
				x = At[i]; y = Bt[j]
				angdist += [ reflectedAngle(x,y, coord=coord) ]
	else:
		for i in range(nx):
			for j in range(ny):
				x = At[i]; y = Bt[j]
				angdist += [ angle(x,y, coord=coord) ]
	return angdist



#
def distanceMatrix(mat,func=eucdistna,baseval=0.):
	"""Computes a distance matrix between column pairs using func"""
	mt = transpose(mat)
	
	dmat = [[baseval for i in range(len(mt))] for j in range(len(mt))]
	# print 'size of distance matrix to make', len(dmat), len(dmat[0])
	for i,j in makeCombinations(len(mt),2):
		# d = func(mt[i],mt[j])
		# if r == ut.nastr: r = 0
		dmat[i][j] = dmat[j][i] = func(mt[i],mt[j])
	
	return dmat



def areaTriangleLen(a, b, c):
	"""Heron's method to compute the area of a triange given 3 side lengths"""
	s = (a+b+c)/2.0
	return math.sqrt(s*(s-a)*(s-b)*(s-c))

def areaTriangleCoords(a,b,c):
	x = norm(a-b); y = norm(b-c); z = norm(c-a)
	return areaTriangleLen(x,y,z)

def signedAreaTriangle(a, b, c):
	"""Returned the signed area of a triangle given 3 2D vertex pairs"""
	x = [a[0],b[0],c[0]]
	y = [a[1],b[1],c[1]]
	# x = [temp[0,t+1], temp[0,t], cent[0]]
	# y = [temp[1,t+1], temp[1,t], cent[1]]
	return (-x[1]*y[0]+x[2]*y[0]+x[0]*y[1]-x[2]*y[1]-x[0]*y[2]+x[1]*y[2])/2.0


def lowess(x, y, F=1/3., NSTEPS=1, tempdir=''):
	"""
	Calls the lowess routine in R via command line.
	"""
	import sio as io
	
	io.printMatrix(enlist(x), file=tempdir+'xdata.txt')
	io.printMatrix(enlist(y), file=tempdir+'ydata.txt')
	pipe = open(tempdir+'lowess_script.R', 'w')
	print >> pipe, 'F = %f' % (F)
	print >> pipe, 'nsteps = %d' % (NSTEPS)
	print >> pipe, 'x = read.delim("%s", sep="\t", header=FALSE)' % (tempdir+'xdata.txt')
	print >> pipe, 'y = read.delim("%s", sep="\t", header=FALSE)' % (tempdir+'ydata.txt')
	print >> pipe, 'x = x[,1]'
	print >> pipe, 'y = y[,1]'
	print >> pipe, 'res = lowess(x=x, y=y, f=F, iter=nsteps)'
	print >> pipe, 'write.table(res$y, file="%s", sep="\t", row.names=FALSE, col.names=FALSE)' % (tempdir+'yfit.txt')
	pipe.close()
	
	os.system('R --vanilla < %s > %s' % (tempdir+'lowess_script.R', tempdir+'lowess.out.R'))
	
	crap,Yfit,c = io.readMatrix(tempdir+'yfit.txt', header=False)
	
	Yfit = map(floatna, Yfit)
	return Yfit

def SVD(X, outdir='', rmtmp=True, LINPACK=False, returnas='numpy', flattenS=True, reuse=False):
	import sio as io, numpy as np
	if not outdir: outdir = os.getcwd()+'/'
	# print 'outdir', outdir
	packval = 'FALSE'
	if LINPACK: packval = 'TRUE'
	
	if reuse: 
		rmtmp = False
		
		# check if files are present
		if not os.access(outdir+'Umat.txt', os.F_OK) or not os.access(outdir+'Smat.txt', os.F_OK) or not os.access(outdir+'Vtmat.txt', os.F_OK):
			reuse = False
		
	
	if not reuse:
		# hook up for using SVD from R
		# print 'Calling SVD in R'
		io.printTable(X,file=outdir+'emat.txt')
		scpt = """X=read.delim(\""""+outdir+'emat.txt'+"""\", header=FALSE)
		s=svd(X,LINPACK="""+packval+""")
		write.table(s$u,file=\""""+outdir+'Umat.txt'+"""\",sep=\"\t\")
		write.table(s$d,file=\""""+outdir+'Smat.txt'+"""\",sep=\"\t\")
		write.table(t(s$v),file=\""""+outdir+'Vtmat.txt'+"""\",sep=\"\t\")
		"""
		fh = open(outdir+'Rscpt.R','w')
		print >> fh, scpt
		fh.close()
		os.system('R --vanilla < "'+outdir+'Rscpt.R'+'" > "'+outdir+'RSVD.out.R'+'"')
		os.system('rm '+outdir+'Rscpt.R'); os.system('rm '+outdir+'RSVD.out.R')
		os.system('rm -f "%s"'%(outdir+'emat.txt'))
		
	try:
		U,crap,crap = io.readMatrix(outdir+'Umat.txt',dtype=float)
		S,crap,crap = io.readMatrix(outdir+'Smat.txt',dtype=float)
		Vt,crap,crap = io.readMatrix(outdir+'Vtmat.txt',dtype=float)
	
		if rmtmp: 
			os.system('rm -f '+outdir+'Umat.txt')
			os.system('rm -f '+outdir+'Smat.txt')
			os.system('rm -f '+outdir+'Vtmat.txt')
	
		U = np.array(U); S = np.array(S); Vt = np.array(Vt)
		if flattenS: S = np.array(flatten(S))
	    # print 'S test R', S
		if returnas == 'numpy': return U, S, Vt
		else: return U.tolist(), S.tolist(), Vt.tolist()
	except IOError:
		print 'cannot access SVD output files in directory %s'%(outdir)
		return -1


def covmat(m,y=None,rowvar=False,bias=False,ret='mat'):
	"""
	Estimate the covariance matrix.
	
	If m is a vector, return the variance. For matrices where each row is an 
	observation, and each column a variable, return the covariance matrix. 
	Note that in this case diag(cov(m)) is a vector of variances for each column.
	
	cov(m) is the same as cov(m, m)
	
	Normalization is by (N-1) where N is the number of observations 
	(unbiased estimate). If bias is True then normalization is by N.
	
	If rowvar is False, then each row is a variable with observations 
	in the columns.
	"""
	# from scipy.stats import stats
	import numpy as np
	mnp = np.array(m)
	if ret == 'list': return np.cov(mnp,y,rowvar,bias).tolist()
	elif ret == 'mat': return np.cov(mnp,y,rowvar,bias)
	return np.matrix(np.cov(mnp,y,rowvar,bias))

def cormat(m, y=None, rowvar=False, bias=False, tolist=False):
	"""Estimate a correlation matrix"""
	import numpy as np
	covxy = covmat(m,y,rowvar=rowvar,bias=bias)
	# divide each by variance products
	sdvector = np.sqrt(np.diagonal(covxy))
	r = covxy/np.outer(sdvector, sdvector)
	if tolist: return r.tolist()
	return r


def weightedCorrelation(X, Y, byRow=True):
	vals = []
	import numpy as np
	tmpX = np.array(X)
	tmpY = np.array(Y)
	for j in range(tmpX.shape[0]): vals += [correlation(tmpX[j,:], tmpY[j,:])]
	corrcoef = avg(vals)
	# print 'correlation between original and shuffled', corrcoef
	return corrcoef
	


def svdByTimeOrig(Xall, nS, T, basis, numAxes=-1, dirs=[], verbose=0, sortvectors=True, fixdirection=False, invert=False):
	"""Breaks the Xall matrix into T groups, where each group has nS columns.\n
	If a basis is provided, project each group onto this basis, determine the
	top numAxes directions of the basis shows greatest variance for that group,
	and return the data coordinates of that group on those axes, and the 
	variance-ranked directions.
	dirs is a list of direction-lists, where the length = T, and each direction-list
	indicates a ranked list of the particular axes of basis that each group's data
	should be projected onto.
	"""
	import numpy as np
	import scipy.linalg as la
	
	goodgenes = Xall.shape[0]
	numSamples = len(T) # number of groups
	numColumns = nS  # number of columns per group
	if invert:
		numSamples = nS
		numColumns = len(T)
	
	Ys = []; 
	if verbose: print 'data has', numSamples, goodgenes, 'x', numColumns, 'matrices'
	
	XbyT = partitionMAbyTime(Xall,nS,T, verbose=0)
	print 'inside, ', len(XbyT), 'numsamples', numSamples
	if basis.tolist():
		basisT = np.transpose(basis)[0:Xall.shape[1],:]
		origBasisT = np.array(basisT) # backup
		if verbose: print 'basisT.shape', basisT.shape
		tmpYs = []; sigma = []
		
		TOTAL = np.zeros([basisT.shape[0],basisT.shape[0]])
		
		# if dirs is given, want to project on each one
		# print 'comparison', len(dirs), len(T)
		addnewdirs = True
		if len(dirs)>0: 
			numSamples = len(dirs)
			sortvectors=False # just use the directions given in dirs
			addnewdirs = False
		print 'numsamples', numSamples
		for t in range(numSamples):
			
			# recover matrix of strain-replicates of this time point / group
			# --------------------------------------------------------------
			X = XbyT[t]
			
			# offset = numColumns*t
			# if not addnewdirs:
				# then we don't actually have more data, keep using same data set but
				# increment the directions
				# offset = 0
			# X = Xall[:,offset:offset+numColumns]
			
			
			if verbose: 
				print '\nt:', t, 'X.shape', X.shape
				print 'X is', X
			
			# print 'basisT is', basisT
			# empirical mean of matrix
			b = []
			for var in range(goodgenes):
				b.append(np.mean(X[var,:]))
				X[var,:] -= b[var]
			
			# project dataset onto global basis axes
			# --------------------------------------
			# print '      basis.T.shape',basisT.shape, 'orig', origBasisT.shape
			if len(dirs) == numSamples:
				print 'yes len(dirs)==numSamples'
				# project onto particular directions
				# newBasisT = np.zeros([len(dirs[t]),basisT.shape[1]])
				# locdir = []
				# for bi in range(len(dirs[t])):
				# 	locdir.append(dirs[t][bi])
				# 	newBasisT[bi,:] = origBasisT[dirs[t][bi],:]
				# basisT = np.array(newBasisT)
			# else: print 'ODD, dirs provided but doesnt match len(T)!'
			
			
			
			if verbose: print '      basis.T.shape',basisT.shape, 'X.shape', X.shape
			Y = np.dot(basisT, X) # 180 x 6300 . 6300 x 10 = 180 x 10
			if verbose: 
				print 'TESTING'
				print ' - Y', Y
				print ' - basisT', basisT
			
			# now want to determine the vector direction of most variance
			# for this time point in this common basis space
			U, S, Vt = SVD(Y)
			# U, S, Vt = la.svd(Y)
			# print 'SVD(Y):'
			print 'SVD', 'U',U.shape, 'x S',S.shape, 'x Vt', Vt.shape
			# get the axis orientations
			
			tmpYs.append(Y)
			if verbose: print '      Y.shape', Y.shape
			# if verbose: print 'Y is ', Y
			
			# covariance matrix of this reduced space data
			# --------------------------------------------
			sigma.append(covmat(Y,rowvar=True))
			if verbose: print '      sigma.shape', sigma[t].shape
			
			# select top numAxes eigenvariables for *this* time point
			# sorted by diagonal variance
			# -------------------------------------------------------
			ST = [sigma[t][i,i] for i in range(sigma[t].shape[1])]
			pairs = [[ST[i],i] for i in range(len(ST))]
			if sortvectors: pairs.sort(); pairs.reverse()
			if verbose: print 'pairs', pairs
			# recover ordering
			perm = [pairs[i][1] for i in range(len(pairs))]
			perm = perm[0:numAxes] # take top numAxes
			if fixdirection: perm = [0,1,2] # constant
			if addnewdirs: dirs.append(perm)
			# print 'Index', t, 'top eigen-variables', perm
			
			Ysub = np.empty([numAxes,Y.shape[1]])
			if verbose: print 'Ysub.shape', Ysub.shape
			# if verbose: print 'Y now?', Y
			for i in range(numAxes):
				if verbose: print 'buildilng', i, perm[i], Y[perm[i],:], Y[i,:]
				Ysub[i,:] = Y[perm[i],:]
			Ys.append(Ysub)
			# det = la.det(Ysub)
			# print 't',t, 'determinant', det
			# print t,'eigenvalues', w
			# if verbose: print 'Ysub'
			# print Ysub
		# restore basis
		# basisT = np.array(origBasisT) # restore original basis
		print
			
		
		
		# # now determine variable permutation to select
		# # top numAxes eigenvariables common over all time points
		# # sorting by total variance summed over time points
		###########################################################
		# ST = [TOTAL[i,i] for i in range(TOTAL.shape[1])]
		# pairs = [[ST[i],i] for i in range(len(ST))]
		# if sortvectors:
		# 	pairs.sort(); pairs.reverse()
		# print 'pairs', pairs
		# # recover ordering
		# perm = [pairs[i][1] for i in range(len(pairs))]
		# perm = perm[0:numAxes] # take top numAxes
		# print 'top eigen-variables', perm
		
		# # select subset of each Y corresponding to vars of top numAxes
		# # variables
		# for t in range(numSamples):
		# 	Ysub = np.empty([numAxes,tmpYs[t].shape[1]])
		# 	print 'Ysub.shape', Ysub.shape
		# 	for i in range(numAxes):
		# 		Ysub[i,:] = tmpYs[t][perm[i],:]
		# 	Ys.append(Ysub)
		# 	# det = la.det(Ysub)
		# 	# print 't',t, 'determinant', det
		# 	# print t,'eigenvalues', w
		# 	print 'Ysub'
		# 	print Ysub
		############################################################
	else:
		# run svd on each matrix individually
		for t in range(numSamples):
			# offset = numColumns*t
			# recover matrix of strain-replicates of this time point / group
			# X = Xall[:,offset:offset+numColumns]
			X = XbyT[t]
			print '\nt:', t, 'X.shape', X.shape
			
			# empirical mean of matrix
			b = []
			for var in range(goodgenes):
				b.append(np.mean(X[var,:]))
				X[var,:] -= b[var]
		
			# PCA transform
			U,S,Vt = SVD(X)
			# U,S,Vt = la.svd(X)
			if verbose: print 'SVD', 'U:', len(U), len(U[0]), 'x S:', len(S), 'x Vt:', len(Vt), len(Vt[0])
			
			# energy distribution
			singsum = 0
			# fractions = S[1:len(S)]**2/sum(S[1:len(S)]**2)
			fractions = (S**2/sum(S**2)).tolist()
			# print 'fractions', fractions
			entropy = -sum([fractions[i]*lg(fractions[i]) for i in range(len(fractions))])
			singsum = singtot = sum(S)
			singprop = singsum/float(singtot)
			print "Entropy (bits): Time %d Entropy %f" % (T[t], entropy)
			print "Sum %f Total %f Prop. %f" % (singsum, singtot, singprop)
	
			# original data projected to |S|-dim space
			Ut = np.transpose(U)[0:X.shape[1],:]
			Y = np.dot(Ut, X)
			print 'Y.shape', Y.shape
			
			# vcm = covmat(Y,rowvar=True)
			# print 'covmat'
			# print vcm
			# sys.exit()
			Ys.append(Y)
	return Ys, dirs


def svdByTimeSingle(Xall, samples, basis, numAxes=-1, dirs=[], verbose=0, sortvectors=True, fixdirection=False):
	"""Breaks the Xall matrix into T groups, where each group has nS columns.
	If a basis is provided, project each group onto this basis, determine the
	top numAxes directions of the basis shows greatest variance for that group,
	and return the data coordinates of that group on those axes, and the 
	variance-ranked directions.
	dirs is a list of direction-lists, where the length = T, and each direction-list
	indicates a ranked list of the particular axes of basis that each group's data
	should be projected onto.
	"""
	import numpy as np
	import scipy.linalg as la
	
	goodgenes = Xall.shape[0]
	numSamples = 1 # number of groups
	numColumns = len(samples)  # number of columns per group
	
	Ys = []; 
	if verbose: print 'data has', numSamples, goodgenes, 'x', numColumns, 'matrices'
	
	XbyT = [Xall]
	print 'inside, ', len(XbyT), 'numsamples', numSamples
	if basis.tolist():
		basisT = np.transpose(basis)[0:Xall.shape[1],:]
		origBasisT = np.array(basisT) # backup
		if verbose: print 'basisT.shape', basisT.shape
		tmpYs = []; sigma = []
		
		TOTAL = np.zeros([basisT.shape[0],basisT.shape[0]])
		
		# if dirs is given, want to project on each one
		# print 'comparison', len(dirs), len(T)
		addnewdirs = True
		if len(dirs)>0: 
			numSamples = len(dirs)
			sortvectors=False # just use the directions given in dirs
			addnewdirs = False
		print 'numsamples', numSamples
		for t in range(numSamples):
			
			# recover matrix of strain-replicates of this time point / group
			# --------------------------------------------------------------
			X = XbyT[t]
			
			# offset = numColumns*t
			# if not addnewdirs:
				# then we don't actually have more data, keep using same data set but
				# increment the directions
				# offset = 0
			# X = Xall[:,offset:offset+numColumns]
			
			
			if verbose: 
				print '\nt:', t, 'X.shape', X.shape
				print 'X is', X
			
			# print 'basisT is', basisT
			# empirical mean of matrix
			b = []
			for var in range(goodgenes):
				b.append(np.mean(X[var,:]))
				X[var,:] -= b[var]
			
			# project dataset onto global basis axes
			# --------------------------------------
			# print '      basis.T.shape',basisT.shape, 'orig', origBasisT.shape
			if len(dirs) == numSamples:
				print 'yes len(dirs)==numSamples'
				# project onto particular directions
				# newBasisT = np.zeros([len(dirs[t]),basisT.shape[1]])
				# locdir = []
				# for bi in range(len(dirs[t])):
				# 	locdir.append(dirs[t][bi])
				# 	newBasisT[bi,:] = origBasisT[dirs[t][bi],:]
				# basisT = np.array(newBasisT)
			# else: print 'ODD, dirs provided but doesnt match len(T)!'
			
			
			
			if verbose: print '      basis.T.shape',basisT.shape, 'X.shape', X.shape
			Y = np.dot(basisT, X) # 180 x 6300 . 6300 x 10 = 180 x 10
			if verbose: 
				print 'TESTING'
				print ' - Y', Y
				print ' - basisT', basisT
			
			# now want to determine the vector direction of most variance
			# for this time point in this common basis space
			U, S, Vt = SVD(Y)
			# U, S, Vt = la.svd(Y)
			# print 'SVD(Y):'
			print 'SVD', 'U',U.shape, 'x S',S.shape, 'x Vt', Vt.shape
			# get the axis orientations
			
			tmpYs.append(Y)
			if verbose: print '      Y.shape', Y.shape
			# if verbose: print 'Y is ', Y
			
			# covariance matrix of this reduced space data
			# --------------------------------------------
			sigma.append(covmat(Y,rowvar=True))
			if verbose: print '      sigma.shape', sigma[t].shape
			
			# select top numAxes eigenvariables for *this* time point
			# sorted by diagonal variance
			# -------------------------------------------------------
			ST = [sigma[t][i,i] for i in range(sigma[t].shape[1])]
			pairs = [[ST[i],i] for i in range(len(ST))]
			if sortvectors: pairs.sort(); pairs.reverse()
			if verbose: print 'pairs', pairs
			# recover ordering
			perm = [pairs[i][1] for i in range(len(pairs))]
			perm = perm[0:numAxes] # take top numAxes
			if fixdirection: perm = [0,1,2] # constant
			if addnewdirs: dirs.append(perm)
			# print 'Index', t, 'top eigen-variables', perm
			
			Ysub = np.empty([numAxes,Y.shape[1]])
			if verbose: print 'Ysub.shape', Ysub.shape
			# if verbose: print 'Y now?', Y
			for i in range(numAxes):
				if verbose: print 'buildilng', i, perm[i], Y[perm[i],:], Y[i,:]
				Ysub[i,:] = Y[perm[i],:]
			Ys.append(Ysub)
			# det = la.det(Ysub)
			# print 't',t, 'determinant', det
			# print t,'eigenvalues', w
			# if verbose: print 'Ysub'
			# print Ysub
		# restore basis
		# basisT = np.array(origBasisT) # restore original basis
		print
			
		
		
		# # now determine variable permutation to select
		# # top numAxes eigenvariables common over all time points
		# # sorting by total variance summed over time points
		###########################################################
		# ST = [TOTAL[i,i] for i in range(TOTAL.shape[1])]
		# pairs = [[ST[i],i] for i in range(len(ST))]
		# if sortvectors:
		# 	pairs.sort(); pairs.reverse()
		# print 'pairs', pairs
		# # recover ordering
		# perm = [pairs[i][1] for i in range(len(pairs))]
		# perm = perm[0:numAxes] # take top numAxes
		# print 'top eigen-variables', perm
		
		# # select subset of each Y corresponding to vars of top numAxes
		# # variables
		# for t in range(numSamples):
		# 	Ysub = np.empty([numAxes,tmpYs[t].shape[1]])
		# 	print 'Ysub.shape', Ysub.shape
		# 	for i in range(numAxes):
		# 		Ysub[i,:] = tmpYs[t][perm[i],:]
		# 	Ys.append(Ysub)
		# 	# det = la.det(Ysub)
		# 	# print 't',t, 'determinant', det
		# 	# print t,'eigenvalues', w
		# 	print 'Ysub'
		# 	print Ysub
		############################################################
	else:
		# run svd on each matrix individually
		for t in range(numSamples):
			# offset = numColumns*t
			# recover matrix of strain-replicates of this time point / group
			# X = Xall[:,offset:offset+numColumns]
			X = XbyT[t]
			print '\nt:', t, 'X.shape', X.shape
			
			# empirical mean of matrix
			b = []
			for var in range(goodgenes):
				b.append(np.mean(X[var,:]))
				X[var,:] -= b[var]
		
			# PCA transform
			U,S,Vt = SVD(X)
			# U,S,Vt = la.svd(X)
			if verbose: print 'SVD', 'U:', len(U), len(U[0]), 'x S:', len(S), 'x Vt:', len(Vt), len(Vt[0])
			
			# energy distribution
			singsum = 0
			# fractions = S[1:len(S)]**2/sum(S[1:len(S)]**2)
			fractions = (S**2/sum(S**2)).tolist()
			# print 'fractions', fractions
			entropy = -sum([fractions[i]*lg(fractions[i]) for i in range(len(fractions))])
			singsum = singtot = sum(S)
			singprop = singsum/float(singtot)
			print "Entropy (bits): Time %d Entropy %f" % (T[t], entropy)
			print "Sum %f Total %f Prop. %f" % (singsum, singtot, singprop)
	
			# original data projected to |S|-dim space
			Ut = np.transpose(U)[0:X.shape[1],:]
			Y = np.dot(Ut, X)
			print 'Y.shape', Y.shape
			
			# vcm = covmat(Y,rowvar=True)
			# print 'covmat'
			# print vcm
			# sys.exit()
			Ys.append(Y)
	return Ys, dirs


def svdByTime(Xall, nS, T, basis, numAxes=-1, dirs=[], verbose=0, sortvectors=True, fixdirection=False, returnVt=False, center=True):
	"""Breaks the Xall matrix into T groups, where each group has nS columns.
	If a basis is provided, project each group onto this basis, determine the
	top numAxes directions of the basis shows greatest variance for that group,
	and return the data coordinates of that group on those axes, and the 
	variance-ranked directions.
	dirs is a list of direction-lists, where the length = T, and each direction-list
	indicates a ranked list of the particular axes of basis that each group's data
	should be projected onto.
	"""
	import numpy as np
	import scipy.linalg as la
	
	goodgenes = Xall.shape[0]
	numSamples = len(T) # number of groups
	numColumns = nS  # number of columns per group
	
	Ys = []; Us = []; topS = []; sUs = []
	if verbose: print 'data has', numSamples, goodgenes, 'x', numColumns, 'matrices'
	
	basisT = np.transpose(basis)[0:Xall.shape[1],:]
	origBasisT = np.array(basisT) # backup
	if verbose: print 'basisT.shape', basisT.shape
	
	XbyT = partitionMAbyTime(Xall,nS,T, verbose=0)
	
	# if dirs is given, want to project on each one
	if verbose: print 'comparison', len(dirs), len(T)
	addnewdirs = True
	if len(dirs)>0: 
		numSamples = len(dirs)
		sortvectors=False # just use the directions given in dirs
		addnewdirs = False
	for t in range(numSamples):
		
		# recover matrix of strain-replicates of this time point / group
		# --------------------------------------------------------------
		
		# offset = numColumns*t
		# if not addnewdirs:
			# then we don't actually have more data, keep using same data set but
			# increment the directions
			# offset = 0
		# X = Xall[:,offset:offset+numColumns]
		
		X = XbyT[t]
		if verbose: 
			print '\nt:', t, 'X.shape', X.shape
			print 'X is', X
		
		# center this matrix
		if center:
			for var in range(goodgenes): X[var,:] -= np.mean(X[var,:])
		
		# project dataset onto global basis axes
		# --------------------------------------
		# print '      basis.T.shape',basisT.shape, 'orig', origBasisT.shape
		if len(dirs) == numSamples:
			# project onto particular directions
			newBasisT = np.zeros([len(dirs[t]),basisT.shape[1]])
			if verbose: print 'newBasisT.shape', newBasisT.shape
			locdir = []
			for bi in range(len(dirs[t])):
				locdir.append(dirs[t][bi])
				# print 'dirs', dirs[t][bi], 'index', bi
				newBasisT[bi,:] = origBasisT[dirs[t][bi],:]
			basisT = np.array(newBasisT)
			# print 'ok, new directions', locdir
		# else: print 'ODD, dirs provided but doesnt match len(T)!'
		
		if verbose: print '      basis.T.shape',basisT.shape
		
		# project this time point's data into global basis
		Y = np.dot(basisT, X) # 180 x 6300 . 6300 x 10 = 180 x 10
		
		# need to recenter the Y?
		# empirical mean of matrix
		# b = []
		# for var in range(Y.shape[0]):
		# 	b.append(np.mean(Y[var,:]))
		# 	Y[var,:] -= b[var]
		
		if verbose: 
			print 'TESTING'
			print 'Y', Y.shape
			print 'basisT', basisT.shape
		
		# now want to determine the vector direction of most variance
		# for this time point in this local basis space
		U, S, Vt = SVD(Y)
		# U, S, Vt = la.svd(Y)
		if verbose: print 'SVD(Y):', 'U',U.shape, 'x S',S.shape, 'x Vt', Vt.shape
		# get the axis orientations U and data projected onto these axes
		newY = np.dot(np.transpose(U), Y)
		if verbose: print 'newY.shape', newY.shape
		Ys.append(Y)
		
		if returnVt: Us.append(Vt)
		else: Us.append(U[:,0:Y.shape[1]]) # don't need the other 162-9 axes
		
		# sUs.append(np.dot(U,S))
		
		# get variance proportion on top direction?
		# print 'S', S.tolist()
		topS.append(S.tolist())
		
		if verbose: print '      Y.shape', Y.shape
		# if verbose: print 'Y is ', Y
		
	print
	
	return Ys, Us, topS


def covmatByTimeIndiv(Xall, numStrains, T, numAxes=-1, rowvar=False, verbose=0):
	""""""
	import numpy as np
	import scipy.linalg as la
	goodgenes = Xall.shape[0]
	superT = T
	
	numSamples = len(superT) # number of groups
	numColumns = numStrains  # number of columns per group
	if verbose: print 'data has ', numSamples, goodgenes, 'x', numColumns, 'matrices'
	
	# XbyT = partitionMAbyTime(Xall,numStrains,T, verbose=1)
	
	SIGMA = []
	for t in range(numSamples):
		# offset = numColumns*t
		# recover matrix of strain-replicates of this time point / group
		# X = Xall[:,offset:offset+numColumns]
		# print '\nt:', t, 'X.shape', X.shape
		X = Xall[t]
		# X = XbyT[t]
		print '\nt:', t, 'X.shape', X.shape
		
		
		# empirical mean of matrix
		b = []
		for var in range(goodgenes):
			b.append(np.mean(X[var,:]))
			X[var,:] -= b[var]
		
		if numAxes > -1:
			# PCA transform
			U,S,Vt = SVD(X)
			# U,S,Vt = la.svd(X)
			if verbose: print 'SVD', 'U:', len(U), len(U[0]), 'x S:', len(S), 'x Vt:', len(Vt), len(Vt[0])
		
			# energy distribution
			singsum = 0
			# fractions = S[1:len(S)]**2/sum(S[1:len(S)]**2)
			fractions = (S**2/sum(S**2)).tolist()
			# print 'fractions', fractions
			entropy = -sum([fractions[i]*lg(fractions[i]) for i in range(len(fractions))])
			singsum = singtot = sum(S)
			if numAxes > 0: singsum = sum(S[0:numAxes])
			singprop = singsum/float(singtot)
			print "Entropy (bits) for %d axes: Time %d Entropy %f" % \
				(numAxes, superT[t], entropy)
			print "Sum %f Total %f Prop. %f" % (singsum, singtot, singprop)
		
			# original data projected to |S|-dim space
			Y = S*Vt 
			# want covariance matrix of this reduced space data
			vcm = []
			if numAxes == 0: vcm = covmat(Y, rowvar=rowvar)
			else: 
				whyy = Y[:,0:numAxes] # want all genes at 7 of 9 strains?
				if rowvar: 
					# if variance across row, axes are rows not cols
					whyy = Y[0:numAxes,:] # all observations of variables
				print 'shape for covmat', whyy.shape
				vcm = covmat(whyy, rowvar=rowvar)
			SIGMA.append(vcm) # c x c matrix, each col holds projected coords. of a strain
		else:
			SIGMA.append(covmat(X,rowvar=rowvar)) # covmat for all genes
		
	return SIGMA

def covmatByTime(Xall, nS, T, basis, numAxes=-1, sort=False, varlist=[], verbose=0):
	""""""
	import numpy as np
	import scipy.linalg as la
	goodgenes = np.array(Xall[0]).shape[0]
	# goodgenes = Xall.shape[0]
	superT = T
	numSamples = len(T) # number of groups
	numColumns = nS  # number of columns per group
	
	# basisT = np.transpose(basis)[0:numAxes,:]
	# basisT = np.transpose(basis)[0:Xall.shape[1],:]
	basisT = np.transpose(basis)
	
	# XbyT = partitionMAbyTime(Xall,nS,T, verbose=1)
	
	if verbose: 
		print 'data has', numSamples, goodgenes, 'x', numColumns, 'matrices'
		print 'basisT.shape', basisT.shape
	
	SIGMA = []; TOTAL = np.zeros([basisT.shape[0],basisT.shape[0]])
	print 'TOTAL.shape',TOTAL.shape
	
	for t in range(numSamples):
		
		# recover matrix of strain-replicates of this time point / group
		# offset = numColumns*t
		# X = Xall[:,offset:offset+numColumns]
		
		# X = XbyT[t]
		X = Xall[t]
		print '\nt:', t, 'X.shape', X.shape
		
		# print 'my X'
		# print X
		
		# U,S,Vt = la.svd(X)
		# print 't',t, S
		
		# center this time point
		# for var in range(goodgenes): X[var,:] -= np.mean(X[var,:])
		
		# dataset projected onto global basis axes
		print '      ', basisT.shape, '.', X.shape
		Y = np.dot(basisT, X) # 180 x 6300 . 6300 x 10 = 180 x 10
		
		# print 'projected'
		# print Y
		
		# covariance matrix of this reduced space data
		vcm = covmat(Y,rowvar=True)
		print '      vcm.shape', vcm.shape
		SIGMA.append(vcm)
		
		TOTAL += vcm # sum of covariances over times
		
		# OLD
		# if numAxes > -1:
		# 	# PCA transform
		# 	U,S,Vt = la.svd(X)
		# 	if verbose: print 'SVD', 'U:', len(U), len(U[0]), 'x S:', len(S), 'x Vt:', len(Vt), len(Vt[0])
		# 
		# 	# energy distribution
		# 	singsum = 0
		# 	# fractions = S[1:len(S)]**2/sum(S[1:len(S)]**2)
		# 	fractions = (S**2/sum(S**2)).tolist()
		# 	# print 'fractions', fractions
		# 	entropy = -sum([fractions[i]*lg(fractions[i]) for i in range(len(fractions))])
		# 	singsum = singtot = sum(S)
		# 	if numAxes > 0: singsum = sum(S[0:numAxes])
		# 	singprop = singsum/float(singtot)
		# 	print "Entropy (bits) for %d axes: Time %d Entropy %f" % \
		# 		(numAxes, superT[t], entropy)
		# 	print "Sum %f Total %f Prop. %f" % (singsum, singtot, singprop)
		# 
		# 	# original data projected to |S|-dim space
		# 	Y = S*Vt 
		# 	# want covariance matrix of this reduced space data
		# 	vcm = []
		# 	if numAxes == 0: vcm = covmat(Y, rowvar=rowvar)
		# 	else: 
		# 		whyy = Y[:,0:numAxes]
		# 		if rowvar: 
		# 			# if variance across row, axes are rows not cols
		# 			whyy = Y[0:numAxes,:] # all observations of variables
		# 		print 'shape for covmat', whyy.shape
		# 		vcm = covmat(whyy, rowvar=rowvar)
		# 	SIGMA.append(vcm) # c x c matrix, each col holds projected coords. of a strain
		# else:
		# 	SIGMA.append(covmat(X,rowvar=rowvar)) # covmat for all genes
		# OLD
		
	perm = range(numAxes)
	if sort:
		# determine variable permutation to select the
		# top numAxes eigenvariables common over all time points
		# sorting by total variance summed over time points
		ST = [TOTAL[i,i] for i in range(TOTAL.shape[1])]
		# ST = [sum(TOTAL[:,i].tolist()) for i in range(TOTAL.shape[1])]
		pairs = [[ST[i],i] for i in range(len(ST))]
		pairs.sort(); pairs.reverse()
		# recover ordering
		perm = [pairs[i][1] for i in range(len(pairs))]
		# take top numAxes
		perm = perm[0:numAxes]
	if varlist: perm = varlist
	
	print 'top eigen-variables', perm
	
	# select subset of each SIGMA corresponding to vars of top numAxes vars
	TOPSIGMA = []
	for t in range(numSamples):
		tmp = np.empty([numAxes,numAxes])
		for i in range(numAxes):
			for j in range(numAxes): tmp[i,j] = SIGMA[t][perm[i],perm[j]]
		TOPSIGMA.append(tmp)
		# print tmp
		det = la.det(tmp)
		print 't',t, 'determinant', det
		
		# check out eigenvalues
		# w,vr = la.eig(tmp)
		# print t,'eigenvalues', w
		
	return TOPSIGMA, perm


def varianceByDirection(X, directions, rank=1):
	"""Project a data set onto each of a list of directions 
	(data and vectors provided must have comparable dimension).
	"""
	import numpy as np
	var = [[nastr for j in range(rank)] for i in range(len(directions))]
	dof = [[nastr for j in range(rank)] for i in range(len(directions))]
	prop = [[nastr for j in range(rank)] for i in range(len(directions))]
	total = [0 for i in range(len(directions))]
	# print 'there are ', len(directions), 'directions'
	# print 'what is rank?', rank
	for i in range(len(directions)): # i.e. time points
		for j in range(rank): # eigenvectors - all 162 of them
		
			d = directions[i][:,j].tolist() # CSA axis at this time i, rank j
			proj = []
			if len(X) > 1:
				for k in range(X[i].shape[1]):
					# project strain k at time i onto direction d
					proj += [dot(X[i][:,k].tolist(), d)]
			else: 
				# print 'odd shape', X[0].shape[1], 'is this 23?'
				for k in range(X[0].shape[1]):
					# project strain k at time i onto direction d
					proj += [dot(X[0][:,k].tolist(), d)]
				# print 'does var have 23 items?', len(proj)
			
			var[i][j] = variance(proj)
			total[i] += var[i][j]
			dof[i][j] = len(proj)
		
		for j in range(rank):
			prop[i][j] = var[i][j]/float(total[i])
			
	return var, prop, dof


entropy = lambda X: -sum([X[i]*math.log(X[i],2) for i in range(len(X))])

# CLUSTERING ALGORITHMS
# ---------------------
euclidean = lambda x,y: math.sqrt(sum([(x[i]-y[i])**2 for i in range(len(x))]))
def qtcluster(file, P=0.001, radius=-1, outdir='', header=True, rownames=True):
	"""Perform a QT clustering on matrix of input data. If radius is not given, 
	compute the allpairs distance distribution and return value at percentile P.
	Z scales the data. Euclidean distance.
	"""
	import sio as io
	intna = lambda x: x == 'NA' and 'NA' or int(x)
	
	print 'Performing QT clustering...'
	
	rownamestr = 'x = x[,2:length(x)]'
	if rownames == False: rownamestr = ''
	
	header2 = 'TRUE'
	if not header: header2 = 'FALSE'
	scpt = """
	library(flexclust)
	require(stats)
	x = read.csv("%s", sep="\t", header=%s)
	%s
	# x = as.numeric(as.matrix(x))
	x = scale(x)
	radius = %s
	if (radius == -1) {
		D = dist(x)
		
		dim(D) <- NULL # make vector
		p = c()
		for (i in seq(25)) {
			p = append(p,quantile(sample(D,100000),%f))
		}
		radius = min(p)
		D = 0
	}
	radius
	graph = qtclust(x, radius=radius)
	clusterInfo = attr(graph,'cluster')
	write(clusterInfo, file="%s", sep="\t")
	
	""" % (file, header2, rownamestr, radius, P, outdir+'QTclusterdat.txt')
	
	fh = open(outdir+'rscript.R', 'w')
	print >> fh, scpt
	fh.close()
	os.system('R --vanilla < %s > %s' % (outdir+'rscript.R', outdir+'rscript.R.out'))
	d,r,c = io.readMatrix(outdir+'QTclusterdat.txt', header=False, rowNames=False)
	
	# -----------------------------
	
	clusterNums = []
	for row in d: clusterNums += row
	# create lists of genes names corresponding to clusters
	d,genes,c = io.readMatrix(file,header=header)
	
	clusterInfo = {}
	for i in range(len(clusterNums)):
		gene = genes[i]
		clusterNum = intna(clusterNums[i])
		# print 'have clusternum', clusterNum
		if clusterNum not in ['NA', nastr]:
			try: clusterInfo[clusterNum] += [gene]
			except KeyError: clusterInfo[clusterNum] = [gene]
		# else: print 'one cluster is NA, skipping...'
	
	# add trans size information to each cluster
	# remainder = clusterTransSize(clusterInfo)
	
	# save to disk
	# if clusterDir:
	# if not clusterFile:
	# 	parts = file.split('/')
	# 	if usesimoladata: 
	# 		crap,name,crap = re.split('.+_(YPS.+)_.+', parts[len(parts)-1])
	# 	else:
	# 		crap,name,crap = re.split('(.+)\.txt', parts[len(parts)-1])
	# 	clusterFile = outdir+'QTclusters_%s.txt' % (name)
	
	# if clusterFile:
	# 	fh = open(clusterFile, 'w')
	# 	keys = clusterInfo.keys(); keys.sort()
	# 	for k in keys:
	# 		print >> fh, k, clusterInfo[k]['genes']
	# 	fh.close()
	
	# os.system('rm -f %s' % (outdir+'rscript.R'))
	# os.system('rm -f %s' % (outdir+'rscript.R.out'))
	# os.system('rm -f %s' % (outdir+'QTclusterdat.txt'))
	
	return {'clusters':clusterInfo, 'radius':radius, 'P':P}

def hcluster(D, method='average'):
	print 'not sure this is good yet...'
	import hcluster as hc
	Z = None
	if method == 'linkage': Z = hc.linkage(D)
	elif method == 'single': Z = hc.single(D)
	elif method == 'complete': Z = hc.complete(D)
	elif method == 'weighted': Z = hc.weighted(D)
	elif method == 'centroid': Z = hc.centroid(D)
	elif method == 'median': Z = hc.median(D)
	elif method == 'ward': Z = hc.ward(D)
	hc.dendrogram(Z)
	

def kmeans(dat, labels=[], header=True, clusters=5, method='Euclidean', maxiters=20, nstart=1, outdir='', quick=False, rmtmp=True):
	import sio as io, os
	intna = lambda x: x == 'NA' and 'NA' or int(x)
	
	# print 'Performing kmeans clustering (k=%s)...'%(clusters)
	
	file = outdir+'kmeansdata.txt'
	if type(dat) == type([]):
		# need to print to file
		if not labels: labels = range(1,len(dat)+1)
		io.printMatrixFast(dat,labels,file=file)
		header=False
	else: file = dat
	
	HAVECENTERS = False
	if type(clusters) == type([]):
		# then we have been provided with the actual centers
		io.printMatrixFast(clusters,file=outdir+'centers.txt')
		HAVECENTERS = True
	
	
	header2 = 'TRUE'
	if not header: header2 = 'FALSE'
	
	EUCLID = 'TRUE'
	# other distance
	if method == "binary" or method == "Binary": EUCLID = 'FALSE'
	
	if not HAVECENTERS:
		scpt = """
		require(stats)
		out = "%s"
		x = read.csv("%s", sep="\t", header=%s, na.strings='"NA"')
		x = x[,2:length(x)]
		if (%s) { kmc = kmeans(x, centers=%s, iter.max=%s, nstart=%s); } else {
			library(amap)
			kmc = Kmeans(x, centers=%s, iter.max=%s, nstart=%s, method="spearman")
		}
		write(kmc$withinss, file=paste(out,"withinss.txt",sep=""), sep="\t")
		write(kmc$centers, file=paste(out,"centers.txt",sep=""), ncolumns=length(x), sep="\t")
		write(kmc$cluster, file=paste(out,"clusterInfo.txt",sep=""), sep="\t")
		""" % (outdir, file, header2, EUCLID, clusters, maxiters, nstart, clusters, maxiters, nstart)
		
	else:
		scpt = """
		require(stats)
		out = "%s"
		x = read.csv("%s", sep="\t", header=%s)
		x = x[,2:length(x)]
		
		# read centers
		c = read.delim("%s",header=FALSE)
		
		if (%s) {
			kmc = kmeans(x, centers=c, iter.max=%s, nstart=%s)
		}
		else {
			library(amap)
			kmc = Kmeans(x, centers=c, iter.max=%s, nstart=%s, method="binary")
		}
		
		write(kmc$withinss, file=paste(out,"withinss.txt",sep=""), sep="\t")
		write(kmc$cluster, file=paste(out,"clusterInfo.txt",sep=""), sep="\t")
		""" % (outdir, file, header2, EUCLID, outdir+'centers.txt', maxiters, nstart, outdir+'centers.txt', maxiters, nstart)
	#, '%*%'
	
	CREATEDDIR = False
	if not os.access(outdir, os.F_OK): 
		io.createdir(outdir)
		CREATEDDIR = True
	
	fh = open(outdir+'rscript.R', 'w')
	print >> fh, scpt
	fh.close()
	os.system('R --vanilla < %s > %s' % (outdir+'rscript.R', outdir+'rscript.R.out'))
	dat = {'clusterInfo':[], 'centers':[], 'withinss':[]}
	for k in dat.keys():
		dat[k] = io.readMatrix(outdir+k+'.txt', header=False, rowNames=False)
	
	error = sum(map(float, flatten(dat['withinss'][0])))
	if quick: 
		if CREATEDDIR and os.access(outdir, os.F_OK): os.system('rm -r "%s"'%(outdir))
		else:
			os.system('rm -f %s' % (outdir+'rscript.R'))
			os.system('rm -f %s' % (outdir+'rscript.R.out'))
			os.system('rm -f %s' % (outdir+'withinss.txt'))
			os.system('rm -f %s' % (outdir+'centers.txt'))
			os.system('rm -f %s' % (outdir+'clusterInfo.txt'))
			if file != dat: os.system('rm -f %s'%(outdir+'kmeansdata.txt'))
		return {'error':error, 'clusterInfo':{}, 'centers':[], 'indicator':[], 'lookup':{}}
	
	clusterNums = []
	for row in dat['clusterInfo'][0]: clusterNums += row
	clusterNums = map(intna,clusterNums)
	
	# print 'len', dat['centers'], 'len', dat['centers'][0]
	centers = map(lambda v: map(float, v), dat['centers'][0])
	# print 'centers size', len(centers), len(centers[0])
	
	# create lists of genes names corresponding to clusters
	if type(dat) != type([]):
		dat,labels,c = io.readMatrix(file,header=header)
	
	clusterInfo = {}
	lookup = {} # keyed by label returning cluster group
	for i in range(len(clusterNums)):
		gene = labels[i]
		# print 'gene', gene
		if clusterNums[i] not in ['NA', nastr]:
			try: clusterInfo[clusterNums[i]] += [gene]
			except KeyError: clusterInfo[clusterNums[i]] = [gene]
			
			lookup[gene] = clusterNums[i]
		# else: print 'one cluster is NA, skipping...'
	
	if CREATEDDIR and os.access(outdir, os.F_OK): os.system('rm -r "%s"'%(outdir))
	elif rmtmp:
		os.system('rm -f %s' % (outdir+'rscript.R'))
		os.system('rm -f %s' % (outdir+'rscript.R.out'))
		os.system('rm -f %s' % (outdir+'withinss.txt'))
		os.system('rm -f %s' % (outdir+'centers.txt'))
		os.system('rm -f %s' % (outdir+'clusterInfo.txt'))
		if file != dat: os.system('rm -f %s'%(outdir+'kmeansdata.txt'))
	
	return {'clusters':clusterInfo, 'error':error, 'clusterIDs':clusterNums, 'lookup':lookup, 'centers':centers}

def kmeansFull(dat, dim=-1, nvars=-1, labels=[], header=True, minclusters=2, maxiter=10, nstart=1, method='Euclidean', penalty='BIC', lamb=1, STABLEMIN=20, outdir='', quick=False, verbosity=1):
	"""kmeansFull: verbosity = 0, 1, or 2"""
	
	import sio as io,os
	
	olabels = labels
	if len(labels) and header == True:
		labels = ['Key'] + labels
	
	CREATEDIR = False
	if not outdir: 
		outdir = os.getcwd()+'/tmp/'
		CREATEDIR = True
	io.createdir(outdir)
	
	# if type(dat) == type(''):
	# 		# then we are given a file
	# 		print 'imputing data...'
	# 		d,labels,c = io.readMatrix(dat,rownames,header,dtype=floatna)
	# 		for i in range(len(d)):
	# 			for j in range(len(d[i])):
	# 				if d[i][j] == nastr: d[i][j] = 0
	# 	
	# 	
	# d = dat; r = labels
	
	# print 'test', len(d), len(d[0])
	
	# dim = len(d[0]) # dimensionality
	# if nvars == -1: nvars = len(d)
	# if verbosity > 0: print 'Data set: dim=%s, nvars=%s'%(dim,nvars)
	
	file = outdir+'kmeansFull.txt'
	if type(dat) == type([]):
		dim = len(dat[0]) # dimensionality
		nvars = len(dat)
		if verbosity > 0: print 'Data set: dim=%s, nvars=%s'%(dim,nvars)
			
		# this will be faster if i save to disk now
		io.printMatrixFast(dat,labels,file=file)
	else: 
		file = dat
		if dim == -1 or nvars == -1: 
			print 'Error: when file is provided dim and nvars must also be specified.'
			return -1
	
	d = None; r = None; c = None
	
	ks = []; errors = []; BICs = []
	
	# iterative, testing for BIC criterion
	MINerror = GLOB_INF
	ATMINIMUM = 0
	MINK = minclusters # starting
	MINcenters = []
	MINRSS = None
	MINgeneclustering = {}
	
	results = None
	E = math.e
	k = MINK
	
	firstRSS = None
	while ATMINIMUM < STABLEMIN: # when to decide we're at minimum
		loutdir = outdir+'kmeans-%s/'%(k)
		io.createdir(loutdir)
		
		n = float(nvars)
		# print 'kmeans on', file
		res = kmeans(file, labels=labels, header=header, clusters=k, maxiters=maxiter, 
					nstart=nstart, method=method, outdir=loutdir, quick=quick)
		
		RSS = res['error']
		
		base = n + n*math.log(2*math.pi, E) + n*math.log(RSS/n, E)
		if firstRSS == None: firstRSS = math.log(RSS/n, E)
		
		fit = base
		
		AIC = base + lamb*(2*(k*dim + 1))
		BIC = base + lamb*(math.log(n, E)*(k*dim + 1))
		np = k*dim
		AICc = base + lamb*(2*(np+1) + (2*np*(np+1))/(n-np-1))
		
		if penalty == 'AIC': fit = AIC
		elif penalty == 'BIC': fit = BIC
		elif penalty == 'AICc': fit = AICc
		
		if fit < 0: fit = MINerror
		
		# want to stop if start getting negative signs
		if fit >= 0.0 and fit < MINerror:
		# if fit < MINerror:
			results = res
			MINerror = fit
			MINRSS = RSS
			MINK = k
			MINcenters = res['centers']
			if not quick: MINgeneclustering = res['clusters']
			ATMINIMUM = 0
		elif fit >= 0.0: ATMINIMUM += 1
		elif firstRSS > 0: ATMINIMUM = STABLEMIN
		
		if verbosity > 0: print '- K=%s, RSS=%s, AIC=%s, AICc=%s, BIC=%s'%(k, RSS, AIC,AICc,BIC), 'atmin', ATMINIMUM
		
		# save the results
		ks += [k]
		errors += [RSS]
		BICs += [fit]
		
		# if not quick:
		# 	for c in res['clusters'].keys():
		# 		io.printList(res['clusters'][c], file=loutdir+'kmeans-%s.txt'%(c))
		
		# if quick: 
		os.system('rmdir "%s"'%(loutdir))
		k += 1
	
	if verbosity > 0: 
		print
		print '--------------------------------'
		print 'Best k found at', MINK, MINerror
		print '--------------------------------'
		print
		if verbosity > 1:
			print 'k =', ks
			print 'RSS =', errors
			print 'BIC =', BICs
			print 'DONE'
	
	os.system('rm -f '+outdir+'kmeansFull.txt')
	
	if CREATEDIR: os.system('rm -r "%s"'%(outdir))
	
	# provide a vector of clusterids in dat input order
	labeling = []
	if labels:
		# turn clusters into dict
		dtmp = {}
		for c in MINgeneclustering.keys():
			for g in MINgeneclustering[c]:
				dtmp[g] = c
		
		labeling = [dtmp[g] for g in olabels]
	
	if not quick:
		loutdir = outdir+'kmeans-%s/'%(MINK)
		io.createdir(loutdir)
		for c in MINgeneclustering.keys():
			io.printList(MINgeneclustering[c], file=loutdir+'kmeans-%s.txt'%(c))
	
	return {'k':MINK, 'error':MINerror, 'RSS':MINRSS, 'centers':MINcenters, 'clusters':MINgeneclustering, 'clusterLabels':labeling}



def resampleMatrix(M,ntimes=100,feature='row',withReplacement=True):
	# standard matrix resampling
	MM = M
	if feature == 'col': MM = transpose(M)
	
	R = [] # list of resampled matrices of same dimension
	for i in range(ntimes):
		Mnew = []
		while len(Mnew) < len(MM):
			Mnew += [random.choice(MM)]
		if feature == 'col': R += [transpose(Mnew)]
		else: R += [Mnew]
	return R

#
def renoiseMatrix(M,ntimes=100,sd=0.01,feature='row'):
	# randomly add or subtract noise to each data point, drawn from gaussian with sd = se
	MM = M
	if feature == 'col': MM = transpose(M)
	
	R = [] # list of resampled matrices of same dimension
	for i in range(ntimes):
		Mnew = []
		for j in range(len(MM)):
			row = map(lambda x: x+random.gauss(0,sd), MM[j])
			row = map(lambda x: x > 0 and x or 0, row)
			Mnew += [row]
		if feature == 'col': R += [transpose(Mnew)]
		else: R += [Mnew]
	return R

def subsampleMatrix(M,ntimes=100,nrows=100):
	# subsample a provided matrix ntimes
	return [random.sample(M,nrows) for i in range(ntimes)]



# def thresholdRescale(dat,niter=1, exp=3, fscale=None):
# 	"""Assumes variable values are symmetric."""
# 	if not fscale: fscale = max(ut.flatten(dat))
# 	if niter < 1: sys.exit("threholdRescale: niter must be a positive integer.")
# 	
# 	mat = dat
# 	# rounds of thresholding
# 	for i in range(niter):
# 		# cube everything to exaggerate values
# 		mat = [map(lambda x: x**exp, mat[g]) for g in range(len(mat))]
# 		# range scale
# 		mat = [map(lambda x: x*(fscale/fscale**exp), mat[g]) for g in range(len(mat))]
# 	return mat
# 

# need an R pvclust clustering method


# Missing data imputation
def imputeMissingData(mat, k=20, KNN=1, outdir='', rmtmp=True):
	"Imputes missing data using K-nearest neighbors approach. First filters rows with more than 50% missing values. Then uses a weighted average of KNN and cubic spline interpolation to impute remaining missing values. Default here is KNN only with k=20."
	
	# R code:
	# Jeremy Weiss, Kim Lab @ UPenn, 2 April 2007
	# modified by DFS for cubic.spl to take in a table with row/col headings
	### These functions require you to install the packages 'fields' and 'SeqKnn'
	#   click on packages --> install packages --> ...
	#
	#   Use the function 'impute' on the data matrix for the optimal found imputation method.
	
	# save matrix to file
	import sio
	fn = 'imputeMissingData.matrix.txt'
	sio.printMatrix(mat,file=outdir+fn)
	
	ncolumns = len(mat[0])
	newfn = 'imputeMissingData.matrix.imputed.txt'
	
	Rcode = """
library(fields)
library(SeqKnn)

removeHighNAgenes <- function(dataM) {
 nas <- rowSums(is.na(dataM))
 removenas <- (nas > (length(dataM[1,])/2))
 dataM[!removenas,]
}

cubic.spl <- function(dataWNA) {
 require(fields)
 missing <- is.na(apply(dataWNA,1,mean)) #find if row has a NA?
 wmiss <- which(missing)
 for(i in 1:sum(missing)) {
   missi <- is.na(dataWNA[wmiss[i],])
   a <- 0:(dim(dataWNA)[2]-1);       b <- dataWNA[wmiss[i],]
   a <- a[-which(missi)]; b <- b[-which(missi)]
   dataWNA[wmiss[i],] <- splint(a,b,0:(dim(dataWNA)[2]-1))
 }
 dataWNA
}

Seq.KNN <- function(matr, k) {
 require(SeqKnn)
 SeqKNN(matr,k)
}

impute <- function(matr) {
 nona <- removeHighNAgenes(matr)
 cs <- cubic.spl(nona)
 if(is.na(sum(cs[,2:dim(cs)[2]])))
   Seq.KNN(nona,%s)
 else
   %s*Seq.KNN(nona, %s) + (1-%s)*cs
}

# Now load the data and run code
setwd("%s")
dat=read.delim("%s", header=FALSE, na.string=%s)
write.table(impute(dat), file="%s", sep="\t", quote=FALSE, na=%s, col.names=NA)

	"""%(k, KNN, k, KNN, outdir, fn, nastr, newfn, nastr)
	
	# save and run script
	sio.printList([Rcode], file=outdir+'Rscpt.imputeMissingData.R')
	os.system('R --vanilla < "%s" > "%s"'%(outdir+'Rscpt.imputeMissingData.R',outdir+'Rscpt.imputeMissingData.out.R'))
	
	# load the matrix
	d,rows,crap = sio.readMatrix(outdir+newfn)
	
	if rmtmp:
		os.system('rm -f "%s"'%(outdir+fn))
		os.system('rm -f "%s"'%(outdir+newfn))
	
	
	return d,rows


# MATRIX FUNCTIONS
# ----------------
def zeros(n): 
	return [0 for i in range(n)]

def centerMatrix(m, func=lambda x: avg(x)):
	""""""
	return [map(subna1(func(m[i])), m[i]) for i in range(len(m))]

def centroid(X):
	"""
	Return the centroid vector of a matrix X[m x n]
	Assumes X[i] is a vector of the ith components
	of each RV X[:,j]. Thus the centroid is a m dim vector.
	"""
	
	b = []
	for i in range(len(m)):
		mu = avg(m[i])
		b.append(mu)
	return b

def npCentroid(X):
	"""Requires numpy array"""
	import numpy as np
	# MEAN CENTER EACH DIMENSION (genes across conditions)
	b = []
	# print 'X.shape', X.shape
	for var in range(X.shape[0]):
		b.append(np.mean(X[var,:]))
	return b


def partitionMA(X,partition, verbose=False):
	# Partition X to recover individual data sets
	try: import numpy as np
	except ImportError:
		try: import Numeric as np
		except ImportError:
			sys.exit('myNeighbors.py: Sorry, neither numpy nor Numeric are installed.')
	
	# if type(partition) == type(5):
	# 	partition = [partition for ]
	
	X = np.array(X)
	off = [0]
	if verbose: 
		print 'Partitioning data matrix:'
		print ' - partition ', 1, off[0], 'to' , off[0]+partition[0]
	
	for i in range(1,len(partition)): 
		off += [ off[i-1]+partition[i-1] ]
		if verbose: print ' - partition ', i+1, off[i], 'to' , off[i]+partition[i]
	
	mX = [X[:,off[i]:off[i]+partition[i]].tolist() for i in range(len(partition))]
	if verbose: 
		print 'Data:',
		for i in range(len(mX)): print i, '(%sx%s)'%(len(mX[i]), len(mX[i][0])),
		print '\n'
	return mX

def partitionMAbyTime(X, n, conditions, verbose=False):
	try: import numpy as np
	except ImportError:
		try: import Numeric as np
		except ImportError:
			sys.exit('myNeighbors.py: Sorry, neither numpy nor Numeric are installed.')
	
	X = np.array(X)
	Xall = []
	if verbose: print 'Partitioning data matrix by time:', 'n=',n, 'conditions=', len(conditions)
	nc = len(conditions)
	
	# group data by condition
	
	# SLOW VERSION
	# for c in range(nc):
	# 	tmp = []
	# 	if verbose: print ' - time index', c
	# 	for i in range(n):
	# 		if verbose: print '   - fetching', i, nc*i+c, 'to', nc*i+c+1
	# 		tmp += [flatten(X[:,n*i+c:n*i+c+1].tolist())]
	# 	tmpt = np.transpose(tmp)
	# 	Xall += [tmpt]
	
	# FAST VERSION USING NUMPY
	for c in range(nc):
		tmp = np.empty([X.shape[0],n])
		if verbose: print ' - time index', c
		for i in range(n):
			if verbose: print '   - fetching', i, nc*i+c, 'to', nc*i+c+1
			tmp[:,i:i+1] = X[:,n*i+c:n*i+c+1]
		Xall += [tmp]
	
	return Xall

def subsetMAData(X, labels, labelSubset):
	Y = []
	G = []
	for lab in labelSubset:
		try: 
			Y += [X[labels.index(lab)]]
			G += [lab]
		except ValueError: pass
	return Y, G


def shuffleMatrix(mat, byRow=True, byCol=True, useNumpy=False):
	"""Takes a data matrix and permutes the columns of each row
	(byRows=True) or the rows of each column (byRows=False)."""
	
	# numpy loader
	if useNumpy:
		try: 
			import numpy as np
			from numpy.random import shuffle
		except ImportError:
			try: 
				import Numeric as np
				from Numeric.random import shuffle
			except ImportError:
				sys.exit('myNeighbors.py: Sorry, neither numpy nor Numeric are installed.')
	
		rmat = np.array(mat)
		if byRow: 
			for i in range(rmat.shape[0]): shuffle(rmat[i,:])
		if byCol: 
			for i in range(rmat.shape[1]): shuffle(rmat[:,i])
		return rmat
	else:
		# use my own shuffler
		tmp = []
		if byRow: tmp = [permute(mat[i]) for i in range(len(mat))]
		if byCol:
			tmp = transpose(tmp)
			tmp = [permute(tmp[i]) for i in range(len(tmp))]
			tmp = transpose(tmp)
		return tmp

def uniqueEntries(mat,useDiagonal=False):
	vals = []
	for i in range(len(mat)):
		for j in range(i+int(useDiagonal)): vals += [mat[i][j]]
	return vals
uniqueMatrixEntries = uniqueEntries

def replaceDiagonal(mat,na=nastr):
	"""Replaces diagonal values in 2D matrix with na"""

	mat2 = [mat[i][:] for i in range(len(mat))]
	for i in range(len(mat2)):
		mat2[i][i] = na
	return mat2


def replaceOffDiagonal(mat,na=nastr,keepDiagonal=False):
	"""Replaces non-unique off-diagonal values in 2D matrix with na"""

	mat2 = [mat[i][:] for i in range(len(mat))]
	for i in range(len(mat2)):
		for j in range(i+int(keepDiagonal)-1):
		# for j in range(i+int(keepDiagonal)):
			mat2[i][j] = na
	return mat2

def matrixAverageByRow(mat):
	"""Returns a list of average per row of mat using unique values across columns."""
	return map(avg, replaceOffDiagonal(mat))

def matrixAverage(mat,useDiagonal=False):
	"""Returns the average of unique off-diagonal values in a 2D matrix."""
	return avg(uniqueEntries(mat,useDiagonal=useDiagonal))

def summarizeMatrix(mat,func=avg,useDiagonal=False):
	"""Returns the average of unique off-diagonal values in a 2D matrix."""
	return func(uniqueEntries(mat,useDiagonal=useDiagonal))


def matrixDifference(D,A,exp=2,W=[]):
	"""Least squares difference between two matrices"""
	import numpy
	n = D.shape[0]#range(D)
	m = D.shape[1]#range(D[0])
	
	if not W: W = [[1.0 for j in range(m)] for i in range(n)]
	
	x = 0.0
	for i in range(n):
		for j in range(m):
			# x += W[i][j]*(D[i][j] - A[i][j])**exp
			x += W[i][j]*abs(D[i,j] - A[i,j])**exp
	return x

def matDifference(D,A,exp=2,W=[]):
	return matrixDifference(D,A,exp,W)

def matrixAngles(D,A,exp=2,W=[]):
	"""Least squares difference between two matrices"""
	import numpy
	n = D.shape[0]#range(D)
	m = D.shape[1]#range(D[0])
	
	x = 0.0
	for i in range(m):
		y = angle(D[:,i],A[:,i])
		if y > 90.0: y = abs(180.0 - y)
		x += y
	return x/float(m)

def matAngles(D,A,exp=2,W=[]):
	return matrixAngles(D,A,exp,W)

# def standardize(X,center=True,scale=True):
# 	"""Mean center and scale by stdev each row in matrix X."""
# 	import numpy as np
# 	Y = np.array(X)
# 	for i in range(Y.shape[0]): 
# 		if center: Y[i,:] -= np.mean(Y[i,:])
# 		if scale: Y[i,:] /= np.std(Y[i,:])
# 	return Y
# 
def standardize(Y, func=lambda m,s: lambda x: divna((x - m), s)):
	"""mean center and scale expression data"""
	return map(func(avg(Y),sd(Y)), Y)


def retypeArray(X,type):
	"""Typecast every element in X to type."""
	X2 = X[:]
	for i in range(len(X)):
		X2[i] = map(type, X[i])
	return X2

def castemat(m, cast):
	return retypeArray(m,cast)

def mat2triples(x):
	triples = []
	for i in range(len(x)):
		for j in range(len(x[i])):
			triples += [ [i,j,x[i][j]] ]
	return triples



# GEOSTATISTICS
# -------------
def correlogram(M, dfunc='euclidean', r2=False, common=False, na=True):
	import numpy as np
	m = v = CI = None
	if common or dfunc=='mahalanobis':
		Mt = transpose(M)
		av = [avg(Mt[i]) for i in range(len(Mt))]
		m = avg(av)
		v = variance(av)
		CI = covmat(M) # compute covariance matrix
		#CI = covmat(M).I # compute covariance matrix
		print >> sys.stderr, 'finished computing C.I'
	print delim.join(['Pair','Cor', 'Dist'])
	corr = []; dist = []
	pairs = makeCombinations(len(M),2,start_idx=0)
	# TODO: replace ints with actual labels....
	for pair in pairs:
		i,j = pair[0],pair[1]
		if hasna(M[i]): print 'na!', M[i]
		if hasna(M[j]): print 'na!', M[j]
		
		c = r2 and correlation(M[i],M[j])**2 or correlation(M[i],M[j])
		d = distance(M[i],M[j], m,v,CI, dfunc=dfunc,na=na)
		#print np.nan, type(d), np.isnan(d)
		if c != nastr:# and not np.isnan(d):
			corr.append( c ); dist.append( d )
			print delim.join([str(pair),str(c),str(d)])
		
	return corr,dist


def hasna(v):
	#print 'bit vec', map(lambda x: x==nastr and 0 or 1, v)
	nacount = sum(map(lambda x: x==nastr and 0 or 1, v))
	if nacount < len(v): return True
	return False


# SORTING ROUTINES
# ----------------
def sortTableByCol(tab,idx=0):
	indexaskey = lambda k: lambda row: row[k]
	return sorted(tab, key=indexaskey(idx))

def quicksortTable(tab,idx=0):
	return sortTable(tab, idx)

def sortTable(tab, idx=0):
	"""Sort a homogeneous list of lists according to idx. the default
	column to sort on is the first column. since there is the chance of using
	large datasets, the inplace quicksort algorithm will be used. Warning: standard recursion limits apply."""
	quicksort(tab, 0, len(tab)-1, idx)

def quicksort(A, p, r, idx=0):
	if p < r:
		q = partition(A, p, r, idx)
		quicksort(A, p, q-1, idx)
		quicksort(A, q+1, r, idx)

def partition(A, p, r, idx=0):
	x = A[r][idx]
	i = p-1
	for j in range(p, r):
		if A[j][idx] <= x:
			i += 1
			temp = A[i]
			A[i] = A[j]
			A[j] = temp
	temp = A[i+1]
	A[i+1] = A[r]
	A[r] = temp
	return i+1


def closestValue(lst, val):
	"""use O(lgn) speed to find nearest value in list"""
	
	mind = nastr
	val = nastr
	for i in range(len(lst)):
		tmp = absna(subna(lst[i], val))
		if tmp != nastr and tmp < mind:
			mind = tmp; val = lst[i]
	return val


def sortVecKeysByValue(dct={}):
	"""Deprecated."""
	return keysSortedByValue(dct)

def keysSortedByValue(dct={}):
	items = dct.items()
	backitems = [[v[1], v[0]] for v in items]
	backitems.sort()
	return [backitems[i][1] for i in range(len(backitems))]

def sortVecByKeys(dct={}): 
	"""Deprecated."""
	return valuesSortedByKey(dct)

def valuesSortedByKey(dct={}):
	keys = dct.keys()
	keys.sort()
	return [dct[k] for k in keys]

def pairsSortedByKey(dct={}):
	"""Returns a vector of key:value pairs sorted by key."""
	items = dct.items()
	items.sort()
	return items

def pairsSortedByValue(dct={}):
	"""Returns a vector of key:value pairs sorted by value."""
	items = dct.items()
	backitems = [[v[1], v[0]] for v in items]
	backitems.sort()
	return [[backitems[i][1], backitems[i][0]] for i in range(len(backitems))]
	
	
	
	# keys = dct.keys()
	# 	values = dct.values()
	# 	svalues = values[:]
	# 	svalues.sort()
	# 	covered = [-1 for i in range(len(values))]
	# 	toreturn = []
	# 	for v in svalues:
	# 		idx = values.index(v)
	# 		covered[idx] += 1
	# 		useidx = idx+covered[idx]
	# 		#print 'v', v, 'idx', idx, 'use', useidx
	# 		toreturn.append([keys[useidx], v])
	# 	return toreturn




# RANDOM SELECTION
# ----------------
def weightedChoice(choices):
	total = sum(w for c, w in choices)
	r = random.uniform(0, total)
	upto = 0
	for c, w in choices:
		if upto + w > r:
			return c
		upto += w
	assert False, "Shouldn't get here"




# OTHER
# -----
def aic(n,k,e,correction=1):
	"""Akaike information content (AIC) takes number of samples (n), 
	number of parameters (k), and the residual squared error of the model.
	There is an optional correction factor, on by default."""
	n = float(n)
	e = float(e)
	k = float(k)
	# watch out for log(0), since error can be 0
	if e == 0: logterm = 0
	else: logterm = log(e/n)
	
	if correction:  return n*logterm+2*k+(2*k*(k+1)/(n-k-1))
	else:			return n*logterm+2*k


def enlist(lst):
	"""Wrap each item in a list in square braces."""
	return [ [lst[i]] for i in range(len(lst)) ]


def fixhyphen(astr):
		ret = astr
		if re.match('^(Y......)(\S)$',astr):
			return astr[0:7]+'-'+astr[7:len(astr)]
		
		return astr

def aliasGeneNames(names, filename='/Users/simola/Documents/UPenn/Thesis/yeast_genomes/annotations/registry.genenames.tab', invert=False, colidx=4):
	# use short name labels
	# load the registry file, then pull the alternate names...probably did this before somewhere
	import sio as io, re
	d,aliases,c = io.readMatrix(filename)
	altmp = map(lambda x: x.lower(), aliases)
	
	akas = vslice(d,0)
	akatmp = map(lambda x: x.lower(), akas)
	if not invert:
		SGDids = vslice(d, colidx)
		
		# add hyphen to names lacking one
		names = map(fixhyphen,names)
		Gnames = []
		for gene in names:
			try: Gnames += [aliases[SGDids.index(gene)]]
			except IndexError: Gnames += [gene]
			except ValueError: Gnames += [gene]
		return Gnames
	else:
		# back-map from aliases to names
		# alias is stored as rowname
		Gnames = []
		for alias in map(lambda x: x.lower(), names):
			# print 'alias', alias
			idx = None
			try:
				idx = altmp.index(alias)
			except ValueError:
				try:
					idx = akatmp.index(alias)
				except ValueError: pass
			sysid = alias
			if idx != None:
				sysid = d[idx][colidx]
			Gnames += [sysid]
		return Gnames
		
		


def getInfo(names, filename='/Users/simola/Documents/UPenn/Thesis/yeast_genomes/annotations/registry.genenames.tab', colidx=4):
	"""Returns a dictionary of information for each gene (name, description, etc)"""
	import sio as io
	d,aliases,c = io.readTable(filename)
	SGDids = vslice(d, colidx)
	info = {}
	
	for gene in names:
		info[gene] = {}
		try:
			idx = SGDids.index(gene)
			alias = aliases[idx]
			rowinfo = d[idx]
			info[gene]['name'] = alias
			for k in range(len(c)):
				info[gene][c[k]] = rowinfo[k]
		except IndexError: pass
		except ValueError: pass
	return info


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

def getRanks(master,query):
	ranks = []
	for q in query:
		ranks += [master.index(q)]
	return ranks
	

nafloat = lambda x: (x==nastr or float(x))
asnumeric = lambda A: map(lambda row: map(nafloat, row), A)
eithereq = lambda s: lambda pair: (pair[0]!=s and pair[1]!=s and True)
tripleeq = lambda s: lambda triple: (triple[0]!=s and triple[1]!=s and triple[2]!=s and True)
filterstr = lambda s: lambda v: filter(eq(s), v)
pairedfilterstr = lambda s: lambda v: filter(eithereq(s), v)
triplefilterstr = lambda s: lambda v: filter(tripleeq(s), v)

def listsub(X,Y):
	X2 = X[:]
	for y in Y:
		try: X2.pop(X2.index(y))
		except KeyError: continue
	return X2

affine = lambda b, a: lambda x: b + a*x

def dictAdd(dct1, dct2, overwrite=False):
	# shallow copy dct1
	fd = {}
	for k in dct1.keys(): fd[k] = dct1[k]
	for k2 in dct2.keys():
		if overwrite:
			fd[k2] = dct2[k2]
		else:
			if k2 not in dct1: fd[k2] = dct2[k2]
	return fd


def argmaxna(v):
	# NOT TESTED
	"""Return argmax of a list"""
	mi = 0
	i = 0
	while v[i] == nastr: i += 1
	mv = v[i]
	
	for i in range(i,len(v)):
		if v[i] != nastr and float(v[i]) > float(mv):
			mv = float(v[i])
			mi = i
	return mi

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
	"""Return argmin of a list"""
	mi = 0
	mv = v[0]
	for i in range(1,len(v)):
		if float(v[i]) < float(mv):
			mv = float(v[i])
			mi = i
	return mi

def argmaxx(v):
	"""Return argmax of a list. Returns a list since multiple maxs are possible."""
	mi = []
	mv = -GLOB_INF
	for i in range(len(v)):
		tv = float(v[i])
		if tv > mv: mv = tv; mi = [i]
		elif tv == mv: mi.append(i)
	return mi



# FORMAT CONVERSION
# ------------------
def pwmscan2gff():
	return 1





# MULTIPLE HYPOTHESIS TESTING
# ---------------------------
def BH(P, L=[], alpha=0.05, dependency='independent', invertPs=False, verbose=False, returncutoff=False, returnP=False):
	"""Given a list of P-values and corresponding labels, returns 
	those hypotheses (Labels) significant at FDR <= alpha using 
	a Benjamini-Hochberg procedure. This involves finding the 
	largest i such that p(i) <= alpha i/m*cm, i in [1:m].
	By default cm = 1, which applies for independent or positively
	correlated p-values
	"""
	if L:
		assert(len(P)==len(L)),\
	      'P-values and Labels must have same length.'
	else: L = range(len(P))
	
	H = pairedfilterstr(nastr)(zip(P,L))
	if invertPs:
		H = [[1-float(H[i][0]),H[i][1]] for i in range(len(H))]
		alpha = 1-alpha
		
	H.sort()
	m = len(H)
	
	cm = 1 # independence of hypotheses
	if dependency == 'positive' or dependency == '+': 
		cm = sum([1/float(j) for j in xrange(1,m+1)])
	elif dependency == 'negative' or dependency == '-': 
		return 'negative dependency not yet implemented'
	
	if verbose: 
		print 'BH: n=%s hypotheses; Dependency=%s c(m)=%s'%(m, dependency, cm)
	
	
	# if alpha < .5: H.reverse() # decreasing pvalues
	if alpha <= .5:
		pvalcut = 1
		
		# top -down
		# i = m
		# while i > 0:
		# 	pvalue = float(H[i-1][0])
		# 	pvalcut = alpha*i/(float(m*cm))
		# 	if verbose: 
		# 		sys.stdout.write('i='+str(i)+': '+str(pvalue)+' < '+str(pvalcut))
		# 		sys.stdout.flush()
		# 	if pvalue <= pvalcut: 
		# 		if verbose: print ' ==> SUCCESS'
		# 		break
		# 	else:
		# 		if verbose: print
		# 	i -= 1
		
		# bottom up
		i = 1
		while i < m:
			pvalue = float(H[i-1][0])
			pvalcut = alpha*i/(float(m*cm))
			if verbose: 
				sys.stdout.write('i='+str(i)+': '+str(pvalue)+' < '+str(pvalcut))
				sys.stdout.flush()
			if pvalue > pvalcut: 
				if verbose: print ' ==> SUCCESS'
				i -= 1
				break
			else:
				if verbose: print
			i += 1
		
	else:
		print 'NEED TO CONVERT PVALUES TO FLOAT'
		K = [1-H[i][0] for i in range(len(H))]
		H = [[K[i],H[i][1]] for i in range(len(H))]
		H.reverse()
		pvalcut = 0
		i = m
		while i > 0:
			pvalue = float(H[i-1][0])
			pvalcut = (1-alpha)*i/float(m*cm)
			if verbose: 
				sys.stdout.write('i='+str(i)+': '+str(1-pvalue)+' > '+str(1-pvalcut))
				sys.stdout.flush()
			if pvalue < pvalcut: 
				if verbose: print ' ==> SUCCESS'
				break
			else:
				if verbose: print
			i -= 1
	
	if verbose:
		if alpha <= .5:
			print 'BH: '+str(i)+'/'+str(m)+' significant hypotheses below FDR alpha', alpha
		else:
			print 'BH: '+str(i)+'/'+str(m)+' significant hypotheses above FDR alpha', 1-alpha
	
	if not returncutoff:
		
		if i==0: return []
		
		if not returnP:
			return [H[n][1] for n in range(len(H))][0:i]
		else:
			return [[H[n][0],H[n][1]] for n in range(len(H))][0:i]
	else:
		return [H[n][1] for n in range(len(H))][0:i], pvalcut

def Bonferroni(P, L=[], alpha=0.05, verbose=False, n=-1, returnP=False):
	if L:
		assert(len(P)==len(L)),\
	      'P-values and Labels must have same length.'
	else: L = range(len(P))
	
	H = [[P[i],L[i]] for i in range(len(P))]
	H = pairedfilterstr(nastr)(H)
	
	m = len(H)
	if n > -1: m = n
	
	if m <=0.0: return []
	BonAlpha = alpha/float(m)
	if alpha > 0.5:
		BonAlpha = (1.0-alpha)/float(m)
		BonAlpha = 1.0-BonAlpha
	# which hypotheses have pvalue < BonAlpha?
	sigPairs = []
	if alpha <= .5:
		sigPairs = filter(lambda x: (float(x[0]) < BonAlpha and True), H)
	else:
		sigPairs = filter(lambda x: (float(x[0]) > BonAlpha and True), H)
	
	if verbose: 
		if alpha <= 0.5:
			print 'Bonferroni:', str(len(sigPairs))+'/'+str(m)+' < FWER alpha = %f' % (BonAlpha)
		else:
			print 'Bonferroni:', str(len(sigPairs))+'/'+str(m)+' > FWER alpha = %f' % (BonAlpha)
	
	if returnP: return sigPairs
	else: return [sigPairs[i][1] for i in range(len(sigPairs))]

def Pcutoff(P,L=[], alpha=.05, returnP=False):
	if L:
		assert(len(P)==len(L)),\
	      'P-values and Labels must have same length.'
	else: L = range(len(P))
	
	H = pairedfilterstr(nastr)(zip(P,L))
	m = len(H)
	
	Pf, Lf = unzip(H)
	sigpairs = filter(lambda x: x[0] < alpha, H)
	# print 'sigpairs', sigpairs
	if len(sigpairs):
		if returnP:
			return sigpairs
		else:
			sigPs, sigLs = unzip(sigpairs)
			return sigLs
	else: return []


def GOenrichment(names, notnames=[], pool='life cycle', dirfile='', verbose=False, pipe=sys.stdout, col=0, useAlias=False):
	"""
	Given a list of gene names, perform a GO enrichment analysis of those names in each group, 
	provided as separate files in a directory dirfile. Each file contains a list of names that are 
	newline delimited.
	"""
	
	if len(dirfile) == 0:
		if pool == 'life cycle':
			dirfile = '/Users/simola/Documents/UPenn/Thesis/Microarrays/yeast/exprVarGroupsSansAll/'
		elif pool == 'hires':
			dirfile = '/Users/simola/Documents/UPenn/Thesis/Microarrays/yeast/GeneGroups/'
		else:
			raise ValueError('Either dirfile or pool arguments must be specified.')
	
	if len(notnames)==0:
		# use genomic background
		import sio as io
		notnames = io.readList('/Users/simola/Documents/UPenn/Thesis/Microarrays/yeast/gene_names.txt')
	
	if useAlias: 
		names = aliasGeneNames(names)
		notnames = aliasGeneNames(notnames)
	
	
	col -= 1 # since we just parsed the first column
	dat = []
	
	import sio as io
	files = io.getFiles(dirfile)
	# split into groups
	results = {}
	metathegenes = []
	metagroup = []
	for f in files:
		parts = f.split('/')
		group = parts[len(parts)-1]
		group = group[0:len(group)-4]
		
		results[group] = {'pvalue':nastr, 'proportion':nastr, 'genes':[], 'count':0}
		
		(contents, groupNames, crap) = io.readMatrix(f, header=False, comments=False)
		
		if col > -1:
			# then parse names from a different column of the file
			# print 'group', f
			for i in range(len(contents)):
				if len(contents[i]) <= 9: print i, contents[i]
			groupNames = vcol(col, contents)
			# split and take last itemm - for AmiGO files
			for i in range(len(groupNames)):
				parts = split(groupNames[i], '|')
				groupNames[i] = parts[len(parts)-1]
		
		# print 'groupNames'
		# print groupNames
		
		# hash the names for O(1) search
		fast = {}
		for n in groupNames: fast[n] = 1
		
		thegenes = []
		groupnames = 0 # a
		for i in range(len(names)):
			if names[i] in fast: 
				groupnames += 1
				thegenes += [names[i]]
			
		groupnotnames = 0 # b
		for i in range(len(notnames)):
			if notnames[i] in fast: groupnotnames += 1
			
		notgroupnames = len(names) - groupnames # c
		notgroupnotnames = len(notnames) - groupnotnames # d
		
		# positive association with group and names
		
		dat += [[[groupnames, groupnotnames], [notgroupnames, notgroupnotnames]]]
		metathegenes += [thegenes]
		metagroup += [group]
	
	pvalues = batchFisherExactTest(dat, tails=1)
	# print 'inside p', pvalues
	for i in range(len(files)):
		group = metagroup[i]
		(groupnames, groupnotnames), (notgroupnames, notgroupnotnames) = dat[i]
		thegenes = metathegenes[i]
		pvalue = pvalues[i]
		results[group]['pvalue'] = pvalue
		
		try: results[group]['proportion'] = groupnames/float(len(names))
		except ZeroDivisionError: results[group]['proportion'] = -1
		results[group]['count'] = groupnames
		results[group]['genes'] = thegenes
		if verbose: print >>pipe, 'GO enrichment on', group
		if verbose:  print >>pipe, 'Gene group:', groupnames, 'Not gene group', groupnotnames
		if verbose:  print >>pipe, 'Gene not in group', notgroupnames, 'Not genes not in group', notgroupnotnames
		if verbose: print >>pipe, 'Fisher\'s Exact test p-value', pvalue
	
	return results


def antGOenrichment(names, notnames=[], gofile='cflo.v3.3.UniProt.GO.txt', dirfile='/Users/simola/Documents/UPenn/BergerLab/data/', verbose=False):
	import sio as io
	# first parse the GO table into a dictionary of GO terms returning genes for each
	dct = {}
	allgenenames = []
	d,r,c = io.readTable(dirfile+gofile, header=False)
	for i in range(len(d)):
		gene = r[i].split('--')[0]
		terms = d[i]
		allgenenames += [gene]
		for gobar in terms:
			for go in gobar.split('|'):
				if len(go)>1: # control for null and \n and \r stray lines
					try: dct[go] += [gene]
					except KeyError: dct[go] = [gene]
	
	# use genomic background?
	if not len(notnames): notnames = allgenenames
	print 'Foreground genes:', len(names), names[:10]
	print 'Background genes:', len(notnames), notnames[:10]
	
	# process enrichment for each GO term
	dat = []
	results = {}
	metathegenes = []
	metagroup = []
	GOterms = sorted(dct.keys())
	
	for group in GOterms:
		results[group] = {'pvalue':nastr, 'proportion':nastr, 'gw-proportion':nastr, 'genes':[], 'count':0, 'group_count':len(dct[group]), 'table':nastr}
		
		# hash the names for O(1) search
		fast = dict( map(lambda x: (x,None), dct[group]) )
		# fast = dict( map(lambda x: (x,None), names) )
		
		# print 'test', names[:50]
		# print 'dict', fast
		
		# unique list of names
		namesU = unique(names)
		notnamesU = unique(notnames)
		# a: genes of interest in group of interest
		thegenes = []
		groupnames = 0 # a
		for i in range(len(namesU)):
			if namesU[i] in fast: 
				groupnames += 1
				thegenes += [namesU[i]]
			# 	print 'OK', namesU[i]
			# else:
			# 	print 'namesU no in dict', namesU[i], fast.items()[:20]
				
		# c: genes of interest not in group of interest
		notgroupnames = len(namesU) - groupnames
		
		# b: background genes in the group of interest
		groupnotnames = 0 # b
		for i in range(len(notnamesU)):
			if notnamesU[i] in fast: groupnotnames += 1
		
		
		# groupnotnames -= groupnames # subtract the genes we are interested in
		
		# d: background genes not in group of interest
		notgroupnotnames = len(notnamesU) - groupnotnames
		
		# positive association with group and names
		dat += [[[groupnames, groupnotnames], [notgroupnames, notgroupnotnames]]]
		metathegenes += [thegenes]
		metagroup += [group]
	
	pvalues = batchFisherExactTest(dat, tails=1)
	
	if not(len(pvalues)): return results
	
	strs = []
	
	i = 0
	for group in GOterms:
		(groupnames, groupnotnames), (notgroupnames, notgroupnotnames) = dat[i]
		thegenes = metathegenes[i]
		pvalue = pvalues[i]
		
		results[group]['count'] = groupnames
		results[group]['genes'] = thegenes
		results[group]['pvalue'] = pvalue
		results[group]['table'] = dat[i]
		
		if groupnames == 0: pvalue = nastr
		
		i+=1
		
		try: results[group]['gw-proportion'] = groupnames/float(len(namesU))
		except ZeroDivisionError: results[group]['gw-proportion'] = nastr
		
		try: results[group]['proportion'] = groupnames/float(groupnames+groupnotnames)
		except ZeroDivisionError: results[group]['proportion'] = nastr
		
		astr = 'GO enrichment on %s\n'%(group)
		astr += 'Gene group=%s, Not gene group=%s\n'%(groupnames,groupnotnames)
		astr += 'Gene not in group=%s, Not genes not in group=%s\n'%(notgroupnames, notgroupnotnames)
		astr += 'Fisher\'s Exact test p-value %s\n'%(pvalue)
		strs += [astr]
	
	if verbose:
		pairs = filter(lambda x: x[0] <= 0.15 and True, zip(pvalues,strs))
		for pval,info in sorted(pairs): print info
	
	return results


def GOsig(results, alpha=.05, type='BH', verbose=False, invertPs=False):
	terms = results.keys()
	pvals = [results[terms[i]]['pvalue'] for i in range(len(terms))]
	# props = [results[terms[i]]['proportion'] for i in range(len(terms))]
	
	if type == 'BH': # Beniamini-Hochberg
		return BH(pvals,terms,alpha=alpha, verbose=verbose, invertPs=invertPs)
	elif type == 'Bonferroni': # Bonferroni
		return Bonferroni(pvals,terms,alpha=alpha, verbose=verbose)
	else:
		if invertPs: pvals = map(lambda x: subna(1,x), pvals)
		return Pcutoff(pvals,terms,alpha=alpha)


def parseGOobo(filename):
	fh = open(filename)
	dat = fh.readlines()
	fh.close()
	
	dct = {}
	
	# erase the header
	header = []
	begin = 0
	for i in range(len(dat)):
		# print dat[i]
		if dat[i] != '[Term]\n': 
			# print 'CONTINUE'
			header += [dat[i]]
		else: 
			# print 'OK!'
			break
		begin = i
		
	# print 'header'
	# print header
	
	dat = '\n'.join(dat[begin:])
	
	for group in dat.split('[Term]\n'):
		info = group.split('\n')
		
		tmp = {'id':'NA', 'name':'NA', 'namespace':'NA', 'def':'NA', 'synonym':'NA', 'is_a':[]}
		
		for line in info:
			for astr in tmp.keys():
				if re.match('^%s:.+'%(astr), line):
					crap, keep, crap = re.split('^%s:(.+)'%(astr), line)
					if astr == 'is_a':
						tmp[astr] += [keep.strip()]
					else:
						tmp[astr] = keep.strip()
						
		if tmp['id'] != 'NA': dct[tmp['id']] = tmp
		# print 'test', tmp['id'], tmp
		
	return dct


def fullGOanalysis(genelst, bglst=[], alpha=.05, method='BH', name='sample', savedir='', godir='/Users/simola/Documents/UPenn/BergerLab/data/', gofile='cflo.v3.3.UniProt.GO.txt', gotable='/Users/simola/Documents/UPenn/BergerLab/data/gene_ontology_ext.obo.txt', gononsig=10, keywords=[],plottop=25):
	
	# method is Bonferroni or BH
	
	import sio
	sio.createdir(savedir)
	
	sio.printList(genelst,file=savedir+'GO.query.genes.txt')
	sio.printList(['# IF empty then ALL genes were used for background']+bglst,file=savedir+'GO.background.genes.txt')
	
	if not len(genelst): 
		print 'Warning: no genes for GO analysis, name=%s'%(name)
		return
	
	akago = parseGOobo(gotable)
	
	res = antGOenrichment(genelst, bglst, gofile=gofile,  dirfile=godir, verbose=0)
	ngoterms = len([res[k]['pvalue'] for k in res.keys()])
	print 'ngoterms', ngoterms
	
	sig = GOsig(res, alpha=alpha, type=method, invertPs=0)
	ps = [res[s]['pvalue'] for s in sig] # make the large values small
	qs = [res[s]['proportion'] for s in sig]
	q2s = [res[s]['gw-proportion'] for s in sig]
	cs = [res[s]['count'] for s in sig]
	sigg = [res[s]['genes'] for s in sig]
	
	# save list of all genes that are GO associated
	flatgogenes = []
	for k in res.keys():
		flatgogenes += res[k]['genes']
	flatgogenes.sort()
	sio.printList(unique(flatgogenes), file=savedir+'all.GO.assoc.genes.txt')
	
	# save list of all genes that are sig GO associated
	flatgogenes = []
	for lst in sigg: flatgogenes += lst
	flatgogenes.sort()
	sio.printList(unique(flatgogenes), file=savedir+'sig.GO.assoc.genes.txt')
	
	# save the significant genes to a directory
	gogenedir = savedir+'genes4sigTerms/'
	sio.createdir(gogenedir)
	
	for s,g in zip(sig,sigg):
		sio.printList(g,file=gogenedir+'%s_genes.txt'%(s))
	
	sigfile = savedir+'%s_GO_%s_%s.txt'%(name,method,alpha)
	gof = open(sigfile,'w')
	
	print '\nGO enrichment(fg=%s, bg=%s, GOterms=%s, %s alpha < %s) => %s sig. terms'%(len(genelst), len(bglst), ngoterms,method,alpha, len(sig))
	print >>gof, 'GO enrichment(fg=%s, bg=%s, GOterms=%s, %s alpha < %s) => %s sig. terms\n\n'%(len(genelst), len(bglst), ngoterms,method,alpha, len(sig))
	print >>gof, '\n\nSignificant terms:'
	print >>gof, '\t'.join(['Name', 'GOID', 'P', 'gw-Q', 'Q', 'Count'])
	
	# replace significant terms with their actual names
	keep = []
	for s,p,q,q2,c in zip(sig,ps,qs,q2s,cs):
		if s in akago: keep += [(p,[akago[s]['name'],s,p,q,q2,c])]
		else: keep += [(p,['NA',s,p,q,q2,c])]
	keep.sort() # most to least significant
	sio.printTable([s for p,s in keep], pipe=gof)
	# sio.printList(['%s\t%s\t%s\t%.3e\t%.3f\t%s'%(s) for p,s in keep], pipe=gof)
	
	# sort all terms by pvalue
	allsort = []
	for k in res.iterkeys():
		allsort += [(res[k]['pvalue'],k)]
	allsort.sort()
	
	# save genes from top 10 categories
	nogogenedir = savedir+'genes4nonsigTerms/'
	sio.createdir(nogogenedir)
	for p,k in allsort[:gononsig]:
		sio.printList(res[k]['genes'],file=nogogenedir+'%s_genes.txt'%(k))
	
	# return list of genes having go term words user specified
	keywordfile = savedir+'keywords-%s.txt'%(','.join(keywords))
	# sio.createdie(keyworddir)
	keywordgenes = {}
	
	
	table = []
	plotdata = [] # record GO label and -log10pvalue
	plotlabels = []
	order = ['pvalue', 'count', 'group_count', 'proportion', 'gw-proportion', 'table']
	print >> gof, '\n\n\n'
	print >>gof, 'All GO terms:'
	for p,k in allsort: 
		row = []
		lab = ''
		if k in akago: lab = akago[k]['name']
		if lab == '': lab = k
		print >> gof, '%s (%s): '%(k,lab),
		row += ['%s (%s): '%(k,lab)]
		
		for kk in keywords:
			if kk in lab:
				for kg in res[k]['genes']:
					if kg in keywordgenes:
						keywordgenes[kg] += ['%s:%s'%(k,lab)]
					else:
						keywordgenes[kg] = ['\t%s:%s'%(k,lab)]
		
		ct = 0
		for kk in order:
			if kk != 'genes':
				if ct == 0: 
					gof.write('%s=%s'%(kk,res[k][kk]))
					row += ['%s=%s'%(kk,res[k][kk])]
				else: 
					gof.write(', %s=%s'%(kk,res[k][kk]))
					row += [', %s=%s'%(kk,res[k][kk])]
				ct += 1
		print >> gof
		table += [row]
		
		# print 'label', lab, 'data', -log(p)
		plotdata += [[-log(p)]]
		plotlabels += [lab]
	
	kgd = [[kg,'\t'.join(keywordgenes[kg])] for kg in sorted(keywordgenes.keys())]
	sio.printTable(kgd,['Gene','Alias','GO terms'],file=keywordfile)
	
	if len(plotdata):
		import plot
		plot.chart(plotdata[:plottop],xlabels=plotlabels[:plottop],file=savedir+'GOsigterms.logp.chart.pdf', ylabel='-log_{10}P-value', fontsize=6, ratio=2,rotateLabels=90)
	
	# return table of significant terms
	return table


def loadOrthoDB(orthofile='/Volumes/antbox/7ant_comparison/work/AntOrthoDB/Formicidae7.OGs.txt',col=1):
	import sio as io
	d,r,c = io.readTable(orthofile)
	
	# map gene name to ogid
	g2og = {}
	og2g = {}
	for i in range(len(d)):
		ogid = r[i]
		target = d[i][col]
		g2og[target] = ogid
		
		if ogid in og2g: og2g[ogid] += [target]
		else: og2g[ogid] = [target]
		
	print >> sys.stderr, '%s ogids'%(len(g2og))
	
	return g2og,og2g


def mapToOrthoDB():
	return 'not yet implemented'


def getGeneForOrthoid(lst,species=None,orthofile=None,orthodb=None):
	if orthofile:
		orthodb = loadOrthoDB(orthofile)[1]
	elif not orthodb:
		orthodb = loadOrthoDB()[1]
	
	ENLISTED = False
	if type(lst) != type([]):
		ENLISTED = True
		lst = [lst]
	
	# list of orthodb ids
	out = []
	for ogid in lst:
		# if species return ids for just that species, otherwise return all entries
		if ogid in orthodb:
			if species:
				out += [filter(lambda x: species in x.lower(), orthodb[ogid])]
			else:
				out += [orthodb[ogid]]
			
		else:
			out += [ogid]
	
	if ENLISTED: return out[0]
	else: return flatten(out)



#
def translateorthodb(lst=None, query='CFLOR', target='OrthoID', dbfile='/Users/simola/Documents/UPenn/BergerLab/data/insect.ortho.map.ODB85.txt', groupparalogs=False, chooseparalog=False):
	import sio as io
	# load the orthodb
	L,foo,Lspecies = io.readTable(dbfile,rownames=0,header=1)
	
	targetidx = filter(lambda x: x[1].lower().count(target.lower()), zip(range(len(Lspecies)), Lspecies))[0][0]
	queryidx = filter(lambda x: x[1].lower().count(query.lower()), zip(range(len(Lspecies)), Lspecies))[0][0]
	
	# targetidx = Lspecies.index(target)
	# queryidx = Lspecies.index(query)
	# print 'query=%s, target=%s'%(queryidx, targetidx)
	# print 'OK', L[0][queryidx], L[0][targetidx]
	
	om = {}
	for row in L:
		if chooseparalog:
			for xx in row[queryidx].split(','): om[xx] = row[targetidx].split(',')[0]
		elif groupparalogs:
			for xx in row[queryidx].split(','): om[xx] = row[targetidx]
		else:
			for xx in row[queryidx].split(','): om[xx] = row[targetidx].split(',')
	
	def lookup(x):
		# print 'test', '"%s"'%(x), x in om and '"%s"'%(om[x]) or -1
		if x in om: return om[x]
		return nastr
	
	if lst == None: lst = om.keys()
	
	return flatten(map(lookup, lst))


def samplesize(x, na=nastr):
	"""Returns the value n which is the number of non-NA data points in a list."""
	return len(filterstr(na)(x))

def matsamplesize(x):
	"""Returns the number of columns in a 2d python array."""
	return len(x[0])


def gaussian(x, m, sd):
	"""Return probability of x under gaussian distribution, mean m, std sd."""
	p = (1./float(sd*math.sqrt(2.*math.pi)))*math.e**(-.5*((x-m)/sd)**2.)
	return p

def BetaStirling(x,y):
	"""Based on Stirling's approximation"""
	
	numer = float(x**(x-0.5)*y**(y-0.5))
	denom = float((x+y)**(x+y-0.5))
	return math.sqrt(2*math.pi)*numer/denom

def Beta(x,y, dt=.0000041):
	"""Integral approximation"""
	beta = 0.0; t = dt
	while t < 1.0: 
		beta += t**(x-1)*(1-t)**(y-1)*dt
		t += dt
	return float(beta)

def incompleteBeta(x,a,b, dt=.0000041):
	"""Approximation to integral based on dt"""
	assert(a>0 and b>0),\
	  'DOF error: a > 0 and b > 0.'
	beta = 0.0; t = dt
	while t < x: 
		beta += t**(a-1)*(1-t)**(b-1)*dt
		t += dt
	return float(beta)

def regularizedIncompleteBetaForIntegers(x,a,b):
	"""Ratio of incomplete beta and complete beta functions for integers a,b.
	
	The incomplete beta function is a generalization of the beta function that 
	replaces the definite integral of the beta function with an indefinite 
	integral. The situation is analogous to the incomplete gamma function being 
	a generalization of the gamma function.
	"""
	
	a = int(a)
	b = int(b)
	eye = 0
	for j in range(a, a+b, 1):
		f = factorial(a+b-1)/(factorial(j)*factorial(a+b-1-j))
		eye += f * x**j * (1-x)**(a+b-1-j)
	return eye

def regularizedIncompleteBeta(x,a,b, dt=1e-4):
	return incompleteBeta(x,a,b, dt=dt)/Beta(a,b)


def BetaStirling(a,b):
	"""Stirling's Beta approximation"""
	
	import math
	
	c = math.sqrt(2*math.pi)
	n = a**(a-.5)*b**(b-.5)
	d = (a+b)**(a+b-.5)
	return c*(n/d)

def TajimaDistPDF(D, n):
	"""Beta distribution PDF following, Tajima (1989).
	User supplies estimate of Tajima D as D, and number
	of taxa used as n
	"""
	import maths
	if n < 2: return nastr
	
	a1 = float(sum(map(inv(1), range(1, n))))
	a2 = float(sum(map(inv(2), range(1, n))))
	b2 = (2*(n**2 + n + 3))/float(9*n*(n - 1))
	c2 = b2 - (n+2)/(a1*n) + a2/a1**2
	e2 = c2/(a1**2 + a2)
	sqrte2 = math.sqrt(e2)
	
	if sqrte2 == 0: return nastr
	
	Dmin = (2/float(n) - 1/a1)/sqrte2
	
	Dmax = nastr
	if n%2 == 0: # even
		Dmax = ( n/float(2*(n-1)) - 1/a1 )/sqrte2
	else:
		Dmax = ( (n+1)/float(2*n) - 1/a1 )/sqrte2
	
	a = Dmin
	b = Dmax
	
	alpha = -(1 + a*b)*b / float(b - a)
	beta =   (1 + a*b)*a / float(b - a)
	
	# this part is not very sensitive
	if nastr in [alpha, beta, b, a]: return nastr
	else:
		phiD = 0.0 # if we've gotten here and can't compute pvalue, it's
		# because it is 0
		try:
			# import scipy.stats as ss
			
			numer = (b - D)**(alpha-1) * (D - a)**(beta-1)
			denom = float((b - a)**(alpha+beta-1))
			# phiD = 1/BetaStirling(alpha,beta) * (numer/denom)
			phiD = 1/maths.beta(alpha,beta) * (numer/denom)
			# phiD = 1/ss.beta(alpha,beta) * (numer/denom)
		except ValueError: pass
		except TypeError: pass
		return phiD


# def BetaDistCDF(x, a, b, dt=1e-4):
# 	"""Ratio of Bx(a,b)/B(a,b). Same as FdistCDF, but
# 	with different parameters.
# 	"""
# 	return regularizedIncompleteBeta(x=x, a=alpha, b=beta, dt=dt)
# 
def BetaDistCDF(x, a, b):
	"""Ratio of Bx(a,b)/B(a,b). Same as FdistCDF, but
	with different parameters.
	"""
	import maths
	return maths.incompleteBeta(x,a,b)



# def FdistCDF(f, df1, df2, dt=1e-4):
# 	if f == nastr: return nastr
# 	x = df1*f/float(df1*f+df2)
# 	return regularizedIncompleteBeta(x=x, a=df1/2.0, b=df2/2.0, dt=dt)
# 
# def FdistCDFApprox(f,df1,df2):
# 	a = int(math.floor(df1/2.0))
# 	b = int(math.floor(df2/2.0))
# 	tdf1 = 2*a; tdf2 = 2*b
# 	x = tdf1*f/float(tdf1*f+tdf2)
# 	pvalue = regularizedIncompleteBetaForIntegers(x=x, a=a, b=b)
# 	if tdf1 == df1 and tdf2 == df2: return pvalue
# 	else:
# 		# approximate function by averaging
# 		if tdf1!=df1:a+=1
# 		if tdf2!=df2:b+=1
# 		tdf1 = 2*a; tdf2 = 2*b
# 		x = tdf1*f/float(tdf1*f+tdf2)
# 		pvalue2 = regularizedIncompleteBetaForIntegers(x=x, a=a, b=b)
# 		return (pvalue+pvalue2)/2.0
# 	# print 'f', f, 'x', x, 'df1',df1, 'df2',df2, 'pval', pvalue
# 
# 



def studentT(lst1=[], lst0=[], mu1=0, mu=0, sderr=None, equalVariance=False, verbose=False):
	"""Return the t-statistic for a student's t distribution.
	Can calculate t-statistic for one or two sample, unpaired.
	"""
	
	# count
	n = len(filterstr(nastr)(lst1))
	n0 = len(filterstr(nastr)(lst0))
	
	if n == 1 and n0 > 1:
		# reverse this or code will error
		n = n0
		n0 = 1
		tmp = lst1[:]
		lst1 = lst0
		lst0 = tmp
		tmp = mu
		mu = mu1
		mu1 = tmp

	if mu==0 and mu1==0 and n==1:
		mu1 = lst1[0]
	elif mu==0 and mu1==0 and n0==1:
		mu = lst0[0]

	if not mu1:
		try:
			# mean
			mu1 = avg(lst1)
			mu0 = mu
			if len(lst0): mu0 = avg(lst0)
			# stdev
			v1 = variance(lst1); s1 = math.sqrt(v1)
			v0 = nastr; s0 = nastr
			
			if len(lst0):
				v0 = variance(lst0); s0 = math.sqrt(v0)
		except TypeError: 
			pass
			# return nastr, nastr
	
		# one sample statistic
		s10 = s1
		v10 = nastr
		df = n-1
		t = nastr
		denom = nastr
	else:
		n0 = 1
		s10 = 0
		mu0 = mu
		df = n-1
	
	# two sample equal variance statistic
	if n0 > 1: 
		if n <= 30 or n0 <= 30:
			v10 = ((n-1)*s1**2 + (n0-1)*s0**2)/float(n + n0 - 2)
			denom = math.sqrt(v10*(1/float(n) + 1/float(n0)))
		else:
			v10 = math.sqrt( v1/float(n) + v0/float(n0) )
			denom = math.sqrt(v10)
		
		t = (mu1-mu0)/denom
		
		if equalVariance:
			df = n+n0-2
		else:
			# use the Welch-Satterthwaite equation to estimate dof
			num = ( v1/float(n) + v0/float(n0) )**2
			den = v1**2/(n**2*(n-1)) + v0**2/(n0**2*(n0-1))
			df = num/float(den)
			# print 'num/den/dof', num, den, df
		
	else:
		# print 'n',n,s10,sderr,mu1,mu0
		# one sample statistic
		denom = (s10/math.sqrt(n))
		if sderr: denom = sderr
		t = divna(mu1 - mu0, denom)
	
	if verbose:
		une = 1
		if len(lst0): une = 2
		print 't-test: '+str(une)+' sample '+str(equalVariance)+' equal variances'
		print 'mu1=%f, mu0=%f' % (mu1,mu0)
		print 'var1=%s, var0=%s' % (v1,v0)
		print 'obs1=%d, obs0=%d' % (n, n0)
		print 'pooled var=%s' % (v10)
		print 'pooled sd=%s' % (denom)
		print 'dof=%.3f' % (df)
		print 't=%s' % (t)
	return t,df

# Ttest = studentT
def TdistCDF(t, df):
	if t == nastr or df == nastr: return nastr
	import maths
	return maths.tdistCDF(t,df)

def TdistPDF(t,df,tails=2):
	"""Only 1 and 2 tailed p-values."""
	if t == nastr or df == nastr: return nastr
	
	abst = abs(t)
	pval = TdistCDF(abst,df)
	if t >= 0: pval = 1 - pval
	if t < 0: pval = 1 - pval
	
	if tails == 2: return pval
	else: return pval/2.0

def Ttest(lst1=[], lst0=[], mu1=0, mu=0, sderr=None, equalVariance=False, verbose=False, tails=2):
	t,df = studentT(lst1,lst0,mu1,mu,sderr,equalVariance,verbose)
	p = TdistPDF(t,df,tails)
	return {'df':df,'t':t, 'T':t, 'P':p, 'tails':tails,'eqvar':equalVariance}

def pairedTtest(lst1, lst2, verbose=False,tails=2):
	pairs = pairedfilterstr(nastr)(zip(lst1,lst2))

	n = len(pairs)
	if n == 0: return {'d':nastr,'se':nastr, 'df':nastr,'t':nastr,'P':nastr, 'n':n, 'tails':tails}
	
	diffs = map(lambda x: subna(x[1],x[0]), pairs)
	dmean = avg(diffs)
	dse = stderr(diffs)
	t = dmean / dse
	df = n-1
	p = TdistPDF(t,df,tails)

	if verbose:
		print 'paired (2 sample) t-test for n=%s observations:'%(n)
		print 'avg d=%f' % (dmean)
		print 'SE d=%f'%(dse)
		print 'dof=%.3f' % (df)
		print 't=%s' % (t)

	return {'d':dmean, 'se':dse, 'df':df, 't':t, 'P':p, 'n':n, 'tails':tails}


# deprecated
tdistCDF = lambda t, df: TdistCDF(t,df)
tdistPDF = lambda t, df, tails: TdistPDF(t, df, tails)

def slopesTest(lr1,lr2):
	import math
	b1 = lr1['m']
	sb1 = lr1['sm']
	n1 = lr1['n']
	b2 = lr2['m']
	sb2 = lr2['sm']
	n2 = lr2['n']
	
	t = (b1 - b2) / math.sqrt(sb1**2 + sb2**2)
	df = n1 + n2 - 4
	p = ut.TdistPDF(t,df)
	return {'t':t,'df':df,'p':p}
	

def FdistCDF(F, df1, df2, dt=0.0000003): 
	"Change the value of dt to make it faster/more accurate if need be."
	import maths
	try: return maths.FdistCDF(F, df1, df2, dt)
	except TypeError: return nastr

def fisherExactTest(a,b=None,c=None,d=None,tails=2,alpha=0.05):
	"""Calculates exact p-value for a Fisher's exact test. 
	Can specify one or two-tailed test.
	"""
	import sio as io
	tmp = 'rtempfile.txt'
	
	if b == None:
		tab = a
		a = tab[0][0]
		b = tab[0][1]
		c = tab[1][0]
		d = tab[1][1]
	
	alt = 'two.sided'
	if tails == 1: alt = 'greater'
	elif tails == -1: alt = 'less'
	
	Rscpt = """
	m = cbind(0,c(2,2))
	m[1,1] = """+str(a)+"""
	m[1,2] = """+str(b)+"""
	m[2,1] = """+str(c)+"""
	m[2,2] = """+str(d)+"""
	
	ft = fisher.test(m,alternative=\""""+alt+"""\",conf.level="""+str(1-float(alpha))+""")
	write.table(ft$p.value,\""""+tmp+"""\")
	"""
	
	fh = open('tmpscript.R','w')
	print >> fh, Rscpt
	fh.close()
	
	os.system('R --vanilla < tmpscript.R > crap.out')
	
	fh = open(tmp)
	d = fh.readlines()
	
	os.system("rm -f "+tmp)
	os.system("rm -f crap.out")
	os.system("rm -f tmpscript.R")
	
	crap,pvalue = d[1].split(' ')
	
	return float(pvalue)


# generalize this to any size matrix
def batchFisherExactTestOLD(dat,tails=2,alpha=0.05):
	"""Compute fisher exact pvalues for each of a list of [[a,b],[c,d]] tables or a,b,c,d tuples"""
	
	import sio as io
	tmp = 'rtempfile.txt'
	
	alt = 'two.sided'
	if tails == 1: alt = 'greater'
	elif tails == -1: alt = 'less'
	
	Rscpt = "r = c()\n"
	
	for tup in dat:
		a,b,c,d = 0,0,0,0
		if len(tup) == 2:
			a = tup[0][0]
			b = tup[0][1]
			c = tup[1][0]
			d = tup[1][1]
		elif len(tup) == 4:
			a,b,c,d = tup
		else:
			raise TypeError, "batchFisherExactTest:table must be 2x2 or 4-tuple."
		
		Rscpt += """
		m = cbind(0,c(2,2))
		m[1,1] = """+str(a)+"""
		m[1,2] = """+str(b)+"""
		m[2,1] = """+str(c)+"""
		m[2,2] = """+str(d)+"""
		ft = fisher.test(m,alternative=\""""+alt+"""\",conf.level="""+str(1-float(alpha))+""")
		r = append(r,ft$p.value)
		"""
	Rscpt += """write.table(r,\""""+tmp+"""\")"""
	
	fh = open('tmpscript.R','w')
	print >> fh, Rscpt
	fh.close()
	
	os.system('R --vanilla < tmpscript.R > crap.out')
	
	dat,ps,crap2 = io.readMatrix(tmp)
	pvalues = []
	for p in ps:
		bad, real = p.split(' ')
		pvalues += [float(real)]
	
	os.system("rm -f "+tmp)
	os.system("rm -f crap.out")
	os.system("rm -f tmpscript.R")
	
	return pvalues


def batchFisherExactTest(dat,tails=2,alpha=0.05, outdir='/tmp/simolacode/'):
	"""Compute fisher exact pvalues for each of a list of [[a,b],[c,d]] tables or a,b,c,d tuples"""
	
	import sio as io
	
	tmp = outdir+'rtempfile.txt'
	
	io.createdir(outdir)
	
	alt = 'two.sided'
	if tails == 1: alt = 'greater'
	elif tails == -1: alt = 'less'
	
	Rscpt = "r = c()\n"
	
	count = 0
	for tup in dat:
		fn = outdir+'tmpfishermat-%s.txt'%(count)
		io.printMatrix(tup,file=fn)
		
		Rscpt += """
		m = read.delim("%s",header=FALSE)
		ft = fisher.test(m,alternative=\""""%(fn)+alt+"""\",conf.level="""+str(1-float(alpha))+""")
		r = append(r, ft$p.value)
		"""
		count += 1
	Rscpt += """write.table(r, sep="\t", row.names=FALSE, col.names=FALSE,\""""+tmp+"""\")"""
	
	
	fh = open('tmpscript.R','w')
	print >> fh, Rscpt
	fh.close()
	
	os.system('R --vanilla < tmpscript.R > crap.out')
	
	pvalues = []#; estimates = []
	pairs = []
	try:
		pairs = io.readList(tmp)
	except IOError: pass
	
	for p in pairs:
		# pval, est = p.split('\t')
		pvalues += [float(p)]
		# estimates += [float(est)]
	
	os.system("rm -f "+tmp)
	os.system("rm -f crap.out")
	os.system("rm -f tmpscript.R")
	
	count = 0
	for i in range(len(dat)):
		fn = outdir+'tmpfishermat-%s.txt'%(count)
		os.system('rm -f "%s"'%(fn))
		count += 1
	
	return pvalues#, estimates


def batchChiSquareTest(dat, Psimulate=False, outdir=''):
	"""Compute fisher exact pvalues for each of a list of [[a,b],[c,d]] tables or a,b,c,d tuples"""
	
	import sio as io
	
	tmp = outdir+'rtempfile.txt'
	
	Rscpt = "r = c()\n"
	
	count = 0
	for tup in dat:
		fn = outdir+'tmpfishermat-%s.txt'%(count)
		io.printMatrix(tup,file=fn)
		
		Rscpt += """
		m = read.delim("%s",header=FALSE)
		ft = chisq.test(m, simulate.p.value=TRUE)"""%(fn)+"""
		r = append(r, ft$p.value)
		"""
		count += 1
	Rscpt += """write.table(r, sep="\t", row.names=FALSE, col.names=FALSE,\""""+tmp+"""\")"""
	
	fh = open('tmpscript.R','w')
	print >> fh, Rscpt
	fh.close()
	
	os.system('R --vanilla < tmpscript.R > crap.out')
	
	pairs = io.readList(tmp)
	pvalues = []#; estimates = []
	for p in pairs:
		# pval, est = p.split('\t')
		pvalues += [float(p)]
		# estimates += [float(est)]
	
	os.system("rm -f "+tmp)
	os.system("rm -f crap.out")
	os.system("rm -f tmpscript.R")
	
	count = 0
	for i in range(len(dat)):
		fn = outdir+'tmpfishermat-%s.txt'%(count)
		os.system('rm -f "%s"'%(fn))
		count += 1
	
	return pvalues#, estimates



def WilksLambdaPDF(rho, n, q):
	"""Return variance ratio (F-statistic)"""
	
	L = (1.0 - rho**2)
	F = ( (1.0 - L)*(n - float(q)) )/(L*float(q))
	import maths
	pval = maths.FdistCDF(F,q,n-q)
	# print 'wilks', 'L', L, 'F', F, 'pval', pval
	return 1 - pval # testing for large correlations

def ccaPDF(rho, n, p, q, i=1):
	"""
	Return Chisquare probability of observing the vector of
	canonical correlations, as a likelihood ratio test of 
	all canonical variables.
	"""
	asum = sum([ln(1.0 - rho[j]**2) for j in range(i-1,len(rho))])
	# asum = sum([ln(1.0 - rho[j]**2) for j in range(i-1,p)])
	# asum = sum([(1-rho[j]**2) for j in range(i,p)])
	x2 = -(n-.5*(p+q+1))*asum
	dof = (p-1)*(q-1)
	# dof = (m-i+1)*(n-i+1)
	# print x2, dof, 'p', p, 'q', q, 'n', n, 'asum', asum
	print 'calling chisquarepdf(%.4f, %d)' % (x2, dof), 'n',n,'.5x', .5*(p+q+1), 'asum', asum
	return ChiSquarePDF(x2, dof)


def ChiSquareTest(tab, tails=2, Yates=True, pseudoCount=0, minP=1e-18):
	"""Apply a Chi square test of association to an RxC contingency table."""
	O = tab; nr = len(O); nc = len(O[0])
	
	if pseudoCount > 0:
		O = [[O[i][j]+pseudoCount for j in range(nc)] for i in range(nr)]
	
	R = [sum(O[i]) for i in range(nr)]
	C = [sum(vslice(O,j)) for j in range(nc)]
	n = sum(R)
	E = [[nastr for j in range(nc)] for i in range(nr)]
	
	X2 = 0.
	if n == 0: return {'X2':X2, 'P':nastr, 'df':nastr}
	df = (nr-1)*(nc-1)
	try:
		for i in range(nr):
			for j in range(nc):
				E[i][j] = R[i]*C[j] / float(n)
				# print 'E(%f) = R(%f)xC(%f)/%d'%(E[i][j],R[i],C[j],n)
				if Yates==True: X2 += (abs(O[i][j]-E[i][j])-.5)**2 / E[i][j]
				else: X2 += (O[i][j]-E[i][j])**2 / E[i][j]
				# print 'O%d%d=%f, E%d%d=%f, X2=%f' % (i,j,O[i][j],i,j,E[i][j],X2)
	except ZeroDivisionError: 
		# print >> sys.stderr, 'Utilities.py::ChiSquareTest: ZERO DIVISION ERROR in E table', E
		return {'X2':0., 'P':nastr, 'df':df, 'O':O, 'E':E}
	
	P = ChiSquarePDF(X2,df,tails=tails)
	if P != nastr and P == 0: P = minP
	# for i in range(nr):
		# for j in range(nc):
			# E[i][j] = '%.2f'%(E[i][j])
	
	return {'X2':X2, 'O':O, 'E':E, 'P':P, 'df':df}


chiSquareTest = ChiSquareTest
def ChiSquareEnrichment(set1,set2,all, tails=2, Yates=True, pseudoCount=0):
	# 
	a = len(setIntersection(set1,set2)) # brain + HT
	b = len(setDifference(set2, set1)) # not brain + HT
	c = len(setDifference(set1, set2)) # brain not HT
	d = len(setDifference(setDifference(all,set1), set2)) # not brain not de
	tab = [[a,b],[c,d]]
	res = ChiSquareTest(tab, tails=tails, Yates=Yates, pseudoCount=pseudoCount)
	return res
	

def ChiSquareFreqTest(obs, exp, tails=2, pseudo=0):
	assert len(obs)==len(exp),\
	  'input vectors must match in length: |obs|=%s, |exp|=%s'%(len(obs), len(exp))
	
	if pseudo > 0:
		obs = map(addna1(pseudo), obs)
		exp = map(addna1(pseudo), exp)
	
	D = map(lambda x,y: (x-y)*(x-y)/y, obs, map(float,exp))
	S = sum(D)
	df = len(obs) - 1
	
	P = ChiSquarePDF(S, df, tails=tails)

	O = map(lambda x: x/float(sum(obs)), obs)
	E = map(lambda x: x/float(sum(exp)), exp)
	D = map(lambda x,y: '%.4f'%(x-y), O,E)

	return {'X2':S, 'P':P, 'E':E, 'O':O, 'D':D, 'df':df}


def ChiSquarePDF(x2, dof=1, tails=2):
	# import maths
	# pval = 1-maths.ChiSquareCDF(x2, dof)
	from scipy.stats import chi2
	pval = 1-chi2.cdf(x2,dof)
	
	if tails == 1: return pval/2.
	return pval
	
	# big and small can be significant
	# if tails == 2 and pval > .5: pval = 1.0-pval
	# # big is significant
	# elif tails == 1: pval = 1.0 - pval
	# # small is significant
	# elif tails == -1: pval = pval
	# return pval

def ChiSquareCDF(x2, dof=1):
	# import maths
	# return maths.ChiSquareCDF(x2, dof)
	from scipy.stats import chi2
	return chi2.cdf(x2,dof)

def optimizeX2(df, alpha, delta=1e-5):
	"""Return the X2 value for a given alpha and dof"""
	
	import maths
	# figure out the X2 value - inverse problem
	X2 = 5.0; X2prev = 10.0
	pval = 1.0
	while abs(pval - alpha) > delta:
		error = pval-alpha #X2 - X2prev
		X2prev = X2
		if error > 0: X2 -= abs(error)/2.0
		else: X2 += abs(error)/2.0
		pval = maths.ChiSquareCDF(X2, df)
		# print 'trial', error, X2, pval
	return X2

def optimizeZ(alpha, delta=1e-5, init=0.):
	import maths
	# figure out the X2 value - inverse problem
	Z = init; Zprev = 10.0
	if Z == Zprev: Zprev = Z + 2.0
	pval = 1.0
	while abs(pval - alpha) > delta:
		error = pval-alpha
		Zprev = Z
		if error > 0: Z -= abs(error)/2.0
		else: Z += abs(error)/2.0
		pval = maths.GaussianCDF(Z, 0., 1.)
		# print 'trial', error, Z, pval
	return Z


def ChiSquareLimits(n, alpha=0.05, method='equal tails', delta=1e-5):
	"""Return a pair (L,H) which specify the lower and upper equal tails confidence limits for the variance."""
	
	if n < 3: return (nastr, nastr)
	X2L = optimizeX2(n-1, alpha/2.0)
	X2H = optimizeX2(n-1, 1-(alpha/2.0))
	return (X2L,X2H)


def meanWithConfidence(x, na=True, ret=nastr, alpha=.05):
	# this is for confidence around statistic
	mean = avg(x)
	# confidence: 1-alpha around mean
	n = samplesize(x)
	
	# find the Z such that 1-alpha lies between -z and +z
	zval = abs(optimizeZ(alpha/2.))
	# print 'at alpha', alpha, 'zval is ', zval
	sem = divna(sd(x),math.sqrt(n))
	L1 = nastr; L2 = nastr
	try:
		L1 = mean - zval*sem; L2 = mean + zval*sem
	except TypeError: pass
	# print 'L1/L2', L1, L2
	return mean, L1, L2
	


def varianceWithConfidence(x, na=True, ret=nastr, alpha=0.05, method='equal tails', delta=1e-5):
	"""Compute sample variance of input data and confidence limits on the estimate."""
	xf = []
	if na: xf = filterstr(nastr)(x)
	else: xf = x[:]
	
	n = len(xf)
	var = variance(xf, na=False, ret=ret)
	low, high = ChiSquareLimits(n, alpha=alpha, method=method, delta=delta)
	L1 = nastr; L2 = nastr
	try:
		L1 = var*(n-1) / high
		L2 = var*(n-1) / low
	except TypeError: pass
	
	return var, L1, L2


def varianceTest(Y, stat=median):
	"Input a list of lists representing sample values for two or more distributions to be compared. median stat is the Brown-Forsythe test, while mean is the Levene test."
	
	if len(Y) == 0:
		print >> sys.stderr, 'varianceTest() error: input list must contain at least 2 lists.'
	
	n = [len(j) for j in Y]
	bigN = sum(n)
	p = len(Y)
	df1 = (p-1)
	df2 = (bigN-p)
	
	# abs deviation from mean/median
	Ytilde = [stat(j) for j in Y]
	Z = [[abs(Y[j][i] - Ytilde[j]) for i in range(n[j])] for j in range(p)]
	
	Zbar = [avg(Z[j]) for j in range(p)]
	Zbarbar = avg(flatten(Z))
	
	B = sum([n[j] * (Zbar[j] - Zbarbar)**2 for j in range(p)])
	W = sum([sum([(Z[j][i] - Zbar[j])**2 for i in range(n[j])]) for j in range(p)])
	
	F = (df2/df1) * (B/W)
	P = FdistCDF(F,df1,df2, dt=0.0000000001)
	
	return {'F':F, 'df1':df1, 'df2':df2, 'P':P}


KScrit = lambda alpha, n1, n2: math.sqrt(-.5*math.log(alpha/2., math.e)) * math.sqrt((n1+n2)/float(n1*n2))# * float(n1*n2)

def KStest(lst1,lst2, alpha=.05):
	"""Input 2 lists representing cumulative count (not frequency) distributions. Performs a 2-tailed KS-test."""
	
	n1 = sum(lst1)
	n2 = sum(lst2)
	
	# convert to CDF
	f1 = map(lambda x: x/float(n1), lst1)
	f2 = map(lambda x: x/float(n2), lst2)
	
	# compute abs deviation
	ds = map(lambda x,y: abs(x-y), f1, f2)
	idx = argmax(ds)
	D = ds[idx]
	# D = avg(ds)
	
	Ka = math.sqrt(-.5*math.log(alpha/2., math.e))
	Da = Ka * math.sqrt((n1+n2)/float(n1*n2))
	P = min(1.0,2.0*math.e**( -2.* D**2. * (n1*n2 / float(n1+n2)) ))
	
	return {'D':D, 'P':P, 'idx':idx, 'n1':n1, 'n2':n2, 'Dcrit':Da, 'alpha':alpha}


def getAlleles(seq, snpref, LAPLACE=0):
	# get the allele frequencies from a samtools pileup
	# don't include N since it is not informative
	AF = {'A':LAPLACE, 'T':LAPLACE, 'C':LAPLACE, 'G':LAPLACE}
	
	# populate read allele freqs
	# parse sequence into vector of alleles, using samtools patterns
	i = 0
	seq = seq.upper()
	while i < len(seq):
		if re.match('\,|\.', seq[i]): AF[snpref] += 1
		elif re.match('[atcgATCG]', seq[i]): AF[seq[i]] += 1
		i+=1
	
	return [AF['A'],AF['T'],AF['C'],AF['G']]




def alleleTest(obs, exp, snpref='', SEQ=True, method='equalize', direction='', LAPLACE=0, VERBOSE=False):
	
	AFo = obs; AFe = exp
	if snpref != '' or (SEQ==True and snpref != ''):
		AFo = getAlleles(obs, snpref, LAPLACE)
		AFe = getAlleles(exp, snpref, LAPLACE)
	elif LAPLACE>0:
		AFo = map(lambda x: x+LAPLACE, AFo)
		AFe = map(lambda x: x+LAPLACE, AFe)
	
	# equalize the counts using the average read number
	if method == 'equalize':
		to = float(sum(AFo))
		te = float(sum(AFe))
		nr = (sum(AFo)+sum(AFe))/2.
		AFo = map(lambda x: x*(nr/to), AFo)
		AFe = map(lambda x: x*(nr/te), AFe)
	elif method == 'normalize':
		# normalize the expected to the observed counts
		to = float(sum(AFo))
		te = float(sum(AFe))
		AFe = map(lambda x: x*(to/te), AFe)
	
	res = ChiSquareFreqTest(AFo,AFe)
	
	if direction != '':
		res2 = ChiSquareFreqTest(AFe,AFo)
		# return the better score, where direction is max or min function
		ret = direction([(res['X2'],res), (res2['X2'],res2)])
		res = ret[1]
		
	if VERBOSE:
		print 'X2=%s, df=%s, P=%s'%(res['X2'], res['df'], res['P'])
		print ' O:', map(str,map(r2, AFo))
		print ' E:', map(str, map(r2, AFe))
		print
	
	return {'P':res['P'], 'AFo':AFo, 'AFe':AFe, 'X2':res['X2']}


frequencyTest = alleleTest

def MHtest(tab, cc=.5):
	"""Mantel-Haenszel test for replicated 2x2 tables
	
	Data must be organized as follows:
	R1	a	b	c	d
	R2	a	b	c	d
	...
	Rn	a	b	c	d
	
	where Ri are the replicated tables
	"""
	
	res = {'omega':nastr, 'omegas':[], 'X2':nastr, 'df':nastr, 'P':nastr, 'N':[]}
	
	tab = castemat(tab,floatna)
	
	N = map(sum,tab) # table sums
	res['N'] = N
	ad = map(lambda x: multna(addna(x[0],cc),addna(x[3],cc)), tab)
	bc = map(lambda x: multna(addna(x[1],cc),addna(x[2],cc)), tab)
	# print 'ad/n', map(lambda x,y: divna(x,y), ad,N)
	# print 'bc/n', map(lambda x,y: divna(x,y), bc,N)
	numers = map(lambda x,y: divna(x,y), ad,N)
	denoms = map(lambda x,y: divna(x,y), bc,N)
	
	res['omega'] = divna(sumna(numers),sumna(denoms))
	res['omegas'] = map(lambda x,y: divna(x,y), numers,denoms)
	
	# X2 statistic
	suma = sumna(map(lambda x: x[0], tab))
	sumabac = sumna(map(lambda x: divna(multna(addna(x[0],x[1]), addna(x[0],x[2])), sumna(x)), tab))
	x2numer = powna(subna(absna(subna(suma,sumabac)), .5))(2)
	x2denom = sumna(map(lambda x: divna(prodna([addna(x[0],x[1]), addna(x[0],x[2]), addna(x[1],x[3]), addna(x[2],x[3])]), subna(sumna(x)**3, sumna(x)**2)), tab))
	res['X2'] = divna(x2numer,x2denom)
	res['df'] = 1
	if res['X2'] != nastr:
		res['P'] = ChiSquarePDF(res['X2'],res['df'])
	return res



def mantelCorrelation(X,Y,permutations=1, outdir='', mantelfile='/Users/simola/bin/mantel.test.R', rmtmp=True):
	import sio as io
	fh = open(mantelfile)
	scpt = ''.join(fh.readlines())
	fh.close()
	
	io.printTable(X,file=outdir+'matrixX.txt')
	io.printTable(Y,file=outdir+'matrixY.txt')
	
	# print 'table printed'
	
	scpt += """
	X = read.delim("%s")
	Y = read.delim("%s")
	#X = X[,2:length(X)]
	#Y = Y[,2:length(Y)]
	mant = mantel.fcn(as.matrix(X), as.matrix(Y), nperm=%d)
	write(c(mant$R,mant$UpperP), file="%s")
	""" % (outdir+'matrixX.txt', outdir+'matrixY.txt', permutations, outdir+'rcorout.txt')
	
	fh = open(outdir+'Rscript.R', 'w')
	print >> fh, scpt
	fh.close()
	os.system('R --vanilla < %s > %s'%(outdir+'Rscript.R', outdir+'Rscript.out.txt'))
	
	d,r,c = io.readMatrix(outdir+'rcorout.txt', header=False)
	
	if rmtmp:
		os.system('rm -f "%s"'%(outdir+'matrixX.txt'))
		os.system('rm -f "%s"'%(outdir+'matrixY.txt'))
	
	vals = map(float, r[0].split(' '))
	R = vals[0]
	P = vals[1]
	return {'R':R, 'P':P, 'permutations':permutations}
	


def runRscript(script, outdir=''):
	fh = open(outdir+'runRscript.R','w')
	print >> fh, script
	fh.close()
	# print 'finished'
	os.system('R --vanilla < "'+outdir+'runRscript.R'+'" > "'+outdir+'runRscript.out.R'+'"')
	# os.system('rm '+outdir+'runRscript.R'); os.system('rm '+outdir+'runRscript.out.R')
	return 1


def regressionR(dat, intercept=True, outdir='', method='pearson', weighted=False, na=True, logtransform=''):
	"""Run a linear regression on pairs of x,y data in R."""
	import sio as io
	
	datr = dat[:]
	xvals = [dat[i][0] for i in range(len(dat))]
	yvals = [dat[i][1] for i in range(len(dat))]
	
	if logtransform == 'x':    xvals = map(log, xvals)
	elif logtransform == 'y':  yvals = map(log, yvals)
	elif logtransform == 'xy':
		xvals = map(log, [dat[i][0] for i in range(len(dat))])
		yvals = map(log, [dat[i][1] for i in range(len(dat))])
	
	if method == 'spearman':
		# rank data
		xvals = rank(xvals,ties='split')
		yvals = rank(yvals,ties='split')
	
	datr = zip(xvals,yvals)
	if na: 
		datr = pairedfilterstr(nastr)(datr)
		xvals = [datr[i][0] for i in range(len(datr))]
		yvals = [datr[i][1] for i in range(len(datr))]
		
	io.printTable(datr, ['X', 'Y'], file=outdir+'regression.dat')
	
	scpt = """setwd("%s")\ndat=read.delim("regression.dat", header=TRUE)\n"""%(outdir)
	
	if intercept: 
		scpt += "model = lm(dat$Y ~ dat$X)\nsummary(model)$coef\n"
		if weighted:
			scpt += "wts <- 1/fitted( lm(abs(residuals(model)) ~ dat$X) )^2\n"
			scpt += "model = lm(dat$Y ~ dat$X, weights=wts)\nsummary(model)$coef\n"
		scpt += """int = summary(model)$coef[1,1]
pvalint = summary(model)$coef[1,4]
slope = summary(model)$coef[2,1]
pvalslope = summary(model)$coef[2,4]"""
	else: 
		scpt += "model = lm(dat$Y ~ 0 + dat$X)\nsummary(model)$coef\n"
		if weighted:
			scpt += "wts <- 1/fitted( lm(abs(residuals(model)) ~ 0 + dat$X) )^2\n"
			scpt += "model = lm(dat$Y ~ 0 + dat$X, weights=wts)\nsummary(model)$coef\n"
		
		scpt += """int = 0
pvalint = 1
slope = summary(model)$coef[1,2]
pvalslope = summary(model)$coef[1,4]"""
	
	scpt += """
fstats = summary(model)$fstatistic
r2 = summary(model)$r.squared
ar2 = summary(model)$adj.r.squared
correlation = cor(dat$X,dat$Y,method="pearson")
spearman = cor(dat$X,dat$Y,method="spearman")
sse = sum(model$residuals^2)
 
tab = c(int,slope,pvalint,pvalslope,r2,ar2,fstats,correlation,spearman,sse)
write.table(tab,file=\"regressionresults.txt\",sep=\"\\t\")\n
	"""
	
	runRscript(scpt, outdir=outdir)
	
	d,r,c = io.readMatrix(outdir+'regressionresults.txt')
	
	results = {}
	results['intercept'] = d[0][0]
	results['slope'] = d[1][0]
	results['pvalintercept'] = d[2][0]
	results['pvalslope'] = d[3][0]
	
	results['Fstatistic'] = d[6][0]
	results['df top'] = d[7][0]
	results['df bottom'] = d[8][0]
	
	results['r'] = d[9][0]
	results['s'] = d[10][0]
	
	results['R2'] = d[4][0]
	results['adjR2'] = d[5][0]
	
	results['SSE'] = d[11][0]
	
	for k in results.keys(): results[k] = float(results[k])
	
	# compute points along regression line using original x values
	fitted = []
	
	# xvals = [datr[i][0] for i in range(len(datr))]
	
	# print >> sys.stderr, 'TEST', results['intercept'], results['slope'] ,xvals

	func = lambda x: results['intercept'] + results['slope'] * float(x)
	
	yvals = map(func, xvals)
	results['x'] = xvals
	results['y'] = yvals
	results['fitted'] = yvals
	
	return results


# fitted_values = lambda b,m: lambda x: b + m*x
#
def linearRegression(X, Y=[], order=1, na=True, stats=1):
	
	# currently only performs a first order regression
	
	# Basic computations to save a little time.
	
	if Y: X = zip(X,Y)
	if na: X = pairedfilterstr(nastr)(X)
	
	Y = map(float,[r[1] for r in X])
	X = map(float,[r[0] for r in X])
	
	length = len(X)
	sum_x = sum(X)
	sum_y = sum(Y)
	
	# Sx^2, and Sxy respectively.
	sum_x_squared = sum(map(lambda a: a * a, X))
	covariance = sum([X[i] * Y[i] for i in range(length)])
	
	# Magic formulae!  
	m = (covariance - (sum_x * sum_y) / length) / (sum_x_squared - ((sum_x ** 2) / length))
	b = (sum_y - m * sum_x) / length
	
	res = {'m':m,'b':b,'p':nastr}
	# significance testing
	if stats:
		import maths
		fit = lambda b,m: lambda x: b + m*x
		
		r = correlation(X,Y)
		syx = 1e-10
		sx = 1e-10
		try:
			syx = stdev(Y) * math.sqrt( (1-r**2)*((length-1)/(length-2)) )
			sx = stdev(X)
		except ValueError: pass # log of 0
		sm = syx / (sx * math.sqrt(length-1))
		res['r'] = r
		res['syx'] = syx
		res['sx'] = sx
		res['sm'] = sm
		res['n'] = length
		
		ybar = avg(Y)
		
		res['df top'] = 1
		res['df bottom'] = res['n'] - (order + 1)
		
		res['fit'] = map(fit(res['b'],res['m']), X)
		res['res'] = map(lambda x,z: x-z, Y, res['fit'])
		
		res['sst'] = sum(map(lambda y: (y - ybar)*(y - ybar), Y))
		res['ssr'] = sum(map(lambda y: (y - ybar)*(y - ybar), res['fit']))
		res['sse'] = sum(map(lambda x: x*x, res['res']))
		res['r2'] = res['ssr']/res['sst']
		
		if res['sse'] == 0: res['sse'] = 1/GLOB_INF
		res['F'] = (res['ssr']/res['df top']) / (res['sse']/res['df bottom'])
		
		# P = 
		# if F >= 1: P = 1-P
		res['p'] = 1. - FdistCDF(res['F'],res['df top'],res['df bottom'])
		
	return res

#
def slopesTest(a,b):
	import math
	
	lr1 = linearRegression(a,stats=1)
	lr2 = linearRegression(b,stats=1)
	
	b1 = lr1['m']
	sb1 = lr1['sm']
	n1 = lr1['n']
	b2 = lr2['m']
	sb2 = lr2['sm']
	n2 = lr2['n']

	t = (b1 - b2) / math.sqrt(sb1**2 + sb2**2)
	df = n1 + n2 - 4
	p = TdistPDF(t,df)
	return {'t':t,'df':df,'p':p, 'm':[lr1['m'],lr2['m']], 'b':[lr1['b'],lr2['b']]}


def slopesTestM(lr1,lr2):
	import math
	b1 = lr1['m']
	sb1 = lr1['sm']
	n1 = lr1['n']
	b2 = lr2['m']
	sb2 = lr2['sm']
	n2 = lr2['n']

	t = (b1 - b2) / math.sqrt(sb1**2 + sb2**2)
	df = n1 + n2 - 4
	p = TdistPDF(t,df)
	return {'t':t,'df':df,'p':p}

	

def regressionOLD(X, Y, order=1, type='full', na=True, nastr=nastr):
	"""Uses the regression.py code."""
	
	import regression as reg, maths
	
	fit = lambda b,m: lambda x: b + m*x
	
	dat = [[Y[i],X[i]] for i in range(len(X))]
	if na: dat = pairedfilterstr(nastr)(dat)
	
	ret = {}
	ret['coefs'] = reg.linearRegression(dat, order)
	
	ret['X'] = [dat[i][1] for i in range(len(dat))]
	ret['Y'] = [dat[i][0] for i in range(len(dat))]
	
	ybar = avg(ret['Y'])
	
	ret['n'] = len(ret['X'])
	ret['df top'] = 1
	ret['df bottom'] = ret['n'] - (order + 1)
	
	ret['fitted values'] = map(fit(ret['coefs'][0], ret['coefs'][1]), ret['X'])
	ret['residuals'] = map(lambda x,z: x-z, ret['Y'], ret['fitted values'])
	
	ret['sst'] = sum(map(lambda y: (y - ybar)*(y - ybar), ret['Y']))
	ret['ssr'] = sum(map(lambda y: (y - ybar)*(y - ybar), ret['fitted values']))
	ret['sse'] = sum(map(lambda x: x*x, ret['residuals']))
	
	ret['R'] = correlation(ret['Y'],ret['X'])
	ret['R2'] = ret['ssr']/ret['sst']
	
	if type == 'full':
		if ret['sse'] == 0: ret['sse'] = 1/GLOB_INF
		ret['F'] = (ret['ssr']/ret['df top']) / (ret['sse']/ret['df bottom'])
		ret['P'] = 1.-maths.FdistCDF(ret['F'], ret['df top'], ret['df bottom'])
	
	return ret



#
def regression(X, Y=None, order=1, type='full', na=True, nastr=nastr):
	"""Uses the regression.py code."""
	
	import regression as reg, maths
	
	fit = lambda b,m: lambda x: b + m*x
	
	if Y:
		dat = [[Y[i],X[i]] for i in range(len(X))]
	else:
		dat = [[X[i][0],X[i][1]] for i in range(len(X))]
	
	if na: dat = pairedfilterstr(nastr)(dat)
	
	ret = {}
	ret['coefs'] = reg.linearRegression(dat, order)
	
	ret['X'] = [dat[i][1] for i in range(len(dat))]
	ret['Y'] = [dat[i][0] for i in range(len(dat))]
	
	ybar = avg(ret['Y'])
	
	ret['n'] = len(ret['X'])
	ret['df top'] = 1
	ret['df bottom'] = ret['n'] - (order + 1)
	
	ret['fitted values'] = map(fit(ret['coefs'][0], ret['coefs'][1]), ret['X'])
	ret['residuals'] = map(lambda x,z: x-z, ret['Y'], ret['fitted values'])
	
	ret['sst'] = sum(map(lambda y: (y - ybar)*(y - ybar), ret['Y']))
	ret['ssr'] = sum(map(lambda y: (y - ybar)*(y - ybar), ret['fitted values']))
	ret['sse'] = sum(map(lambda x: x*x, ret['residuals']))
	
	ret['R'] = correlation(ret['Y'],ret['X'])
	ret['R2'] = ret['ssr']/ret['sst']
	
	if type == 'full':
		if ret['sse'] == 0: ret['sse'] = 1/GLOB_INF
		ret['F'] = (ret['ssr']/ret['df top']) / (ret['sse']/ret['df bottom'])
		ret['P'] = 1.-maths.FdistCDF(ret['F'], ret['df top'], ret['df bottom'])
	
	return ret


def anova(X):
	# the first dimension of mat defines groups
	# compute an anova
	
	a = len(X) # number of groups
	n = [len(X[i]) for i in range(len(X))] # data points per group
	bign = sumna(n) # total data points
	
	# means
	xbar = [avg(X[i]) for i in range(a)]
	bigxbar = divna(sumna(flatten(X)),bign)
	
	B = sum([n[i]*(xbar[i] - bigxbar)*(xbar[i] - bigxbar) for i in range(a)])
	W = [0 for i in range(a)]
	for i in range(a):
		for j in range(n[i]):
			W[i] += (X[i][j] - xbar[i])*(X[i][j] - xbar[i])
	
	F = (B/float(a-1))/(sum(W)/float(bign-a))
	
	P = FdistCDF(F,a-1,bign-a)
	if F >= 1: P = 1-P
	
	return {'F':F, 'B':B, 'W':sum(W), 'within':W, 'df2':bign-a, 'df1':a-1, 'p':P}



# PWM position weight matrix
def PFMtoPWM(mat, prior=[0.25, 0.25, 0.25, 0.25]):
	p = prior
	lg = lambda x: math.log(x,2)
	N = [sum(pos) for pos in mat]
	rN = [math.sqrt(sum(pos)) for pos in mat]
	
	w = map(lambda j: map(lambda i: lg( ( mat[j][i] + rN[j] * p[i] ) / ( N[j] + rN[j] ) / p[i] ), xrange(4)), range(len(mat)))
	
	return w





# various functions for heterochrony analysis
# -------------------------------------------
modulate = lambda period: lambda X: map(lambda xi: xi%period, X)

def makePeriodic(Y, perIdx):
	"""Modulate Y at perIdx and average remaining data with Y[0:perIdx+len(Y)]"""
	Z = [Y[i] for i in range(perIdx)]
	for i in range(len(Y)-perIdx): Z[i] = (Z[i]+Y[perIdx+i])/2.
	return Z

skewFunc = []
skewFunc += [lambda Z: lambda t: t] 
skewFunc += [lambda Z: lambda t: Z[0]*t] 
skewFunc += [lambda Z: lambda t: Z[0]*t + Z[1]] 
skewFunc += [lambda Z: lambda t: Z[0]*t + Z[1] + Z[2]*t**2]

# also define a model with phase and power parameters
powFunc = lambda theta: lambda x: x**theta[0] + theta[1]
powFunc2 = lambda theta: lambda x: theta[0] + x**theta[1]

# import scipy.stats as ss
# betaFunc = lambda theta: lambda x: ss.beta.cdf(x,theta[1],theta[2])

sse = lambda x,y: sum(map(lambda a,b: (a-b)**2, x,y))

lsq = lambda x,y: sum(map(lambda a,b: (a-b)**2, x,y))
wlsq = lambda x,y: sum(map(lambda a,b: (1/float(abs(a-b)))*(a-b)**2, x,y))

def bruteForce(func, limits, grid=True, args=()):
	"""Brute force search method for 3 parameters"""
	opt = GLOB_INF
	optParms = None
	
	range0 = limits[0]
	range1 = limits[1]
	range2 = limits[2]
	if not grid:
		range0 = makeRange(*limits[0])
		range1 = makeRange(*limits[1])
		range2 = makeRange(*limits[2])
		
	# print 'within n=%s'%(len(range0)*len(range1)*len(range2))
	# print range0
	
	for d0 in range0:
		for d1 in range1:
			for d2 in range2:
				parms = [d0,d1,d2]
				val = func(parms, *args)
				if val < opt:
					opt = val
					optParms = parms
	return (optParms, opt)



def chronify(X, Y, parms, knots=[]):
	from scipy.interpolate import splev, splrep
	splEq = None
	if knots: splEq = splrep(X, Y, k=3, s=0, per=1, t=knots)
	else:     splEq = splrep(X, Y, k=3, s=0, per=1, task=0)
	period = X[len(X)-1]
	order = len(parms)
	newX = modulate(period)(map(skewFunc[order](parms), X))
	return splev(newX, splEq)

# exponential model
def chronifyExp(X, Y, parms, knots=[], split=0):
	from scipy.interpolate import splrep, splev
	splEq = None
	if knots: splEq = splrep(X, Y, k=3, s=0, per=1, t=knots)
	else:     splEq = splrep(X, Y, k=3, s=0, per=1, task=0)
	period = float(X[len(X)-1])
	
	# parms[1] /= period # fix parms so that alpha is scaled
	newX = map(lambda x: x/period, X)    # scale down
	newX = map(powFunc(parms), newX)     # evaluate
	newX = map(lambda x: x*period, newX) # scale up
	newX = modulate(period)(newX)        # modulate
	return splev(newX, splEq)


def chronifyExp2point(X, Y, parms, knots=[], split=0):
	from scipy.interpolate import splrep, splev
	splEq = None
	if knots: splEq = splrep(X, Y, k=3, s=0, per=1, t=knots)
	else:     splEq = splrep(X, Y, k=3, s=0, per=1, task=0)
	period = float(X[len(X)-1])
	
	newX = map(lambda x: x/period, X)    # scale down
	
	# -----------------------
	theta1,theta2 = (parms[0],parms[1]), (parms[0],parms[2])
	split = int(split*len(X))
	
	if split > 0 and split < len(X):
		X1 = map(powFunc2(theta1), X[0:split])
		X2 = map(powFunc2(theta2), X[split-1:len(X)])
	
		# match up the knot point
		dev = float(min(X2)-max(X1))
		X2 = map(lambda x: x-dev, X2)
	
		# scale everything to 0,1
		newX = X1+X2[1:len(X2)]
		maxx = newX[len(newX)-1]
		newX = map(lambda x: x/maxx, newX)
	elif split == 0:
		X1 = map(powFunc2(theta2), X)
		maxx = X1[len(X1)-1]
		newX = map(lambda x: x/maxx, X1)
	elif split == len(X):
		X1 = map(powFunc2(theta1), X)
		maxx = X1[len(X1)-1]
		newX = map(lambda x: x/maxx, X1)
		
	# -----------------------
	
	newX = map(lambda x: x*period, newX) # scale up
	newX = modulate(period)(newX)        # modulate
	
	# make periodic again
	newY = splev(newX, splEq).tolist()
	splEq = None
	if knots: splEq = splrep(X, newY, k=3, s=0, per=1, t=knots)
	else:     splEq = splrep(X, newY, k=3, s=0, per=1, task=0)
	return splev(X, splEq)



# ===============================
def chronifyBeta(X, Y, parms, knots=[], resolution=None, Xnew=[]):
	from scipy.interpolate import splrep, splev
	import scipy.stats as ss
	splEq = None
	if knots: splEq = splrep(X, Y, k=3, s=0, per=1, t=knots)
	else:     splEq = splrep(X, Y, k=3, s=0, per=1, task=0)
	
	if Xnew: X = Xnew
	elif not Xnew and resolution: X = makeRange(min(X),max(X),n=resolution)
	# else use existing X
	
	period = float(X[len(X)-1])
	# newX = map(lambda x: x/period, X)      # scale down
	# newX = map(lambda x: ss.beta.cdf(x,parms[1],parms[2]), newX) # evaluate
	# newX = map(lambda x: x*period, newX)   # scale up
	# newX = map(lambda x: x+parms[0], newX) # translate
	# newX = modulate(period)(newX)          # modulate
	newX = modulate(period)(map(lambda x: parms[0] + period*ss.beta.cdf(x/period,parms[1],parms[2]), X)) # evaluate
	
	return splev(newX, splEq)

def periodify(X, Y, perIdx, parms=(0,1,1), knots=[], resolution=None, Xnew=[], mapFunc=chronifyBeta):
	"""Given an expression trajectory, make it periodic and cover exactly 
	1 period at new resolution, with possible time shift.
	"""
	Y0 = makePeriodic(kernelsmooth(Y,w=3), perIdx)
	X0 = X[0:perIdx]
	return mapFunc(X0, Y0, parms, knots=knots, resolution=resolution, Xnew=Xnew)

def subSample(Y, X): return [Y[i] for i in X]


def heterochrony_optimization_control(X, Y1, Y2, knots, period, perIdx, resolution=100, limits=(), errorFunc=sse, mapFunc=chronifyBeta, method='full', x0=(0,1,1), ftol=1e-5, scale=False):
	"""Same as het_opt but without the estimation - just a control"""
	
	import scipy.optimize as so
	
	def heterochronyObjFunc(parms, T, Y, Z, knots=[]): 
		invparms = (-parms[0], 1/parms[1], 1/parms[2])
		return errorFunc(Z, mapFunc(T, Y, parms, knots))+errorFunc(Y, mapFunc(T, Z, invparms, knots))
	
	
	def heterochronyObjFunc2(parms, gamma, T, Y, Z, knots=[]): 
		tparms = [gamma, parms[0], parms[1]]
		invparms = (-gamma, 1/parms[0], 1/parms[1])
		return errorFunc(Z, mapFunc(T, Y, tparms, knots))+errorFunc(Y, mapFunc(T, Z, invparms, knots))
	
	
	# Prepare data
	# ------------
	T = X; Y0 = Y1; Yt = Y2
	Tnew = makeRange(0,period,n=period)
	Y0 = periodify(T,Y0,perIdx,(0,1,1),knots,resolution,Xnew=Tnew)
	Yt = periodify(T,Yt,perIdx,(0,1,1),knots,resolution,Xnew=Tnew)
	
	# this is where you do estimation
	# ---------------------
	xopt = x0
	xoptinv = [-x0[0],1/float(x0[1]), 1/float(x0[2])]
	# ---------------------
	
	# now apply the transformation
	# Yopt = Yt
	# Yoptinv = Y0
	Yopt = mapFunc(Tnew, Y0, xopt, knots)
	Yoptinv = mapFunc(Tnew, Yt, xoptinv, knots)
	
	if scale:
		Yopt = standardize(Yopt); Yoptinv = standardize(Yoptinv)
		Y0 = standardize(Y0); Yt = standardize(Yt)
		
	# now compute final fit of this model
	# --------------------------------------
	# have parms - need regression in both directions
	
	res = regression(Y0, Yt, type='basic')
	err_h0 = res['sse']; r2_h0 = res['R2']
	
	res = regression(Yopt, Yt, type='basic')
	Aerr_h1 = res['sse']; Ar2_h1 = res['R2']
	An = res['n']
	res = regression(Yoptinv, Y0, type='basic')
	Berr_h1 = res['sse']; Br2_h1 = res['R2']
	Bn = res['n']
	
	print 'R2 A', Ar2_h1, 'R2 B', Br2_h1
	
	res_n = An+Bn # data points involved
	terr_h0 = err_h0; tr2_h0 = r2_h0 # i.e., 2*err_h0/2
	
	# average R2 is equivalent to the sum of SSRs / sum of sstot
	terr_h1 = (Aerr_h1 + Berr_h1)/2.
	tr2_h1 = (Ar2_h1 + Br2_h1)/2.
	
	# return the 18 pt data - output like input
	fixme = lambda Y: subSample(Y[0:len(Y)-1]+Y[0:max(X)-period+2], T)
	Y0 = fixme(Y0); Yt = fixme(Yt)
	Yopt = fixme(Yopt); Yoptinv = fixme(Yoptinv)
	
	
	return {'Y0':Y0, 'Ytarget':Yt, 'Yopt':Yopt, 'Yoptinv':Yoptinv, 'x0':x0, 'xopt':xopt, 'xoptinv':xoptinv, 'T':T, 
			'SSE_h0':err_h0, 'R2_h0':r2_h0, 'n':res_n,
			'SSE_h1':Aerr_h1, 'R2_h1':Ar2_h1, 
			'SSE_h1inv':Berr_h1, 'R2_h1inv':Br2_h1,
			
			'SSE_tot_h0':terr_h0, 'R2_tot_h0':tr2_h0,
			'SSE_tot_h1':terr_h1, 'R2_tot_h1':tr2_h1
			}

# ===============================


def splinesmooth(X,Y,smooth=.1, per=0):
	return resample(X, Y, X, splSmooth=smooth, per=per)


def resample(X, Y, newX, window=0, splSmooth=0, task=0, per=0, knots=[]):
	"""Assumes input X, Y are periodic, but not the newX"""
	import numpy as np
	from scipy.interpolate import splev, splrep
	
	pairs = [[X[i],Y[i]] for i in range(len(X))]
	pairs = pairedfilterstr(nastr)(pairs)
	T = [float(pairs[i][0]) for i in range(len(pairs))]
	V = [float(pairs[i][1]) for i in range(len(pairs))]
	V = kernelsmoothlist(V, window)
	
	tck = None
	if knots: tck = splrep(T,V, k=3, s=splSmooth, per=per, t=knots)
	else: tck = splrep(T,V, k=3, s=splSmooth, task=task, per=per)
	newY = np.transpose(splev(newX, tck))
	return newY.tolist()

def resampleIter(X, Y, newX, window=0, iterations=0, splSmooth=0, task=0, per=0, knots=[]):
	"""Assumes input X, Y are periodic, but not the newX"""
	import numpy as np
	import scipy.linalg as la
	# scipy.optimize as so
	from scipy.interpolate import splev, splrep
	pairs = [[X[i],Y[i]] for i in range(len(X))]
	pairs = pairedfilterstr(nastr)(pairs)
	T = [float(pairs[i][0]) for i in range(len(pairs))]
	V = [float(pairs[i][1]) for i in range(len(pairs))]
	
	if window: 
		for it in range(iterations): V = kernelsmoothlist(V, window)
		
	tck = None
	# print 'resample', len(T), len(V)
	
	if knots: tck = splrep(T,V, k=3, s=splSmooth, per=per, t=knots)
	else: tck = splrep(T,V, k=3, s=splSmooth, task=task, per=per)
	newY = np.transpose(splev(newX, tck))
	return newY.tolist()



def timeTransform(X, Y, period, perIdx, parms=(), knots=[], resolution=100, mapFunc=chronifyBeta, newX=[]):
	
	Tuse = modulate(period)(X)
	Y0 = makePeriodic(Y, perIdx)
	X0 = X[0:perIdx]
	
	T = None
	if newX: T = newX; Tuse = T
	else:
		resolution = (resolution == None and len(X)) or resolution
		step = period / float(resolution)
		step += step/float(resolution)
		T = [min(X) + i*step for i in range(resolution)]
	
	T0 = modulate(period)(T)
	Y0 = resample(X0, Y0, T0, window=3, per=1, knots=knots)
	if len(parms) > 0: Y0 = mapFunc(T, Y0, parms, knots, per=1)
	return Y0


def timeTransformNew(X, Y, period, perIdx, parms=(), knots=[], resolution=100, mapFunc=chronifyBeta, newX=[], periodify=False):
	"""Generate the time-domain transform of a gene expression curve given the original
	data and time parameters
	"""
	
	Tuse = modulate(period)(X)
	T = None
	
	if newX: T = newX; Tuse = T
	else:
		resolution = (resolution == None and len(X)) or resolution
		step = period / float(resolution)
		step += step/float(resolution)
		T = [min(X) + i*step for i in range(resolution)]
	
	if periodify:
		Y0 = makePeriodic(Y, perIdx)
		X0 = X[0:perIdx]
		T0 = modulate(period)(T)
		Y0 = resample(X0, Y0, T0, window=3, per=1, knots=knots)
	else:
		Y0 = resample(X, Y, modulate(period)(T), window=3, per=1, knots=knots)
	
	if len(parms) > 0: Y0 = mapFunc(T, Y0, parms, knots, per=1)
	
	return Y0
	
	# if len(newX): return Y0
	# return resample(T, Y0, Tuse, window=3, per=1, knots=knots)


def heterochrony_optimization_fminb(X, Y1, Y2, knots, period, perIdx, resolution=None, limits=(), errorFunc=sse, mapFunc=chronifyBeta, ftol=10, nstart=10, approx_grad=True):
	import scipy.optimize as so
	
	def heterochronyObjFunc(parms, T, Y, Z, knots=[]): 
		return errorFunc(Z, mapFunc(T, Y, parms, knots))
	
	x0 = (0,1,1) # beta, linear mapping, no phase
	lenexoh = len(x0)
	if not limits: limits = [(None,None) for i in range(len(x0))]
	
	# Prepare data
	# ------------
	Tpre = X[:]
	Torig = Tpre[0:perIdx] # correct time
	Tuse = modulate(period)(Tpre)
	
	# make periodic, reeval to 18 points, and standardize
	Y0 = makePeriodic(Y1, perIdx)
	Yt = makePeriodic(Y2, perIdx)
	
	# re-sample at original 18 points
	resolution = (resolution == None and len(X)) or resolution
	step = period / float(resolution)
	step += step/float(resolution)
	T = [min(Torig) + i*step for i in range(resolution)]
	
	Y0 = resample(Torig, Y0, T, window=3, per=1, knots=knots)
	Yt = resample(Torig, Yt, T, window=3, per=1, knots=knots)
	
	# standardize so that optimization fits the pattern
	Y0 = standardize(Y0)
	Yt = standardize(Yt)
	# ------------
	
	args = (T, Y0, Yt, knots)
	
	# constrained
	import random
	xopt = x0; fopt = GLOB_INF; resdct = {}
	for xx in range(nstart):
		x0t = (xx > 0 and [0]+[random.randint(0,5) for xy in range(1,lenexoh)]) or x0
		results = so.fmin_l_bfgs_b(heterochronyObjFunc, x0t, args=args, factr=ftol, bounds=limits, 
									approx_grad=True, iprint=-1, maxfun=15000, pgtol=1e-5, epsilon=1e-8, m=10)
		(Txopt, Tfopt, Tresdct) = results
		
		if xx == 0: (xopt,fopt,resdct) = (Txopt,Tfopt,Tresdct)
		elif xx > 0 and Tfopt < fopt: (xopt,fopt,resdct) = (Txopt,Tfopt,Tresdct)
	# print 'done', xopt, fopt
	# transform original expression data with null and opt parameters
	Y0 = mapFunc(T, Y0, x0, knots)
	Yopt = mapFunc(T, Y0, xopt, knots)
	Yt = mapFunc(T, Yt, x0, knots)
	
	error = errorFunc(Yopt,Yt); error0 = errorFunc(Y0, Yt)
	return {'Y0':Y0, 'Yopt':Yopt, 'Ytarget':Yt, 'x0':x0, 'xopt':xopt, 'T':T, 'error':error, 'error0':error0}


def heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution=100, limits=(), errorFunc=sse, mapFunc=chronifyBeta, method='full', x0=(0,1,1), ftol=1e-5, scale=False):
	import scipy.optimize as so
	
	def heterochronyObjFunc(parms, T, Y, Z, knots=[]): 
		invparms = (parms[0], 1/parms[1], 1/parms[2])
		return errorFunc(Z, mapFunc(T, Y, parms, knots))+errorFunc(Y, mapFunc(T, Z, invparms, knots))
	
	def heterochronyObjFunc2(parms, gamma, T, Y, Z, knots=[]): 
		tparms = [gamma, parms[0], parms[1]]
		invparms = (gamma, 1/parms[0], 1/parms[1])
		return errorFunc(Z, mapFunc(T, Y, tparms, knots))+errorFunc(Y, mapFunc(T, Z, invparms, knots))
	
	
	# Prepare data
	# ------------
	T = X; Y0 = Y1; Yt = Y2
	
	T = makeRange(min(T),max(T),n=resolution)
	Y0 = timeTransform(X, Y0, period, perIdx, x0, knots=knots, newX=T, periodify=True)
	Yt = timeTransform(X, Yt, period, perIdx, x0, knots=knots, newX=T, periodify=True)
	
	# Y0 = timeTransform(X, Y0, period, perIdx, (), knots=knots)
	# Yt = timeTransform(X, Yt, period, perIdx, (), knots=knots)
	# ------------
	
	# joint estimation search problem for the 2 symmetric cases
	# -------------------------------------------------------
	# estimate the 3 parameters such that these and their inverses yield the lowest total error
	args = (T, standardize(Y0), standardize(Yt), knots)
	err0 = heterochronyObjFunc(x0, *args)
	xopt = None; fopt = None
	
	# downhill simplex algorithm on shape, then grid phase with final fmin
	if method == 'brute':
		if not limits: limits = [(None,None) for i in range(len(x0))]
		args = (x0[0], T, Y0, Yt, knots)
		
		results = so.brute(heterochronyObjFunc2, ranges=limits, args=args, full_output=1)
		(xopt, fopt, grid, Jout) = results
		xopt = [x0[0], xopt[0], xopt[1]]
		
	elif method == 'simplex':
		results = so.fmin(heterochronyObjFunc, x0, args=args, 
		                  full_output=1, disp=0, ftol=ftol)
		(xopt, fopt, iter, funcalls, warnflag) = results
	else:
		
		# first local optimal search
		results = so.fmin(heterochronyObjFunc, x0, args=args, 
		                  full_output=1, disp=0, ftol=ftol)
		(xopt, fopt, iter, funcalls, warnflag) = results
		
		errat = fopt/err0
		
		# try random search
		import random, math
		fib = [1,1,2,3,5,8,13,21,34,55,89,144,233]#,377]
		nrand = map(multna1(10),fib)
		searchcutoff = .4
		KEEPSEARCHING = errat >= searchcutoff; count = 0; ctier = 0
		
		while KEEPSEARCHING and ctier < len(nrand) and count < nrand[ctier]:
			gamma = random.uniform(-133.5,133.5)
			alpha = math.e**random.uniform(-1.61,1.61)
			beta = math.e**random.uniform(-1.61,1.61)
			xc = (gamma, alpha, beta)
			
			error = heterochronyObjFunc(xc, *args)
			if abs(error/err0 - fopt/err0) <= .1: 
				results = so.fmin(heterochronyObjFunc, xc, args=args, 
				                  full_output=1, disp=0, ftol=ftol)
				(txopt, tfopt, iter, funcalls, warnflag) = results
				if tfopt < fopt: fopt = tfopt; xopt = txopt
			
			count += 1
			
			# print ' - tier', ctier, 'count', count, 'continue', KEEPSEARCHING
			if count == nrand[ctier]:
				ctier += 1; count = 0
				# adaptive - start harsh 
				if fopt/err0 == errat: searchcutoff += (1-searchcutoff)/float(len(nrand))
				KEEPSEARCHING = fopt/err0 >= searchcutoff
				# print 'tier', ctier, 'count', count, 'continue', KEEPSEARCHING
	# -------------------------------------------------------
	
	xopt = map(r3, xopt) # round
	xoptinv = map(r3, [-xopt[0],1/xopt[1],1/xopt[2]])
	
	# Yopt = timeTransform(X, Y0, period, perIdx, xopt, knots=knots)
	# Yoptinv = timeTransform(X, Yt, period, perIdx, xoptinv, knots=knots)
	
	Yopt = mapFunc(T, Y0, xopt, knots)
	Yoptinv = mapFunc(T, Yt, xoptinv, knots)
	
	
	if scale:
		Yopt = standardize(Yopt); Yoptinv = standardize(Yoptinv)
		Y0 = standardize(Y0); Yt = standardize(Yt)
		
	# now compute final fit of this model
	# --------------------------------------
	# have parms - need regression in both directions
	
	res = regression(Y0, Yt, type='basic')
	err_h0 = res['sse']; r2_h0 = res['R2']
	
	res = regression(Yopt, Yt, type='basic')
	Aerr_h1 = res['sse']; Ar2_h1 = res['R2']
	An = res['n']
	res = regression(Yoptinv, Y0, type='basic')
	Berr_h1 = res['sse']; Br2_h1 = res['R2']
	Bn = res['n']
	
	res_n = An+Bn # data points involved
	terr_h0 = err_h0; tr2_h0 = r2_h0 # i.e., 2*err_h0/2
	
	# average R2 is equivalent to the sum of SSRs / sum of sstot
	terr_h1 = (Aerr_h1 + Berr_h1)/2.
	tr2_h1 = (Ar2_h1 + Br2_h1)/2.
	
	return {'Y0':Y0, 'Ytarget':Yt, 'Yopt':Yopt, 'Yoptinv':Yoptinv, 'x0':x0, 'xopt':xopt, 'xoptinv':xoptinv, 'T':T, 
			'SSE_h0':err_h0, 'R2_h0':r2_h0, 'n':res_n,
			'SSE_h1':Aerr_h1, 'R2_h1':Ar2_h1, 
			'SSE_h1inv':Berr_h1, 'R2_h1inv':Br2_h1,
			
			'SSE_tot_h0':terr_h0, 'R2_tot_h0':tr2_h0,
			'SSE_tot_h1':terr_h1, 'R2_tot_h1':tr2_h1
			}



def heterochrony_optimization_complete(X, Y1, Y2, knots, period, perIdx, resolution=None, limits=(), errorFunc=sse, mapFunc=chronifyBeta, x0=(0,1,1), ftol=1e-7, bruteLimits=(2,2,2), scale=False, nreps=5):
	"""Wrapper function"""
	import numpy as np
	
	# this is the problem: even if given phase, how to identify this as optimal?
	
	# given phase, estimate shape by brute force?
	# simplex doesn't do it.
	
	# fix timeTransform
	# limits = np.s_[-45:45.1:1, 1/5.:5.+.1:.1, 1/5.:5.+.1:.1]
	
	limits = np.s_[1/5.:5.+.1:.1, 1/5.:5.+.1:.1]
	res = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
				limits, errorFunc, mapFunc, 'brute', x0, ftol, scale)
	
	# res = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
				# limits, errorFunc, mapFunc, 'simplex', x0, ftol, scale)
	
	if 0:
		# first search for fmin around each phase
		# res = {'SSE_tot_h1':GLOB_INF}
		gamma = res['xopt'][0]
		bounds = [max(-133.5,gamma-(133.5/2.)), min(133.5,gamma+(133.5/2.))]
		for gamma in makeRange(bounds[0],bounds[1],n=11): # n=133, 2 minute resolution
			xc = [gamma,x0[1],x0[2]]
			tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
						limits, errorFunc, mapFunc, 'simplex', xc, ftol, scale)
			if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: res = tmp
		
		# choice: if improvement then repeat at finer resolution
		# else, do random parameters - 'full'
		
		# repeat at finer resolution around gamma
		gamma = res['xopt'][0]
		bounds = [max(-133.5,gamma-(133.5/2.)), min(133.5,gamma+(133.5/2.))]
		print 'first round', res['xopt'], 'bounding', bounds
	
		for gamma in makeRange(bounds[0],bounds[1],n=11): # n=133, 2 minute resolution
			xc = [gamma,x0[1],x0[2]]
			tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
						limits, errorFunc, mapFunc, 'simplex', xc, ftol, scale)
			if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: res = tmp
	
		gamma = res['xopt'][0]
		bounds = [max(-133.5,gamma-(133.5/2.)), min(133.5,gamma+(133.5/2.))]
		print 'second round', res['xopt'], 'bounding', bounds
	
		for gamma in makeRange(bounds[0],bounds[1],n=11): # n=133, 2 minute resolution
			xc = [gamma,x0[1],x0[2]]
			tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
						limits, errorFunc, mapFunc, 'simplex', xc, ftol, scale)
			if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: res = tmp
		print 'third round', res['xopt']
	
	
	
	# each rep takes worst time 1 minute
	try: errat = res['SSE_tot_h1']/res['SSE_tot_h0']
	except ZeroDivisionError: errat = 50
	# errat = 0
	if errat >= bruteLimits[0] and errat <= bruteLimits[1]:
		NEWER = 0 # how many iterations ago did we improve?
		for i in range(nreps):
			# print 'opt iter', i, res['SSE_tot_h1'], res['R2_tot_h1']
			tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
						limits, errorFunc, mapFunc, False, x0, ftol, scale)
			NEWER += 1
			if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: 
				res = tmp
				errat = res['SSE_tot_h1']/res['SSE_tot_h0']
				NEWER = 0
				if errat < bruteLimits[0]: break
			if NEWER >= 1: break
	
	elif bruteLimits[3] == True and errat >= bruteLimits[1] and errat <= bruteLimits[2]:
		print '=> BRUTE: possible local trap %.2f vs %2.f (%.2f) xopt'%(res['SSE_tot_h1'], res['SSE_tot_h0'], errat)
		tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution, 
					limits, errorFunc, mapFunc, BRUTE=True)
		if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: res = tmp
	
	return res



def heterochrony_optimization_complete_rand(X, Y1, Y2, knots, period, perIdx, resolution=None, limits=(), errorFunc=sse, mapFunc=chronifyBeta, x0=(0,1,1), ftol=1e-5, bruteLimits=(2,2,2), scale=False, nreps=5):
	"""Wrapper function"""
	
	# i took the hardest search i could find, and fount P(max) = .1
	res = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
				limits, errorFunc, mapFunc, False, x0, ftol, scale)
	# each rep takes worst time 1 minute
	errat = res['SSE_tot_h1']/res['SSE_tot_h0']
	if errat >= bruteLimits[0] and errat <= bruteLimits[1]:
		NEWER = 0 # how many iterations ago did we improve?
		for i in range(nreps):
			# print 'opt iter', i, res['SSE_tot_h1'], res['R2_tot_h1']
			tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution,
						limits, errorFunc, mapFunc, False, x0, ftol, scale)
			NEWER += 1
			if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: 
				res = tmp
				errat = res['SSE_tot_h1']/res['SSE_tot_h0']
				NEWER = 0
				if errat < bruteLimits[0]: break
			if NEWER >= 1: break
	
	elif bruteLimits[3] == True and errat >= bruteLimits[1] and errat <= bruteLimits[2]:
		print '=> BRUTE: possible local trap %.2f vs %2.f (%.2f) xopt'%(res['SSE_tot_h1'], res['SSE_tot_h0'], errat)
		tmp = heterochrony_optimization(X, Y1, Y2, knots, period, perIdx, resolution, 
					limits, errorFunc, mapFunc, BRUTE=True)
		if tmp['SSE_tot_h1'] < res['SSE_tot_h1']: res = tmp
	
	return res


# OLD!
def heterochrony_optimization_lastgoodone(X, Y1, Y2, knots, period, perIdx, resolution=100, ftol=1e2, limits=(), stages=(0), stagedLimits=0, errorFunc=sse, approx_grad=True, hires=False, mapFunc=chronifyBeta):
	# from scipy import *
	import numpy as np
	import scipy.linalg as la, scipy.optimize as so
	from scipy.interpolate import splev, splrep
	
	# def heterochronyObjFunc(parms, T, Y, Z, knots=[], split=0): 
	# 	return errorFunc(Z, mapFunc(T, Y, parms, knots, split))
	
	def heterochronyObjFunc(parms, T, Y, Z, knots=[]): 
		return errorFunc(Z, mapFunc(T, Y, parms, knots))
	
	# x0 = (1,0) # initial guess, linear mapping (exp, phase)
	x0 = (0,1,1) # beta, linear mapping, no phase
	if not limits: limits = [(None,None) for i in range(len(x0))]
	
	# if mapFunc == chronifyExp2point: 
	# 	x0 = (0,1,1) # (split,phase, exp1, exp2)
	# 	if not limits: limits = [(-period/2.,period/2.),(0,None),(0,None)]
	
	# need to perform various manipulations of the data
	# -------------------------------------------------
	Z1 = Y1[:]; Z2 = Y2[:]  # copy data
	
	Tpre = X[:]
	Torig = Tpre[0:perIdx] # correct time
	Tuse = modulate(period)(Tpre)
	
	# make periodic
	Z1per = makePeriodic(Z1, perIdx)
	Z2per = makePeriodic(Z2, perIdx)
	
	# need some number of discrete samples between min and max T
	start = min(Torig)
	step = period / float(resolution)
	step += step/float(resolution)
	T = [start + i*step for i in range(resolution)]
	Thires = T[:]
	
	# re-sample at reasonable resolution for smooth optimization
	window = 3; iters = 1
	Y0 = resampleIter(Torig, Z1per, T, window=window, iterations=iters, per=1, knots=knots)
	Ytarget = resampleIter(Torig, Z2per, T, window=window, iterations=iters, per=1, knots=knots)
	
	# standardize so that optimization fits the pattern
	Y0s = standardize(Y0)
	Ytargets = standardize(Ytarget)
	args = (T, Y0s, Ytargets, knots)
	
	# unconstrained
	# results = so.fmin(heterochronyObjFunc, x0, args=args, 
	#                   full_output=1, disp=0, ftol=ftol)
	# (xopt, fopt, iter, funcalls, warnflag) = results
	
	# constrained
	results = so.fmin_l_bfgs_b(heterochronyObjFunc, x0, args=args, factr=ftol, bounds=limits, approx_grad=True, iprint=-1)
	(xopt, fopt, resdct) = results
	# error0 = heterochronyObjFunc(xopt, T, Y0, Ytarget, knots)
	# error = heterochronyObjFunc(xopt, T, Y0, Ytarget, knots)
	
	# constrained
	# xopt = x0; split = 0
	# if mapFunc == chronifyExp2point:
	# 	xopt = x0; error = GLOB_INF
	# 	for spl in stages:
	# 		args = (T, Y0, Ytarget, knots, spl)
	# 		
	# 		# stage-specific limits
	# 		if stagedLimits > 0:
	# 			f1 = spl*stagedLimits; f2 = (1-spl)*stagedLimits
	# 			if f1 == 0: f1 = f2
	# 			elif f2 == 0: f2 = f1
	# 			climits = (limits[0], (1/f1, f1), (1/f2, f2))
	# 			# print spl, 'climit', climits
	# 		else: climits = limits
	# 		results = so.fmin_l_bfgs_b(heterochronyObjFunc, x0, args=args, factr=ftol, bounds=climits, approx_grad=True, iprint=-1)
	# 		(copt, fopt, resdct) = results
	# 		cerr = heterochronyObjFunc(copt, T, Y0, Ytarget, knots, spl)
	# 		if cerr < error:
	# 			error = cerr; xopt = copt; split = spl
	# else:
	# 	args = (T, Y0, Ytarget, knots)
	# 	results = so.fmin_l_bfgs_b(heterochronyObjFunc, x0, args=args, factr=ftol, bounds=limits, approx_grad=True, iprint=-1)
	# 	(xopt, fopt, resdct) = results
	# 	error = heterochronyObjFunc(xopt, T, Y0, Ytarget, knots)
	
	# transform original expression data with null and opt parameters
	Y0 = mapFunc(T, Y0, x0, knots)
	Yopt = mapFunc(T, Y0, xopt, knots)
	Ytarget = mapFunc(T, Ytarget, x0, knots)
	# Yopt = mapFunc(T, Y0, xopt, knots, split)
	
	# if optimization is worse than original, return original
	error = lsq(Yopt,Ytarget); error0 = lsq(Y0, Ytarget)
	# if error > error0:
		# xopt = x0; Yopt = Y0; error = error0
	
	# return new Y1, Y2
	# xopt = [split]+xopt.tolist()
	
	if not hires: 
		# resample these to the original (periodic) timepoints
		Y0s = resampleIter(T, Y0, Tuse, window=3, iterations=1, per=1, knots=knots)
		Yopts = resampleIter(T, Yopt, Tuse, window=3, iterations=1, per=1, knots=knots)
		Ytargets = resampleIter(T, Ytarget, Tuse, window=3, iterations=1, per=1, knots=knots)
		return {'Y0':Y0s, 'Yopt':Yopts, 'Ytarget':Ytargets, 'x0':x0, 'xopt':xopt, 'T':Tpre, 'error':error, 'error0':error0}
	else: 
		return {'Y0':Y0, 'Yopt':Yopt, 'Ytarget':Ytarget, 'x0':x0, 'xopt':xopt, 'T':T, 'error':error, 'error0':error0}



def heterochrony_optimization_old(X, Y1, Y2, order, knots, period, perIdx, resolution, ftol=1e-8, limits=(), approx_grad=False, keepres=False, type='phase2', errorFunc=lsq):
	# from scipy import *
	import numpy as np
	import scipy.linalg as la, scipy.optimize as so
	from scipy.interpolate import splprep, splev, splrep
	
	# ftol: tolerance of optimization function
	
	# helper methods
	# -------------------------
	skewFunc = []
	if type == 'phase2':
		skewFunc += [lambda Z: lambda t: t] 
		skewFunc += [lambda Z: lambda t: Z[0]*t] 
		skewFunc += [lambda Z: lambda t: Z[0]*t + Z[1]] 
		skewFunc += [lambda Z: lambda t: Z[0]*t + Z[1] + Z[2]*t**2] 
	else:
		skewFunc += [lambda Z: lambda t: t] 
		skewFunc += [lambda Z: lambda t: Z[0] + t] 
		skewFunc += [lambda Z: lambda t: Z[0] + Z[1]*t] 
		skewFunc += [lambda Z: lambda t: Z[0] + Z[1]*t + Z[2]*t**2] 
	
	def heterochronyObjFunc(parms, T, Y, Z, knots=[]): 
		return errorFunc(Z, chronify(T, Y, parms, knots))
		# return wlsq(Z, chronify(T, Y, parms, knots))
		# return 10*(1-cor(Z, chronify(T, Y, parms, knots)))
	
		# -------------------------
	
	x0 = (0,1,0) # initial guess, quadratic skew
	if type == 'phase2':
		if order == 0: x0 = ()
		elif order == 1: x0 = (1)
		elif order == 2: x0 = (1,0)
		elif order == 3: x0 = (1,0,0)
	else:
		if order == 0: x0 = ()
		elif order == 1: x0 = (0)
		elif order == 2: x0 = (0,1)
		elif order == 3: x0 = (0,1,0)
	
	
	# need to perform various manipulations of the data
	# -------------------------------------------------
	Z1 = Y1[:]; Z2 = Y2[:]  # copy data
	
	Tpre = X[:]
	Torig = Tpre[0:perIdx] # correct time
	Tuse = modulate(period)(Tpre)
	
	# make periodic
	Z1per = makePeriodic(Z1, perIdx)
	Z2per = makePeriodic(Z2, perIdx)
	
	# need some number of discrete samples between min and max T
	start = min(Torig)
	step = period / float(resolution)
	step += step/float(resolution)
	T = [start + i*step for i in range(resolution)]
	Thires = T[:]
	# print 'period', period, 'resolution', resolution, 'step', step, 'last', tail(T)
	
	# re-sample at reasonable resolution for smooth optimization
	window = 3; iters = 1
	Z1 = resample(Torig, Z1per, T, window=window, iterations=iters, per=1, knots=knots)
	Z2 = resample(Torig, Z2per, T, window=window, iterations=iters, per=1, knots=knots)
	
	# standardize
	Y0 = standardize(Z1)
	Ytarget = standardize(Z2)
	
	
	args = (T, Y0, Ytarget, knots)
	
	# unconstrained
	# results = so.fmin(heterochronyObjFunc, x0, args=args, 
	#                   full_output=1, disp=0, ftol=ftol)
	# (xopt, fopt, iter, funcalls, warnflag) = results
	
	# constrained
	# if type(ftol) == type(list):
	# 		pgtol
	results = so.fmin_l_bfgs_b(heterochronyObjFunc, x0, args=args, factr=ftol, 
	                           approx_grad=approx_grad, bounds=limits)
	(xopt, fopt, resdct) = results
	
	error = heterochronyObjFunc(xopt, T, Y0, Ytarget, knots)
	
	# now we have the parameters
	# we need to again resample the data, but at high resolution this time
	# for smooth heterochronic curve
	resolution = period # changing the resolution
	
	# apply the transformation, and then subsample and standardize 
	# to recover the proper data
	start = min(Torig)
	step = period / float(resolution)
	step += step/float(resolution)
	T = [start + i*step for i in range(resolution)]
	
	Y0 = resample(Torig, Z1per, T, window=3, iterations=1, per=1, knots=knots)
	# transform the query gene using transformation parameters
	Yopt = chronify(T, Y0, xopt, knots)
	Ytarget = resample(Torig, Z2per, T, window=3, iterations=1, per=1, knots=knots)
	
	# now sub sample to original 18 times
	Y0s = subSample(Y0, Tuse)
	Yopts = subSample(Yopt, Tuse)
	Ytargets = subSample(Ytarget, Tuse)
	
	# standardize data for comparison
	Y0 = standardize(Y0)
	Yopt = standardize(Yopt)
	Ytarget = standardize(Ytarget)
	
	Y0s = standardize(Y0s)
	Yopts = standardize(Yopts)
	Ytargets = standardize(Ytargets)
	
	# return new Y1, Y2
	if not keepres: return {'Y0':Y0s, 'Yopt':Yopts, 'Ytarget':Ytargets, 'xopt':xopt, 'T':Tuse, 'error':error}# Y0s, Yopts, Ytargets, xopt, Tuse
	else: return {'Y0':Y0, 'Yopt':Yopt, 'Ytarget':Ytarget, 'xopt':xopt, 'T':T, 'error':error}


def heterochrony_optimization_choice(X, Y1, Y2, order, knots, period, perIdx, resolution, ftol=1e-8, search_method='simplex', limits=(), Nsamples=20, keepres=False, type='phase2'):
	# from scipy import *
	import numpy as np
	import scipy.linalg as la, scipy.optimize as so
	from scipy.interpolate import splprep, splev, splrep
	
	# ftol: tolerance of optimization function
	
	# helper methods
	# -------------------------
	skewFunc = []
	if type == 'phase2':
		skewFunc += [lambda Z: lambda t: t] 
		skewFunc += [lambda Z: lambda t: Z[0]*t] 
		skewFunc += [lambda Z: lambda t: Z[0]*t + Z[1]] 
		skewFunc += [lambda Z: lambda t: Z[0]*t + Z[1] + Z[2]*t**2] 
	else:
		skewFunc += [lambda Z: lambda t: t] 
		skewFunc += [lambda Z: lambda t: Z[0] + t] 
		skewFunc += [lambda Z: lambda t: Z[0] + Z[1]*t] 
		skewFunc += [lambda Z: lambda t: Z[0] + Z[1]*t + Z[2]*t**2] 
	
	lsq = lambda x,y: sum(map(lambda a,b: (a-b)**2, x,y))
	wlsq = lambda x,y: sum(map(lambda a,b: (1/float(abs(a-b)))*(a-b)**2, x,y))
	def heterochronyObjFunc(parms, T, Y, Z, knots=[]): 
		return lsq(Z, chronify(T, Y, parms, knots))
		# return wlsq(Z, chronify(T, Y, parms, knots))
		# return 10*(1-cor(Z, chronify(T, Y, parms, knots)))
	
		# -------------------------
	
	x0 = (0,1,0) # initial guess, quadratic skew
	if type == 'phase2':
		if order == 0: x0 = ()
		elif order == 1: x0 = (1)
		elif order == 2: x0 = (1,0)
		elif order == 3: x0 = (1,0,0)
	else:
		if order == 0: x0 = ()
		elif order == 1: x0 = (0)
		elif order == 2: x0 = (0,1)
		elif order == 3: x0 = (0,1,0)
	
	
	# need to perform various manipulations of the data
	# -------------------------------------------------
	Z1 = Y1[:]; Z2 = Y2[:]  # copy data
	
	Tpre = X[:]
	Torig = Tpre[0:perIdx] # correct time
	Tuse = modulate(period)(Tpre)
	
	# make periodic
	Z1per = makePeriodic(Z1, perIdx)
	Z2per = makePeriodic(Z2, perIdx)
	
	# need some number of discrete samples between min and max T
	start = min(Torig)
	step = period / float(resolution)
	step += step/float(resolution)
	T = [start + i*step for i in range(resolution)]
	Thires = T[:]
	# print 'period', period, 'resolution', resolution, 'step', step, 'last', tail(T)
	
	# re-sample at reasonable resolution for smooth optimization
	window = 3; iters = 1
	Z1 = resample(Torig, Z1per, T, window=window, iterations=iters, per=1, knots=knots)
	Z2 = resample(Torig, Z2per, T, window=window, iterations=iters, per=1, knots=knots)
	
	# standardize
	Y0 = standardize(Z1)
	Ytarget = standardize(Z2)
	
	# estimate transformation parameters given order
	args = (T, Y0, Ytarget, knots)
	
	# results = so.fmin(heterochronyObjFunc, x0, args=args, 
	#                   full_output=1, disp=0, ftol=ftol)
	# (xopt, fopt, iter, funcalls, warnflag) = results
	
	# choice of search methods
	# ------------------------
	if search_method == 'simplex':
		# downhill simplex algorithm
		results = so.fmin(heterochronyObjFunc, x0, args=args, 
		                  full_output=1, disp=0, ftol=ftol)
		(xopt, fopt, iter, funcalls, warnflag) = results
	elif search_method == 'fmin_bfgs':
		# Broyden-Fletcher-Goldfarb-Shanno algorithm
		results = so.fmin_bfgs(heterochronyObjFunc, x0, args=args, 
		                  full_output=1, disp=0)
		(xopt, fopt, iter, funcalls, warnflag) = results
		
	elif search_method == 'fmin_l_bfgs_b':
		xopt, fopt, dct = so.fmin_l_bfgs_b(heterochronyObjFunc, x0, args=args,
		approx_grad=True, bounds=limits)
		
	elif search_method == 'scipy-brute':
		# scipy's brute force
		results = so.brute(heterochronyObjFunc, limits, args=args, 
		                   full_output=1, Ns=Nsamples)
		(xopt, fopt, grid, Jout) = results
	elif search_method =='anneal':
	# simulated annealing
		results = so.anneal(heterochronyObjFunc, x0, 
		                      args=args, full_output=1,
		                      schedule='fast', boltzmann=1, dwell=50, Tf=1)
		(xopt, fopt, temp, feval, iters, accept, flag) = results
	elif search_method == 'brute':
		# my brute force
		(xopt, fopt) = bruteForceSearch(heterochronyObjFunc, limits, args=args)
	else: sys.exit('bad search method.')
	# ------------------------
	
	error = heterochronyObjFunc(xopt, T, Y0, Ytarget, knots)
	
	# now we have the parameters
	# we need to again resample the data, but at high resolution this time
	# for smooth heterochronic curve
	resolution = period # changing the resolution
	
	# apply the transformation, and then subsample and standardize 
	# to recover the proper data
	start = min(Torig)
	step = period / float(resolution)
	step += step/float(resolution)
	T = [start + i*step for i in range(resolution)]
	
	Y0 = resample(Torig, Z1per, T, window=3, iterations=1, per=1, knots=knots)
	# transform the query gene using transformation parameters
	Yopt = chronify(T, Y0, xopt, knots)
	Ytarget = resample(Torig, Z2per, T, window=3, iterations=1, per=1, knots=knots)
	
	# now sub sample to original 18 times
	Y0s = subSample(Y0, Tuse)
	Yopts = subSample(Yopt, Tuse)
	Ytargets = subSample(Ytarget, Tuse)
	
	# standardize data for comparison
	Y0 = standardize(Y0)
	Yopt = standardize(Yopt)
	Ytarget = standardize(Ytarget)
	
	Y0s = standardize(Y0s)
	Yopts = standardize(Yopts)
	Ytargets = standardize(Ytargets)
	
	# return new Y1, Y2
	if not keepres: return {'Y0':Y0s, 'Yopt':Yopts, 'Ytarget':Ytargets, 'xopt':xopt, 'T':Tuse, 'error':error}# Y0s, Yopts, Ytargets, xopt, Tuse
	else: return {'Y0':Y0, 'Yopt':Yopt, 'Ytarget':Ytarget, 'xopt':xopt, 'T':T, 'error':error}


# this needs to be updated with new model
def het_tests(X, Y1, Y2, order, knots, period, perIdx, resolution, ftol=1e-8, keepres=False):
	"""Performs a set of nested hypothesis tests about the order complexity of
	the transformation relationship between two temporal trajectories."""
	
	const = 2
	R2 = []
	np = []
	df = []
	N = []
	F = []
	P = []
	parms = []
	for i in range(order+1):
		
		# prepare the temporal transformation given current order
		args = (i,knots,period,perIdx,resolution,ftol,keepres)
		Z0, Z1, Z2, xopt, Tused = heterochrony_optimization(X,Y1,Y2,*args)
		ret = regression(Z1,Z2)
		np += [const+i]
		N += [ret['n']]
		R2 += [ret['R2']]
		parms += [xopt.tolist()]
		# perform nested F tests
		if i == 0:
			# significant regression?
			F += [ret['F']]
			P += [ret['P']]
			df += [[ret['df top'],ret['df bottom']]]
		else:
			# relative improvement test
			df1 = float(np[i]-np[i-1])
			df2 = float(N[i]-np[i]-1)
			F += [((R2[i] - R2[i-1])/df1) / ((1-R2[i])/df2)]
			# looking for large F, so reverse
			P += [1.-FdistCDF(F[i], df1, df2)]
			df += [[df1,df2]]
			
	
	i = order
	jump = {}
	if i == 0:
		jdf1 = df[i][0]
		jdf2 = df[i][1]
		jF = F[i]
		jP = P[i]
		jump = {'F':jF, 'P':jP, 'df1':jdf1, 'df2':jdf2}
	else:
		# one final jump-up F value
		jdf1 = float(np[i]-np[0]) # delta dof for numerator
		jdf2 = float(N[i]-np[i]-1)
		jF = ((R2[i] - R2[0])/jdf1) / ((1-R2[i])/jdf2)
		jP = 1.-FdistCDF(jF, jdf1, jdf2)
		jump = {'F':jF, 'P':jP, 'df1':jdf1, 'df2':jdf2, 'R2':R2[i], 'R2_H0':R2[0]}
	
	# and a LLR evaluation
	jumpllr = {}
	if i == 0:
		res = regression(Z0, Z2, type='basic')
		err_h0 = res['sse']; r2_h0 = res['R2']
		jumpllr = {'LLR':0, 'P':1, 'df':0, 'R2':nastr, 'R2_H0':r2_h0}
	else:
		res = regression(Z0, Z2, type='basic')
		err_h0 = res['sse']; r2_h0 = res['R2']
		res = regression(Z1, Z2, type='basic')
		err_h1 = res['sse']; r2_h1 = res['R2']
		
		llr = err_h0 - err_h1
		rdf = np[i]-np[0]
		x2p = llr == abs(llr) and ChiSquarePDF(llr, rdf) or ChiSquarePDF(0, rdf)
		jumpllr = {'LLR':llr, 'P':x2p, 'df':rdf, 'R2':r2_h1, 'R2_H0':r2_h0}
	
	# organize results
	tuples = []
	for i in range(len(F)):
		tuples += [{'F':F[i], 'P':P[i], 'R2':R2[i], 'df':df[i], 'N':N[i], 'parms':parms[i]}]
	# print tuples
	results = dict([[np[i]-const,tuples[i]] for i in range(len(F))])
	results['F'] = jump
	results['llr'] = jumpllr
	results['T'] = Tused
	results['orig'] = Z0
	results['opt'] = Z1
	results['target'] = Z2
	return results


###########################

if __name__ == "__main__":
	### UNIT TESTING ###
	####################
	# print 'mean', 100000/500.
	ncells = 500000
	withmark = .25
	nwm = ncells*withmark
	print nwm, 'binomial', 1-binomialPDF(nwm,200,1/500.,approx=1)
	sys.exit()
	
	# compute probability of sequencing the same DNA template more than once
	cycles = 18
	cells = 50e6
	reads = 50e6
	
	
	print 'fact n', factorial(1), factorial(10), factorial(100), factorial(1000)
	sys.exit()
	
	
	def compute(cycles,cells,reads):
		# this is probability of selecting a unique start site molecule
		# P = binomial(cells, 1, (2**cycles)/float(cells),approx=1)
		
		Pus = 1 - cells/float(cells*2**cycles)
		P = binomial(cells, 0, Pus, approx=1)
		print 'Pus', Pus, 'P', P
		# print 'Prob of seeing same molecule more than once', 1-P
		
		# this is prob of sequencing same start site more than once
		return 1 - binomialPDF(reads, 1, P, approx=1)
	
	
	print 'testing dual reads'
	# for cycles in [0,1,2,5,10,15,18,20,21,22,24,25]:
	for cycles in [0.1, 1, 2, 3, 4]:
		print 'cycles=%s, cells=%.2e, reads=%.2e => Pr(multiple reads)=%s'%(cycles,cells,reads, compute(cycles,cells,reads))
	print
	sys.exit()
	print 'Varying reads does not affect probability'
	for reads in [1e6,10e6,25e6,50e6,100e6]:
		cycles = 18
		print 'cycles=%s, cells=%.d, reads=%.2e => Pr(multiple reads)=%s'%(cycles,cells,reads, compute(cycles,cells,reads))
	print
	print 'Varying cells does not affect'
	for cells in [10e3,100e3,250e3,500e3,1e6,5e6,1e8,1e9]:
		cycles = 18
		reads = 50e6
		print 'cycles=%s, cells=%.d, reads=%.2e => Pr(multiple reads)=%s'%(cycles,cells,reads, compute(cycles,cells,reads))
	print
	print 'but need to vary ratio of cells to cycles'
	reads = 50e6
	cycles = 10; cells = 10e6
	print 'molecules=%.3e, cycles=%s, cells=%.3f, reads=%.2e => Pr(multiple reads)=%s'%(cells*2**cycles, cycles,cells,reads, compute(cycles,cells,reads))
	cycles = 1; cells = 10e6
	print 'molecules=%.3e, cycles=%s, cells=%.3f, reads=%.2e => Pr(multiple reads)=%s'%(cells*2**cycles, cycles,cells,reads, compute(cycles,cells,reads))
	
	sys.exit()
	
	
	n = 50
	k = 10
	p = .5
	
	# mu = n*p
	# sd = mu*(1-p)
	# print 'mu', mu, 'sd', sd
	
	# print 'test', sum([binomial(50,k,.5) for k in range(50)])
	# print 'test', sum([binomial(50,k,.5,1) for k in range(50)])
	import maths
	print 'test', binomial(n,k,p,approx=0), binomial(n,k,p,approx=1)#, maths.GaussianPDF(k,n*p,n*p*(1-p))
	print 'test', binomialPDF(n,k,p,approx=0), binomialPDF(n,k,p,approx=1)
	sys.exit()
	
	a = histogram([1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7], bins=[0,1,2,3,4,5,6,7,8,9,10])
	print 'a', a
	print 'new bins', [1,2,3]
	print reHist(a,bins=[1,2])
	
	# print binomialPDF(9,9,.7), binomial(9,9,.7)
	sys.exit()
	
	print 'test', ChiSquareCDF(8147.17, 10320)
	print 'test', ChiSquareCDF(1000,900), ChiSquarePDF(1000,900)
	sys.exit()
	
	AFo = hom1 = [10,1,1,1]
	AFe = hom2 = [25,2,2,2]
	
	het1 = [10,1,8,1]
	het2 = [13,1,7,1]
	het3 = [13,13,1,2]
	
	print '1 hom/hom'
	alleleTest(hom1, hom2, 'A', SEQ=False, VERBOSE=1)
	print '2 hom/hom'
	alleleTest(hom2, hom1, 'A', SEQ=False, VERBOSE=1)
	print 'hom/het'
	alleleTest(hom1, het1, 'A', SEQ=False, VERBOSE=1, direction=min)
	print 'het/hom'
	alleleTest(het1, hom1, 'A', SEQ=False, VERBOSE=1, direction=min)
	print '1 het/het'
	alleleTest(het1, het2, 'A', SEQ=False, VERBOSE=1)
	print '2 het/het'
	alleleTest(het2, het1, 'A', SEQ=False, VERBOSE=1)
	print '3 het/het*'
	alleleTest(het2, het3, 'A', SEQ=False, VERBOSE=1)
	alleleTest(het3, het2, 'A', SEQ=False, VERBOSE=1)
	
	sys.exit()
	
	print 'binom', binomialPDF(10,3,.5)
	sys.exit()
	
	a = [63,31,28,12,39,16,40,12]
	b = [67.78135, 22.59375, 22.59375, 7.53125, 45.18750, 15.06250, 45.18750, 15.06250]
	print 'X2', ChiSquareFreqTest(a,b), ChiSquarePDF(14.067, 7)
	sys.exit()
	
	for n in range(5,30):
		for d in range(-100,100,10):
			print d/10., n, TajimaDistPDF(d/10.,n)
	
	
	for n in [5,10,20,30]:
		for d in range(-50,50,10):
			print d/10., n, TajimaDistPDF(d/10.,n)
	
	sys.exit()
	
	
	
	
	print optimizeZ(0.05)
	print optimizeZ(0.05/2.)
	sys.exit()
	
	import random
	for i in range(100):
		print random.randint(0,44)
	sys.exit()
	
	# print 'multinomial', multinomial([.5,.5],[3,8])
	print multinomial([1/15. for i in range(15)],[6 for i in range(15)])
	
	
	print BetaDistCDF(.5,1,1)
	sys.exit()
	
	import sio as io
	lds = '/Users/simola/Desktop/time_domain/chronicFOCI-.3/LDS_modules/'
	files = io.getFiles(lds)
	for f in files:
		genes = io.readList(f)
		
		res = GOenrichment(genes)
		print 'sig', GOsig(res, alpha=.05, type='BH')
	sys.exit()
	
	dat = [ [10,5,4,2],[10,5,2,34],[110,200,3,3] ]
	print 'pvalues', batchFisherExactTest(dat)
	sys.exit()
	
	for i in range(100,110): print 1-FdistCDF(i,1,88)
	sys.exit()
	
	tab = [[147,9],[810,158]]
	print ChiSquareTest(tab)
	sys.exit()
	
	
	x = [0,24,48,63,87,111,135,152,176,194,218,227,251,260,284,301,325,345]
	y1 = [0.429561790397,0.005779913836,-0.610865582554,-0.0569982139213,-0.222701823544,-0.475650004298,-0.177082806469,-0.0699324367053,0.21111063154,-0.0570539193644,0.2290898954,-0.403391315057,-0.253786935942,-0.204859910677,-0.123269412586,-0.338414244498,-0.0573740345495,-0.238554286451]
	y2 = [-0.19985150453,-0.109426504756,-0.388747783111,-0.208569833308,-0.121576302385,-0.339192118717,-0.234488527784,-0.151191647877,-0.0476659828524,-0.0171644419154,-0.123176108417,-0.228510605629,-0.394963672192,-0.0643263009683,-0.577817563904,-0.114664847302,-0.0546669605207,-0.241756625377]
	
	order = 3
	knots = [24,48,63,87,111,135,152,176,194,218,227,251]
	period = 267
	perIdx = 14 # index closest to period, 260 min.
	resolution = 100 # number of interpolated time points for optimization
	
	args = (order,knots,period,perIdx,resolution)
	# z1,z2 = heterochrony_optimization(x,y1,y2,*args)
	results = het_tests(x,y1,y2,*args)
	print 'test', results[1]['R2']
	
	for order,rem in results.items():
		print 'order', order, 'results', rem
	sys.exit()
	
	
	print 'x2', ChiSquarePDF(-1, 10)
	
	print 't-tests'
	tval = -1.15
	df = 15
	print 't-test: 1', tdistPDF(tval, df=df, tails='1')
	print 't-test: 1', tdistPDF(-tval, df=df, tails='1')
	print 't-test: 2' , tdistPDF(tval, df=df, tails='2')
	print 't-test: 2' , tdistPDF(-tval, df=df, tails='2')
	
	
	
	lst = [1,2,3,4,2,2,7,4,5,8]
	print 'valus', lst
	print 'ranks', rank(lst)
	
	
	r = .86519
	n = 12
	tr = r * math.sqrt((n-2)/(1-r**2))
	# obtain a pvalue
	pval = tdistPDF(tr,n-2)
	print 'results', tr,n, {'r':r,'pvalue': pval}
	
	
	import maths
	
	print 'test ZdistPDF', maths.GaussianPDF(1.96, 0, 1)
	print 'my test ZdistPDF', GaussianPDF(2.0, 0., 1.)
	print 'tes gaussianCDF', maths.GaussianCDF(1.96, 0., 1.)
	
	
	a = [1,2,3,4,3,2,3,4,1,2,2,3,1,2,3]
	print 'mean with', meanWithConfidence(a)
	
	
	print 'test', maths.ChiSquareCDF(.484, 4)
	print 'optimize', optimizeX2(4, .05/2.0)
	
	print 'sokal test'
	ss = 13.52
	alpha = .05
	low, high = varianceConfidenceLimits(ss, 5, alpha=alpha)
	print ss, low, high
	print
	
	Ps = [maths.FdistCDF(12,8,22), FdistCDF(15,8,22,dt=1e-4)]
	print 'F = %f, Ps[0] = %f, Ps[1] = %f' % (1.36343708661e+24, Ps[0], Ps[1])
	# Ps = map(lambda x: 1-x, Ps)
	out = Bonferroni(Ps, range(len(Ps)), alpha=0.998611111111, n=4973,verbose=True)
	print 'before', 0.998611111111, 'actual', 1-(1-0.998611111111)/4973.
	print 'Bonferroni out', out
	print 'eh?', 1.0-(.05/36.0)
	
	print 'fishers test', fisherExactTest(0,10,12,2,24)
	
	
	print 'factorial comp'
	for i in range(10,20,5):
		print i, factorial(i), fact(i)
	
	testTajimaPDF = 0
	testFCDF = 0
	
	abi1 = [0.9,1,1.1,0.7,0.7,1.3,1.2,1,0.9,1.4,3.3,0.9,1,1.2,0.9,0.7,1.2,3.6,1,1,1,0.9,1,2.4,0.9,3.8,0.8,1.5,1,1,1.1,0.8,1.3,1.6,0.7,1.1,1,7,10]
	abi2 = [1.3,1.2,4.2,2.7,1.5,1.6,0.8,0.8,0.8,2.9,1.6,1.3,3.4,1.2,1.3,1.3,0.9,1.4,0.6,1.4,0.7,1.3,2,1,1.1,1.3,0.9,10]
	
	print 't-test'
	t,df= studentT(abi1,abi2, equalVariance=0, verbose=True)
	print 'df=%.3f, P(T<=t)=%f' % (df, tdistPDF(t,df, tails=2))
	
	####################
	L = [1,2,3,4,5,-5,6,7,1e10]
	i = argmin(L)
	print 'argmin', i, L[i]
	
	# j = permute([1,2,3,4,5,6,7])
	# print 'test', j
	
	# sys.exit()
	
	# print 9/10.0*avg([3,4,5,4,3,2,3,4,5]) + 1/10.0*avg([3])
	# print avg([3,4,5,4,3,2,3,4,5,3])
	import time_profiler as tp
	import maths
	
	X2 = 16800
	dof = 1634
	print 'cdf', maths.ChiSquareCDF(X2,dof)
	sys.exit()
	
	profiler = tp.timeprofile()
	profiler.mark('begin')
	
	print 'GO!'
	
	print maths.beta(.5,.5)
	print 'gamma', maths.gamma(.1), maths.gamma(.5), maths.gamma(1), maths.gamma(5)
	# print 'gamma2', maths.gammaRegularized(1,1), maths.gammaRegularized(2,1), maths.gammaRegularized(3,1)
	print 'Chi2', 1-maths.ChiSquareCDF(1.01, 1), maths.ChiSquareCDF(.85, 1), 1-maths.ChiSquareCDF(3.84, 1), 1-maths.ChiSquareCDF(6.63, 1)
	
	
	# for i in range(100): maths.beta(.5,.5)
	# print 'done C', profiler.lastdiff(); profiler.mark()
	
	# for i in range(10): Beta(.5,.5)
	# print 'done Python', profiler.lastdiff(); profiler.mark()
	
	# print 'Beta(.5,.5)', maths.beta(0.5,0.5), 'goal is 3.14159'
	# print 'Beta(.5,.5)', maths.incompleteBeta(1,0.5,0.5), 'goal is 3.14159'
	
	# print 'Beta(.5,.5)', Beta(.5,.5, .0000041)
	# print 'Fdist 3, 1, 1', 1-FdistCDF(3,10,10,.0000041)
	
	# print 'Fdist 3, 1, 1', 1-maths.FdistCDF(3,2,2)
	# print 'done', profiler.lastdiff(); profiler.mark()
	
	# Testing of Tajima D PDF p-values (cf Tajima 1989)
	if testTajimaPDF:
		a = 100.0
		for n in [5,10,20,30]:
			print 'n = ', n
			for mx in [-3,-2,-1,0,1,2,3,4]:
				pool = []
				# rnge = range(int(mx*a),int((mx-1)*a),1)
				rnge = range(int((mx-1)*a),int(mx*a),1)
				# print mx, mx*a, (mx-1)*a, min(rnge), max(rnge)
				for Db in rnge:
					D = Db/a
					phi = TajimaDistPDF(D,n)
					# print 'TajimaD PDF phi('+str(D)+',n='+str(n)+') = ',  phi
					pool.append(phi)
				print 'average D ['+str(mx)+':'+str((mx-1))+']', avg(pool)
	print 'done', profiler.lastdiff(); profiler.mark()
	# x = .1
	# a = 8
	# b = 24

	#D = 1.9
	#n = 5
	#print 'TajimaD PDF phi('+str(D)+',n='+str(n)+') = ',  TajimaDistPDF(D,n)
	
	# print 'App. Fdist('+str(x)+','+str(a)+','+str(b)+') = ', FdistCDFApprox(x,a,b)
	# print 'Real Fdist('+str(x)+','+str(a)+','+str(b)+') = ', FdistCDF(x,a,b)
	# print 'incompleteBeta('+str(x)+','+str(a)+','+str(b)+') = ', incompleteBeta(x,a,b)
	# print 'Beta('+str(a)+','+str(b)+') = ', Beta(a,b)
	# print 'regularizedIncompleteBeta('+str(x)+','+str(a)+','+str(b)+') = ', regularizedIncompleteBeta(x,a,b, dt=1e-4)
	
	# import sio
	# d,r,c = io.readMatrix('/Users/simola/Desktop/exprcrap/pvalues.txt', key=True)
	# dt = transpose(d)
	# a = dt[5]
	# b = r
	
	a = [.025, .05, .001, .005, .3, .4, .5]
	# a = a + [1-.025, 1-.05, 1-.001, 1-.005, 1-.3, 1-.4, 1-.5]
	a = a + [.025, .05, .001, .005, .3, .4, .5]
	# # a = [1-.05, 1-.025, 1-.01, 1-.001, 1-.0001, 1-.00001, 1-.000001]
	b = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'a2', 'b2', 'c2', 'd2', 'e2', 'f2', 'g2']
	print 'both', BH(a, b, alpha=.05, verbose=True)
	# print 'low', BH(a, b, alpha=.025, verbose=True)
	# print 'high', BH(a, b, alpha=1-.025, verbose=True)
	# Bonferroni(a, b, alpha=1-.001/2.0, verbose=True)


