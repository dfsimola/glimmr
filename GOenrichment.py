#! /usr/bin/env python

import re, sys, os, stat
from glimmrAccessories import *

version = """\nGOenrichment.py, version %s
 
 Requires Unix/Mac OS X/CYGWIN with Python 2.5+
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


# LOCAL METHODS
# ------------------------------------------------------------------

nastr = '"NA"'

roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))
r3 = lambda x: roundna(3)(x)

#
def unique(lst):
	"""Returns a subset of lst containing the unique elements."""
	d = {}
	for i in lst:
		if i not in d: d[i] = 0
		else: d[i] += 1
	e = d.keys(); e.sort()
	return e
	


#
def getLabel(astr):
	return '.'.join(astr.split('/')[-1].split('.')[:-1])
	


#
def getFiles(indirfile, include=[], exclude=[], literal=False):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	return getContents(indirfile, typ='files', include=include, exclude=exclude, literal=literal)


#
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


#
def slash(astr):
	if not len(astr): return '/'
	elif astr[-1] !='/': astr+='/'
	elif astr[len(astr)-2:len(astr)] == '//': astr = astr[:-1]
	return astr


#
def readList(fname, dtype=None, delim='\n'):
	fh = open(fname)
	lst = None
	if dtype: lst = map(lambda x: dtype(x[:-1]), fh)
	else: lst = map(lambda x: x[:-1], fh)
	fh.close()
	
	# lst = lst[0].split('\r')
	
	if delim != '\n' and len(lst)==1: return lst[0].split(delim)
	return lst


#
def printList(lst, title="", delim="\n", pipe=sys.stdout, file='', sort=False, newline='\n', ioMethod='w'):
	if file: pipe = open(file,ioMethod)
	if sort: lst.sort()
	if title: lst = [title]+lst
	for i in range(len(lst)):
		if i == len(lst)-1: pipe.write(str(lst[i]))
		else: pipe.write(str(lst[i])+delim)
	
	if newline: pipe.write("\n")
	
	if file: pipe.close()


#
def trimComments(dat):
	# parse comments
	ndat = []
	for line in dat:
		# if not line.strip(delim).strip('\n'): continue
		# only care about non commented lines
		if not re.match('(^\s*#.*$)|(^\s$)', line):
			# parse any commentary on this line
			cleanme = re.split('^(.*)#.*$', line)
			if len(cleanme) > 1: cleanme.pop(0)
			ndat.append(cleanme[0].strip('\n'))
	return ndat


#
def readTable(filename, header=True, rownames=True, delim="\t", comments=True, keepKey=False, fill=nastr):
	fh = open(filename)
	table = fh.readlines()
	fh.close()
	if not len(table): return [],[],[]
	if comments: table = trimComments(table)
	if not rownames: keepKey = True
	tabHeader = []; rowIDs = []; data = []
	if header: 
		tabHeader = table.pop(0).strip('\n').split(delim)
		if not keepKey: tabHeader.pop(0)
		# else: tabHeader = table.pop(0).strip('\n').split(delim)
	for line in table:
		sline = line.strip('\n').split(delim)
		if fill != None:
			sline = map(lambda x: x=='' and fill or x, sline)
		# print 'sline', sline
		if rownames:
			rowIDs.append(sline.pop(0))
		data += [sline]
	return data, rowIDs, tabHeader


#
def printTable(tab, header=[], rows=[], delim="\t", newline='\n', file='', pipe=sys.stdout, ioMethod='w'):
	
	if file: pipe = open(file,ioMethod)
	
	if header:
		printList(header, delim=delim, pipe=pipe, newline=newline)
	if len(rows):
		assert len(rows) == len(tab), "Rows and tab must have same length."
		for ri in range(len(tab)):
			r = [rows[ri]]+list(tab[ri])
			printList(r, delim=delim, pipe=pipe, newline=newline)
	else:
		for r in tab:
			printList(r, delim=delim, pipe=pipe, newline=newline)
	
	if file: pipe.close()


#
def createdir(dirname):
	try:
		parent = '/'.join(dirname.split('/')[:-1])
		if not dirname: return 0
		if not os.access(dirname, os.F_OK):
			os.mkdir(dirname)
			return 1
		else: return 0
	except OSError:
		print >> sys.stderr, 'sio.createdir(): Cannot access %s'%(parent)


#
def setDifference(a=[], b=[]):
	"""Setwise a - b: returns a list containing those items unique to a.
	"""
	# return filter(lambda item: item not in b, a)
	
	# faster version for big items
	bd = {}
	for item in b: bd[item] = 1
	return filter(lambda item: item not in bd, a)


#
def batchFisherExactTest(dat,tails=2,alpha=0.05, outdir='/var/tmp/'):
	"""Compute fisher exact pvalues for each of a list of [[a,b],[c,d]] tables or a,b,c,d tuples"""
	
	
	tmp = outdir+'rtempfile.txt'
	
	createdir(outdir)
	
	alt = 'two.sided'
	if tails == 1: alt = 'greater'
	elif tails == -1: alt = 'less'
	
	Rscpt = "r = c()\n"
	
	count = 0
	for tup in dat:
		
		# print count, tup
		
		fn = outdir+'tmpfishermat-%s.txt'%(count)
		printTable(tup,file=fn)
		
		Rscpt += """
		m = read.delim("%s",header=FALSE)
		ft = fisher.test(m,alternative=\""""%(fn)+alt+"""\",conf.level="""+str(1-float(alpha))+""")
		r = append(r, ft$p.value)
		"""
		count += 1
	Rscpt += """write.table(r, sep="\t", row.names=FALSE, col.names=FALSE,\""""+tmp+"""\")"""
	
	
	fh = open(outdir+'tmpscript.R','w')
	print >> fh, Rscpt
	fh.close()
	
	os.system('R --vanilla < "'+outdir+'tmpscript.R" > "'+outdir+'crap.out"')
	
	pvalues = []#; estimates = []
	pairs = []
	try: pairs = readList(tmp)
	except IOError: pass
	
	# print 'pairs', pairs, len(pairs)
	
	for p in pairs:
		# pval, est = p.split('\t')
		pvalues += [float(p)]
		# estimates += [float(est)]
	
	os.system("rm -f "+tmp)
	os.system('rm -f "%s"'%(outdir+'crap.out'))
	os.system('rm -f "%s"'%(outdir+'tmpscript.R'))
	
	count = 0
	for i in range(len(dat)):
		fn = outdir+'tmpfishermat-%s.txt'%(count)
		os.system('rm -f "%s"'%(fn))
		count += 1
	
	return pvalues#, estimates


#
def antGOenrichment(names, notnames=[], gofile='cflo.v3.3.UniProt.GO.txt', dirfile='./', verbose=False, splitnames=False):
	print '\nnum query genes', len(names), names[:5]
	# first parse the GO table into a dictionary of GO terms returning genes for each
	dct = {}
	allgenenames = []
	d,r,c = readTable(dirfile+gofile, header=False)
	for i in range(len(d)):
		gene = r[i]
		if splitnames: gene = gene.split('--')[0].strip()
		terms = d[i]
		allgenenames += [gene]
		for gobar in terms:
			for go in gobar.split('|'):
				if len(go)>1: # control for null and \n and \r stray lines
					try: dct[go] += [gene]
					except KeyError: dct[go] = [gene]
	
	# use genomic background?
	if not len(notnames): notnames = allgenenames
	
	# print len(names), len(notnames)
	if splitnames: 
		names = map(lambda x: x.split('--')[0].strip(), names)
		notnames = map(lambda x: x.split('--')[0].strip(), notnames)
	
	# unique list of names
	namesU = unique(names)
	notnamesU = unique(notnames)
	
	# print 'sample', namesU[:10], notnamesU[:10], allgenenames[:10]
	
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
		
		# a: genes of interest in group of interest
		thegenes = []
		groupnames = 0 # a
		for i in range(len(namesU)):
			if namesU[i] in fast: 
				groupnames += 1
				thegenes += [namesU[i]]
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
		ttab = [[groupnames, groupnotnames], [notgroupnames, notgroupnotnames]]
		dat += [ttab]
		metathegenes += [thegenes]
		metagroup += [group]
		
		if verbose: print 'group', group, ttab
		
	
	pvalues = batchFisherExactTest(dat, tails=1)
	
	if not(len(pvalues)): return results
	
	strs = []
	
	i = 0
	for group in GOterms:
		(groupnames, groupnotnames), (notgroupnames, notgroupnotnames) = dat[i]
		thegenes = metathegenes[i]
		
		pvalue = nastr
		try:
			pvalue = pvalues[i]
		except IndexError:
			print 'P-value error for %s: i/%s'%(group,len(pvalues))
		
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


#
eithereq = lambda s: lambda pair: (pair[0]!=s and pair[1]!=s and True)
filterstr = lambda s: lambda v: filter(eq(s), v)
pairedfilterstr = lambda s: lambda v: filter(eithereq(s), v)
#
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


#
def Bonferroni(P, L=[], alpha=0.05, verbose=False, n=-1):
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
	return [sigPairs[i][1] for i in range(len(sigPairs))]

#
def Pcutoff(P,L=[], alpha=.05):
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
		sigPs, sigLs = unzip(sigpairs)
		return sigLs
	else: return []


#
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


#
def fullGOanalysis(genelst, bglst=[], alpha=.05, method='Bonferroni', name='sample', savedir='', godir='', gofile='cflo.v3.3.UniProt.GO.txt', gotable='gene_ontology_ext.obo.txt', splitnames=False, verbose=True):
	
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
	
	
	createdir(savedir)
	# print 'HEY', godir+gotable
	akago = parseGOobo(godir+gotable)
	
	res = antGOenrichment(genelst, bglst, gofile=gofile,  dirfile=godir, splitnames=splitnames, verbose=verbose)
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
	printList(unique(flatgogenes), file=savedir+'all.GO.assoc.genes.txt')
	
	# save list of all genes that are sig GO associated
	flatgogenes = []
	for lst in sigg: flatgogenes += lst
	flatgogenes.sort()
	printList(unique(flatgogenes), file=savedir+'sig.GO.assoc.genes.txt')
	
	
	# save the significant genes to a directory
	gogenedir = savedir+'genes4sigTerms/'
	createdir(gogenedir)
	for s,g in zip(sig,sigg):
		printList(g,file=gogenedir+'%s_genes.txt'%(s))
	
	sigfile = savedir+'%s_GO_%s_%s.txt'%(name,method,alpha)
	gof = open(sigfile,'w')
	
	print '\nGO enrichment (n=%s, %s alpha < %s)'%(ngoterms,method,alpha), len(sig)
	print >>gof, '\n\nSignificant terms:'
	print >>gof, '\t'.join(['Name', 'GOID', 'P', 'Q(GO)', 'Q(FG)', 'Count'])
	
	# replace significant terms with their actual names
	keep = []
	for s,p,q,q2,c in zip(sig,ps,qs,q2s,cs):
		if s in akago: keep += [(p,[akago[s]['name'],s,p,r3(q),r3(q2),c])]
		else: keep += [(p,['NA',s,p,r3(q),r3(q2),c])]
	keep.sort() # most to least significant
	printTable([s for p,s in keep], pipe=gof)
	
	# sort all terms by pvalue
	allsort = []
	for k in res.iterkeys():
		allsort += [(res[k]['pvalue'],k)]
	allsort.sort()
	
	# save genes from top 10 categories
	nogogenedir = savedir+'genes4nonsigTerms/'
	createdir(nogogenedir)
	for p,k in allsort[:10]:
		printList(res[k]['genes'],file=nogogenedir+'%s_genes.txt'%(k))
		
	order = ['pvalue', 'count', 'group_count', 'proportion', 'gw-proportion', 'table']
	print >> gof, '\n\n\n'
	print >>gof, 'All GO terms:'
	for p,k in allsort: 
		lab = ''
		if k in akago: lab = akago[k]['name']
		print >> gof, '%s (%s): '%(k,lab),
		ct = 0
		for kk in order:
			if kk != 'genes':
				if ct == 0: gof.write('%s=%s'%(kk,res[k][kk]))
				else: gof.write(', %s=%s'%(kk,res[k][kk]))
				ct += 1
		print >> gof
	return



# USER PREFS
# -------------------------------------------------------------------
indirfile = ''
readfunc = lambda x: readList(x,delim='\r')
bgfile = None

outdir = ''

alpha = .01
method = 'BH'

godir = ''
gofile = godir+'cflo.v3.3.UniProt.GO.txt'
gotable = godir+'gene_ontology_ext.obo.txt'

delim = '\n'

NKEEP = -1

splitnames = False # remove --CRAP suffix to gene names in GO table
verbose = False	

# CLI ARGS
# -------------------------------------------------------------------
help = """\nUsage: GOenrichment.py -i <indir>|gene_lst.txt [-o <outdir> -alpha P-value -method BH|Bonferroni -go gofile.txt -golookup table.txt]
"""
nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

if len(args)==1: sys.exit('\n'.join(help.split('\n')))

while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('in|^i$|^f$|^fg$', arg): indirfile = val; nhelps += 1
	elif re.match('^background$|^b$|^bg$', arg): bgfile = val
	elif re.match('out|^o$', arg): outdir = slash(val)
	elif arg == 'alpha': alpha = float(val)
	elif arg == 'method': method = val
	elif arg == 'go': gofile = val
	elif arg == 'golookup': gotable = val
	elif arg == 'godir': godir = slash(val)
	elif arg == 'table': readfunc = lambda x: readTable(x,header=0)[1]; ai-=1
	elif arg == 'split': splitnames = True; ai-=1
	elif arg == 'verbose': verbose = True; ai-=1
	elif arg == 'keep': NKEEP = int(val)
	elif re.match('^help|h$', arg.lower()): sys.exit(help)
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2

if method != 'BH' and method != 'Bonferroni': sys.exit('Error: method must be one of {BH, Bonferroni}.')


# if godir:
# 	gofile = godir+'cflo.v3.3.UniProt.GO.txt'
# 	gotable = godir+'gene_ontology_ext.obo.txt'



# I/O WORK
# -------------------------------------------------------------------
createdir(outdir)

# DO WORK
# -------------------------------------------------------------------
for f in getFiles(indirfile):
	
	name = f.split('/')[-1]
	asdf = 'cflo.v3.3.UniProt.GO.txt'
	if bgfile: asdf = getLabel(bgfile)
	sys.stdout.write('- GO enrichment, fg=%s, bg=%s...'%(name,asdf)); sys.stdout.flush()
	bginfo = []
	
	keygo = []
	import utilities as ut
	dat = readfunc(f)
	fginfo = dat[:NKEEP]
	
	if bgfile: bginfo = readfunc(bgfile)
	elif NKEEP > 0: bginfo = dat[NKEEP:]
	# res = fullGOanalysis(readfunc(f), bginfo, alpha=alpha, method=method, name=name, savedir=outdir, godir=godir, gofile=gofile, gotable=gotable, splitnames=splitnames, verbose=verbose)
	
	# subtract out fg genes
	bginfo = ut.setDifference(bginfo,fginfo)
	
	res = ut.fullGOanalysis(fginfo, bginfo, alpha=alpha, method=method, name=name, savedir=outdir, godir=godir, gofile=gofile, gotable=gotable, gononsig=100, keywords=keygo)
	
	sys.stdout.write('done\n'); sys.stdout.flush()
















