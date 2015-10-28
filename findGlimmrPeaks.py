#! /usr/bin/env python

import re, sys, os, stat, math, glob, random
from glimmrAccessories import *

version = """\nfindGlimmrPeaks.py, version %s
 
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

PLOTINCLUDE = False
# import plot

# Important note: this only works up to a certain stringency Q ~ 200, beyond which anything will be considered significant


# LOCAL METHODS
# ------------------------------------------------------------------
def buildRegions2(X,Y, cutoff, gap=0, minlength=None, maxlength=None):
	peaks = zip(X,Y)
	nallpeaks = len(peaks)
	peaks = filter(lambda x: x[1]>=cutoff, peaks) # count peaks above some threshold
	npeaks = len(peaks)
	
	# cluster segments
	clusters = []
	i = 0
	curr = []
	
	# SAW_THE_VOID = False
	# INSIDE_PEAK2 = False
	
	while i < len(peaks):
		# if len(curr): 
		# 	if abs(peaks[i][0] - curr[-1]) <= 10: SAW_THE_VOID = False
		# 	else: SAW_THE_VOID = True
		# else: SAW_THE_VOID = False
		
		# if len(curr) and SAW_THE_VOID:
			# if abs(peaks[i][0] - curr[-1]) <= 10: INSIDE_PEAK2 = True
			# else: INSIDE_PEAK2 = False
			
		
		# print 'i', i, 'lasti', lasti
		if len(curr) and abs(peaks[i][0] - curr[-1]) > gap or (SAW_THE_VOID and INSIDE_PEAK2):
			# print 'Completed cluster', len(curr)
			clusters += [[curr[0], curr[-1]]]
			# print 'DONE CLUSTER', i, lasti, curr[-1]
			curr = []
			# SAW_THE_VOID = False
		else: 
			curr += [peaks[i][0]]
			# SAW_THE_VOID = False
		i += 1
	
	if minlength: clusters = filter(lambda x: abs(x[0]-x[1])+1 >= minlength, clusters)
	if maxlength: clusters = filter(lambda x: abs(x[0]-x[1])+1 < maxlength, clusters)
	
	return {'nallpeaks':nallpeaks, 'npeaks':npeaks, 'nclusters':len(clusters), 'clusters':clusters}



def buildRegions(X, fragmentLength, cutoff=None, direction=lambda x,y: x >= y, gap=0, minlength=None, maxlength=None, loffset=0, hoffset=0, scaffoldlength=100e9):
	nallpeaks = len(X)
	
	# matrix method
	# X = zip(pos,Q,[Qthresh for x in pos],Q)
	# original global fixed cutoff
	# peaks = filter(lambda x: x[1]>=cutoff, X) # count peaks above some threshold
	
	# adaptive threshold
	# X = [pos,Q,cutoff,Q]
	# X[1] is current score, X[2] is cutoff for this position
	peaks = filter(lambda x: direction(x[1],x[2]), X) # count peaks at or above some threshold
	npeaks = len(peaks)
	
	# Q score dictionary
	Qdct = dict(map(lambda x: (x[0],x[3]), peaks))
	
	# filter P scores for profile
	peaks = map(lambda x: (x[0],x[1]), peaks)
	
	# print 'before/after', len(X), len(peaks)
	clusters = []
	i = 0
	
	SAW_PEAK1 = False
	SAW_PEAK2 = False
	curr = []
	
	while i < npeaks:
		
		if not SAW_PEAK1:
			# start a new bound region
			curr = [peaks[i][0]]
			SAW_PEAK1 = True
			# print i, peaks[i], 'new region'
			i += 1
			
		elif SAW_PEAK1 and SAW_PEAK2 and (abs(peaks[i][0] - curr[-1]) > gap):
			# close the region
			# print i, peaks[i], curr[0], curr[-1], abs(peaks[i][0] - curr[-1]), 'vs', gap, '...ending peak', len(clusters)
			# bed format is 0-offset for the start and 1-offset for stop
			clusters += [[max(0,curr[0]-loffset),min(curr[-1]+1+hoffset, scaffoldlength)]]
			
			# reset
			curr = []
			SAW_PEAK1 = False
			SAW_PEAK2 = False
			# print 'i', i
			
		#
		elif SAW_PEAK1 and not SAW_PEAK2 and abs(peaks[i][0] - curr[-1]) > gap:
			# print i, peaks[i][0], 'state', SAW_PEAK1, SAW_PEAK2, abs(peaks[i][0] - curr[-1]), curr
			
			# bed format is 0-offset for the start and 1-offset for stop
			# clusters += [[curr[0]-loffset,curr[-1]+1+hoffset]]
			clusters += [[max(0,curr[0]-loffset),min(curr[-1]+1+hoffset, scaffoldlength)]]
			# reset
			curr = []
			SAW_PEAK1 = False
			SAW_PEAK2 = False
			i+=1
			
		elif SAW_PEAK1 and abs(peaks[i][0] - curr[-1]) <= gap:
			# criteria for peak 2
			# end the peak if current position is a jump from previous
			# or if the existing cluster is long enough already
			if abs(peaks[i][0] - curr[-1]) > 1 or abs(peaks[i][0] - curr[0]) >= fragmentlength: 
				SAW_PEAK2 = True
				
				# if curr[-1]-curr[0] > 600: print 'saw peak2', curr[-1]-curr[0]
			
			curr += [peaks[i][0]]
			# print i, peaks[i], 'within the gap', peaks[i][0], curr[-1], abs(peaks[i][0] - curr[-1]), gap
			i += 1
			
		
	
	# print '- n=%s clusters before size filtering'%(len(clusters))
	# print 'clusters', map(lambda x: abs(x[0]-x[1])+1, clusters)
	
	if minlength: clusters = filter(lambda x: abs(x[0]-x[1]) >= minlength, clusters)
	if maxlength: clusters = filter(lambda x: abs(x[0]-x[1]) < maxlength, clusters)
	
	# Find peak centers
	# -----------------
	# trim back the binding boundaries by finding the peak signal
	# clusters = map(lambda x: int(round(x[0]+(abs(x[0]-x[1])+1)/2.)), clusters)
	# clusters = map(lambda x: [x-int(regionlength/2.), x+1+int(regionlength/2.)], clusters)
	
	# or, find center based on the model position
	pdid = dict(peaks)
	cl2 = []
	for clo,chi in clusters:
		# get binding probability profile over this peak range
		profile = map(lambda x: x in pdid and pdid[x] or 0, range(clo,chi))
		# print 'profile',clo,chi
		# # compute the overall peak score as the average binding P-value
		# score = ut.avg(profile)
		
		# compute overall peak score as average Q-value
		score = avg(map(lambda x: x in Qdct and Qdct[x] or nastr, range(clo,chi)))
		
		# kernel smooth this peak
		profile = kernelsmooth(profile,10,1)
		# print 'pre profile', len(profile), profile[:10]
		# Robustly find the center of the peak
		
		# get top most part of the peak profile
		cut = percentile(profile, .9) # this is the upper 90% cutoff
		profile = filter(lambda x: x[1]>=cut, zip(xrange(clo,chi),profile))
		# profile now contains only positions and scores above cutoff
		# print 'post filter', len(profile), profile[:10]
		# profile = zip(xrange(clo,chi),profile)
		if len(profile): 
			posx, profile = unzip(profile)
			
			# 2 alternatives
			# -------------------------------
			# take the center of this upper 90%
			idx = int( round((posx[-1] - posx[0])/2.) )
		
			# find the first position with greatest score
			# idx = int(round(ut.argmax(profile)))
			# -------------------------------
		
		
			# rprofile = [profile[len(profile)-i-1] for i in range(len(profile))]
			# idx2 = int(round((ut.argmax(profile) + ut.argmax(rprofile))/2.))
		
			# average of the mode and median
			# midx = (posx[1]-posx[0])/2.
			# idx = int(round((.6*idx + .3*midx + .1*idx2)))
		
			# custom = 'set arrow nohead from %s,0 to %s,1 lt 1'%(idx,idx)
			# custom += '; set arrow nohead from %s,0 to %s,1 lt 2'%(idx2,idx2)
			# plot.scatter(zip(range(len(profile)), profile), style='linespoints', custom=custom)
		
			center = posx[0]+idx
			# center = (posx[0] + idx+1 + int(regionlength/2.)) - (posx[0] + idx - int(regionlength/2.))
		
			# cl2 += [[posx[0] + idx - int(regionlength/2.), posx[0] + idx+1 + int(regionlength/2.), posx[0]+idx]]
			# cl2 += [[posx[0] + idx - int(regionlength/2.), posx[0] + idx+1 + int(regionlength/2.), posx[0]+idx]]
		
			cl2 += [[clo,chi,center,score]]
		
	clusters = cl2
	
	return {'nallpeaks':nallpeaks, 'npeaks':npeaks, 'nclusters':len(clusters), 'clusters':clusters}


P2Q = lambda p: p>0 and -10*math.log(p,10) or 0
Q2P = lambda e: (e == 0 and 1-1e-16) or 10**(-e/10.) # avoid 1

logscore = lambda x: lg(addna(x,1))

cfunc = floatna

# USER PREFS
# -------------------------------------------------------------------
method = 'glimmr'
# key = 'ptm' # ptm or nuc
AU = 'ALL'
indirfile = ''
outdir = ''
reffile = ''

goldfile = ''
chrom = ''

PLOTDISTANCEONLY = False
name = ''
SUPERMETA = False

# ---- plot parameters ----
PLOT = False
SPATIALPROFILE = False
SPATIALPLOT = False

diameter = 50 # window size for smoothing
# -------------------------

# ------ PARAMETERS ------
Qthresh = 20.0 # minimum model confidence to begin peak call
dirstr = '';#'greater'
CUTDIRECTION = None; lambda x,y: x >= y

# regionlength = 147

# values based on bioanalyzer results
fragmentlength = 147
fragmentlengthsd = 30
# ------------------------

MINDEPTH = 0

bins = [0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000]

# glimmr starts reporting binding signal about a read length before the true binding start
# so correct for this
loffset = 36 # set this to about the read length
hoffset = 0

# single peak model
MAXGAP = 1 # how many sub-significant positions to allow before breaking the peak?

# what is minimum size of a peak?
# this is determined by resolution offered by fragmentlength and stringency (since this trims in the edges)
# and whether TF or histone

# the minimum is probably approximately the smallest fragment length given 2*sd
# MINLENGTH = fragmentlength - 2*fragmentlengthsd

# MINLENGTH = regionlength/2. # must be slightly bigger for PTM at high redundancy
# MINLENGTH = (2/3.)*regionlength # how does this look?

# but also include fragsize parameters

# MINLENGTH = fragmentlength - 2*fragmentlengthsd # this is fine for nuc only too stingent for PTM

# good for histones
# MINLENGTH = fragmentlength - fragmentlengthsd # this is about right to prevent PTM FP on highly degenerate sequence

# for TFs

# this is now universally good
MINLENGTH = fragmentlength - 2*fragmentlengthsd

MAXLENGTH = None # no cap on maximum size

DISTANCEHISTOGRAM = False
VERBOSE = False
CLEARFILES = 0

# CLI ARGS
# -------------------------------------------------------------------
help = """\nHELP for findGlimmrPeaks.py

Identify significant regions of binding using glimmr.py posterior files.
"""
nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('in|^i$', arg): indirfile = val; nhelps += 1
	elif re.match('out|^o$', arg): outdir = ut.slash(val)
	elif arg == 'ref': reffile = val
	elif arg == 'scaffold': chrom = val
	elif arg == 'gold': goldfile = val
	elif arg == 'method': method = val
	elif re.match('^threshold$|^Q$', arg): Qthresh = float(val)
	elif arg == 'le': CUTDIRECTION = lambda x,y: x <= y; dirstr = 'lesser.'; ai-=1
	elif arg == 'ge': CUTDIRECTION = lambda x,y: x >= y; dirstr = 'greater.'; ai-=1
	elif arg == 'clear': CLEARFILES = True; ai-=1
	elif re.match('^diameter$|^f$', arg): fragmentlength = float(val)
	elif re.match('^diametersd$|^diamsd$', arg): fragmentlengthsd = float(val)
	elif arg == 'minlength' or arg == 'minlen': MINLENGTH = int(val)
	elif re.match('^readlength$|^readlen$', arg): loffset = int(val)
	elif re.match('^maxgap$|^gap$', arg): MAXGAP = int(val)
	elif arg == 'logscore': cfunc = logscore; ai-=1
	elif re.match('verbose|^v$', arg): VERBOSE = True; ai -=1
	elif arg == 'spatialplot': SPATIALPROFILE = True; ai-=1
	elif arg == 'plotdist': PLOTDISTANCEONLY = val
	elif arg == 'supermeta': SUPERMETA = True; ai-=1
	elif arg == 'name': name = val
	elif arg == 'plot': PLOT = True; ai-=1
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)

# I/O WORK
# -------------------------------------------------------------------
createdir(outdir)
sizedir = outdir+'size distributions/'
createdir(sizedir)

# determine name from existing file
if not name:
	name = '.'.join(indirfile.split('/')[-1].split('.')[:-1])
# if VERBOSE: print 'Name:', name


# DO WORK
# -------------------------------------------------------------------

if PLOTDISTANCEONLY and SUPERMETA:
	print 'SUPERMETA!!!'
	dirfiles = glob.glob(PLOTDISTANCEONLY)
	# print 'recovered', dirfiles
	
	try:
		labels = []
		metasizes = []
		metaLS = []
		for f in dirfiles:
			lab = f.split('/')[-1].split('.')[0]
			dirstr = None
			if 'greater' in f:
				lab += '.greater'
				dirstr = 'greater.'
			elif 'lesser' in f:
				lab += '.lesser'
				dirstr = 'lesser.'
			
			print '- Lab', lab
			labels += [lab]
			
			tab,r,c = readTable(f)
			sizes = map(lambda x: abs(int(x[0])-int(x[1])), tab)
			# if not len(sizes): metasizes += [[0]]
			metasizes += sizes
			
			# also make a scatterplot of score vs length
			LS = map(lambda x: (abs(int(x[0])-int(x[1])), cfunc(x[3])), tab)
			metaLS += LS
			
		fo = sizedir+'supermeta peak size distribution freq.pdf'
		
		if PLOTINCLUDE:
			plot.hist(metasizes, bins=[100,125,150,175,200,225,250,275,300], file=fo, custom='set size ratio 1; set yrange [0:*]', xlabel='Peak size (nt)', ylabel='Frequency', yfreq=1)
			
			plot.scatter(metaLS, xlabel='ROI length (nt)', ylabel='ROI score', file=sizedir+'supermeta_score_vs_size.pdf', logscale='x', custom='set grid; set xrange [%s:*]'%(MINLENGTH))
		
	except IOError: sys.exit('Cannot access bedfile %s'%(bedfile))
	
	sys.exit()



if PLOTDISTANCEONLY:
	print 'PLOT DISTANCE ONLY'
	
	# dirfiles = [PLOTDISTANCEONLY]
	# if not os.access(PLOTDISTANCEONLY, os.F_OK):
	# 	dirfiles = glob.glob(PLOTDISTANCEONLY)
	#
	#
	if isdir(PLOTDISTANCEONLY) or os.access(PLOTDISTANCEONLY,os.F_OK):
		dirfiles = getFiles(PLOTDISTANCEONLY)
	else:
		dirfiles = glob.glob(PLOTDISTANCEONLY)
	
	
	print 'recovered', dirfiles
	try:
		labels = []
		metasizes = []
		metaLS = []
		for f in dirfiles:
			# lab = f.split('/')[-1].split('.')[0]
			lab = getLabel(f)
			
			dirstr = None
			if 'greater' in f:
				lab += '.greater'
				dirstr = 'greater.'
			elif 'lesser' in f:
				lab += '.lesser'
				dirstr = 'lesser.'
			
			print '- Lab', lab
			labels += [lab]
			
			tab,r,c = readTable(f)
			sizes = map(lambda x: abs(int(x[0])-int(x[1])), tab)
			# if not len(sizes): metasizes += [[0]]
			metasizes += [sizes]
			
			hist = histogram(sizes, bins=bins)
			hist += [('Total', len(sizes))]
			printTable(hist, file=sizedir+'%s peak size distribution Q%s.%stxt'%(name,Qthresh,dirstr))
			
			
			# also make a scatterplot of score vs length
			LS = map(lambda x: (abs(int(x[0])-int(x[1])), cfunc(x[3])), tab)
			metaLS += [LS]
			if PLOTINCLUDE:
				plot.scatter(LS, xlabel='ROI length (nt)', ylabel='ROI score', file=sizedir+'%s_size_vs_score.%spdf'%(name,dirstr), logscale='x', custom='set grid; set xrange [%s:*]'%(MINLENGTH))
			
		labels = map(lambda x: x.replace('_', '-'), labels)
		
		if PLOTINCLUDE:
			fo = sizedir+'%s peak size distribution freq Q%s.pdf'%(name, Qthresh)
			plot.hist(metasizes, bins=[50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=fo, custom='set size ratio .2; set yrange [0:.4]', xlabel='Peak size (nt)', ylabel='Frequency', yfreq=1, legends=labels, style='histogram')
			
			plot.scatter(metaLS, xlabel='ROI length (nt)', ylabel='ROI relative score difference', file=sizedir+'meta_%s_size_vs_score.pdf'%(name), logscale='x', custom='set grid; set xrange [%s:*]; set key top'%(MINLENGTH), legends=labels)
		
		# count not frequency
		# plot.hist(sizes, bins=[0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=sizedir+'%s peak size distribution Q%s.pdf'%(name, Qthresh), custom='set size ratio .3; set yrange [0:*]', xlabel='Peak size (nt)', ylabel='Count', yfreq=0)#'; set yrange [0:%s]'%(rnge+.05*rnge))
	
	
	except IOError: sys.exit('Cannot access bedfile %s'%(bedfile))
	
	# plot.hist(sizes, bins=[0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=sizedir+'%s peak size distribution Q%s.pdf'%(name, Qthresh), custom='set size ratio .3', xlabel='Peak size (nt)', ylabel='Count')#'; set yrange [0:%s]'%(rnge+.05*rnge))
	
	sys.exit('DONE REPLOT')
############################################


if method == 'glimmr':
	Pthresh = 1-(10**(-float(Qthresh)/10.))
	print 'Name', name, 'Q/thresh', Qthresh, Qthresh

# ---------------------------------------
G = {}
chromlen = {}
nucprofile = {}
ptmprofile = {}
if SPATIALPROFILE:
	print 'Loading reference genome'
	G = readFasta(reffile)
	if chrom:
		for k in G.keys():
			if k != chrom: del G[k]
	print 'DONE'
	
	chromlen = dict( map(lambda x: (x,len(G[x]['seq'])), G.keys()) )
	
	print 'working on spatial profile...'
	for chrom in G.keys():
		# chromlen[chrom] = len(G[chrom]['seq'])
		# print chrom, chromlen[chrom]
		nucprofile[chrom] = [(i,0) for i in range(chromlen[chrom])]
		ptmprofile[chrom] = [(i,0) for i in range(chromlen[chrom])]
	print 'blank loaded into memory.'
# ---------------------------------------

# open pipe for output
# bedfile = None
# bedfile = outdir+'%s peaks Q%s.txt'%(name,Qthresh)
# bedfile = outdir+'%s peaks Q%s.%s.txt'%(name,Qthresh,dirstr)
bedfile = outdir+'%s.ROI.%sbed'%(name,dirstr)
if VERBOSE: print 'Output:', bedfile
if CLEARFILES: os.system('rm -f "%s"'%(bedfile))

pf = open(bedfile,'w')
print >> pf, '# CUTOFF=%s, DIRECTION=%s'%(Qthresh,dirstr)
print >> pf, '# '+'\t'.join(['Chrom', 'Start', 'Stop', 'Center', 'Score', 'Direction'])

# variables for peak finding
nkeep = 0; nlocikeep = 0
nscaffolds = 0 # number of chromosomes with at least 1 peak
npeaks = 0 # number of individual peaks
scaffolds = {} # record number of lines corresponding to each scaffold

# ---------------------------
if method == 'glimmr':
	# Load data from posteriors file
	# this should be done by chromosome to save memory
	# use the same idea as in glimmrToMatrix.py
	print '\nLoading input data...', indirfile
	
	currScaf = None
	profile = [] # store the binding profile for current scaffold
	scafcount = 1
	# variables for peaks themselves
	clusters = {} # temporary storage for bound regions
	
	# process the input data
	fh = open(indirfile)
	for line in fh:
		# if nkeep % 1000000 == 0: sys.stdout.write(',%s'%(ct)); sys.stdout.flush()
		# if nkeep >= 500000: break
		
		# just add info to current scaffold until we reach the next scaffold
		row = line[:-1].split('\t')
		
		# include degeneracy for multiple test correction
		chrom, pos, ref, issig, Q, D, r1, r2 = row[:8]
		info = dict(map(lambda x: x.split(':'), row[8:]))
		dmean = []
		for k in info:
			# print 'test', info[k].split(',')
			lik,dep,mul,tmul,d = info[k].split(',')
			try: dmean += [float(d)] # harmonic
			# try: dmean += [float(tmul)/float(dep)] # arithmetic
			except ValueError:pass
			
		udalpha = avg(dmean)
		pos = int(pos)
		D = float(D)
		aQthresh = 1 - ((1-Pthresh) / (2*fragmentlength * (1 + udalpha)))
		
		P = 1-(10**(-float(Q)/10.))
		
		try: scaffolds[chrom] += 1
		except KeyError: scaffolds[chrom] = 1
		nkeep += 1
		
		if not currScaf: currScaf = chrom; currPileup = {}
		
		# when we have processed all lines for current scaffold/chromosome
		if currScaf and chrom != currScaf:
			scafcount += 1
			# now we find peaks for this scaffold
			if VERBOSE: print '- lines=%s | %s (%s) | filtering (Q>%s)'%(nkeep, currScaf, scafcount,Qthresh),
			# print '- Finding multi-locus bound regions...'
			res = buildRegions(profile, fragmentlength, gap=MAXGAP, minlength=MINLENGTH, maxlength=MAXLENGTH, loffset=loffset, hoffset=hoffset)
			# clusters[chrom] = res['clusters']
			if VERBOSE: print '=> %s/%s => %s peaks'%(res['npeaks'], res['nallpeaks'], res['nclusters'])
			
			nlocikeep += res['npeaks']
			if res['nclusters'] > 0: nscaffolds += 1
			
			# process each peak
			for x in res['clusters']:
				#                                     start, stop, center, Q-score
				print >> pf, '\t'.join(map(str,[cuffScaf, x[0], x[1], x[2], round(x[3],3)]))
				npeaks += 1
			# now we still need to save the current line's data
			# which corresponds to a new scaffold
			currScaf = chrom
			profile = [(pos,P,aQthresh,Q)] # reset information for new scaffold
		
		if D >= MINDEPTH: profile += [(pos,P,aQthresh,Q)]
		
		# construct spatial probability profile
		if SPATIALPROFILE and issig == 'nuc': nucprofile[chrom][pos] = (pos,P)
		elif SPATIALPROFILE and issig == 'ptm': ptmprofile[chrom][pos] = (pos,P)

	# deal with the last chromosome
	scafcount += 1
	# now we find peaks for this scaffold
	if VERBOSE: print '- lines=%s | %s (%s) | filtering (Q>%s)'%(nkeep, currScaf, scafcount,Qthresh),
	# print '- Finding multi-locus bound regions...'
	
	res = buildRegions(profile, fragmentlength, cutoff=Qthresh, gap=MAXGAP, minlength=MINLENGTH, maxlength=MAXLENGTH, loffset=loffset, hoffset=hoffset)
	# clusters[chrom] = res['clusters']
	if VERBOSE: print '=> %s/%s => %s peaks'%(res['npeaks'], res['nallpeaks'], res['nclusters'])

	nlocikeep += res['npeaks']
	if res['nclusters'] > 0: nscaffolds += 1

	# process each peak
	for x in res['clusters']:
		#                                     start, stop, center, score
		print >> pf, '\t'.join(map(str,[chrom, x[0], x[1], x[2], x[3]]))
		npeaks += 1
		
	# now we still need to save the current line's data
	# which corresponds to a new scaffold
	# currScaf = chrom
	# profile = [(pos,P,aQthresh)] # reset information for new scaffold
	
	fh.close()

elif method == 'matrix':
	# Load data from matrix file
	# this should be done by chromosome to save memory
	print '\n [[Matrix]] Loading input data...', indirfile
	currScaf = None; scafcount = 1
	clusters = {} # temporary storage for bound regions
	
	# process the input data
	fh = open(indirfile)
	for line in fh:
		# if nkeep % 1000000 == 0: sys.stdout.write(',%s'%(ct)); sys.stdout.flush()
		# if nkeep >= 500000: break
		
		# just add info to current scaffold until we reach the next scaffold
		row = line[:-1].split('\t')
		
		chrom = row.pop(0)
		pos = xrange(len(row))
		Q = map(floatna,row)
		# cut = [Qthresh for x in pos]
		# profile =
		
		try: scaffolds[chrom] += 1
		except KeyError: scaffolds[chrom] = 1
		nkeep += 1
		
		if not currScaf: currScaf = chrom; currPileup = {}
		
		scafcount += 1
		# now we find peaks for this scaffold
		if VERBOSE: print '- lines=%s | %s (%s) | filtering (Q>%s)'%(nkeep, currScaf, scafcount,Qthresh),
		# print '- Finding multi-locus bound regions...'
		resP = buildRegions(zip(pos,Q,[Qthresh for x in pos],Q), fragmentlength, direction=lambda x,y: x >= y, gap=MAXGAP, minlength=MINLENGTH, maxlength=MAXLENGTH, loffset=loffset, hoffset=hoffset, scaffoldlength=len(pos))
		resN = buildRegions(zip(pos,Q,[-Qthresh for x in pos],Q), fragmentlength, direction=lambda x,y: x < y, gap=MAXGAP, minlength=MINLENGTH, maxlength=MAXLENGTH, loffset=loffset, hoffset=hoffset, scaffoldlength=len(pos))
		
		# merge positive and negative
		tPeaks = resP['npeaks']+resN['npeaks']
		taPeaks = resP['nallpeaks']+resN['nallpeaks']
		tClusters = resP['nclusters']+resN['nclusters']
		if VERBOSE: print '=> %s/%s => %s peaks'%(tPeaks,taPeaks,tClusters)
		
		nlocikeep += tPeaks
		if tClusters > 0: nscaffolds += 1
		
		# process each peak
		for p,strand in [(resP['clusters'],'+'), (resN['clusters'],'-')]:
			for x in p:
				#                                     start, stop, center, Q-score
				print >> pf, '\t'.join(map(str,[chrom, x[0], x[1], x[2], round(x[3],3), strand]))
				npeaks += 1
		
		# construct spatial probability profile
		# if SPATIALPROFILE and issig == 'nuc': nucprofile[chrom][pos] = (pos,P)
		# elif SPATIALPROFILE and issig == 'ptm': ptmprofile[chrom][pos] = (pos,P)
	
	fh.close()
# ---------------------------

# close the bed output file
pf.close()

# DONE PEAK FINDING: Summarize results
# ------------------------------------
print '- Retaining %s/%s loci (post-threshold)'%(nlocikeep,nkeep)
print '- n=%s bound regions (after size filter)'%(npeaks)
print '- Number of scaffolds with at least 1 peak:', nscaffolds

tab,r,c = readTable(bedfile)
sizes = map(lambda x: abs(int(x[0])-int(x[1])), tab)
hist = histogram(sizes, bins=bins)
hist += [('Total', len(sizes))]
printTable(hist, file=sizedir+'%s peak size distribution Q%s.txt'%(name,Qthresh))

tab,r,c = readTable(bedfile)
sizes = map(lambda x: abs(int(x[0])-int(x[1])), tab)
hist = histogram(sizes, bins=bins)
hist += [('Total', len(sizes))]

printTable(hist, file=sizedir+'%s peak size distribution Q%s.%s.txt'%(name,Qthresh,dirstr))

if PLOTINCLUDE and PLOT:
	plot.hist(sizes, bins=[0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=sizedir+'%s peak size distribution freq Q%s.pdf'%(name, Qthresh), custom='set size ratio .3; set yrange [0:.4]', xlabel='Peak size (nt)', ylabel='Frequency', yfreq=1)#'; set yrange [0:%s]'%(rnge+.05*rnge))
	
	plot.hist(sizes, bins=[0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=sizedir+'%s peak size distribution Q%s.pdf'%(name, Qthresh), custom='set size ratio .3; set yrange [0:*]', xlabel='Peak size (nt)', ylabel='Count', yfreq=0)#'; set yrange [0:%s]'%(rnge+.05*rnge))
	
	# rnge = ut.listRange(sizes)
	plot.hist(sizes, bins=bins, file=sizedir+'%s peak size distribution Q%s.png'%(name,Qthresh), custom='set size ratio .3', xlabel='Peak size (nt)', ylabel='Count')#'; set yrange [0:%s]'%(rnge+.05*rnge))
	# plot.hist(sizes, bins=[1,5,10,25,50,75,100,125,150,175,200,225,250,275,300], file=outdir+'Size distro %s.pdf'%(chrom))
# ------------------------------------


# -------------------------------------------- VALIDATION --------------------------------------------
# ----------------------------------------------------------------------------------------------------
# compute distance between predicted and known binding sites

if not goldfile: sys.exit('No true binding file provided. Results stored in\n%s'%(outdir))

# load up true results
mbc = {}
# markedhash = {}
fh = open(goldfile)
curreg = None
currmin = 0
currmax = 0
for line in fh:
	row = line[:-1].split('\t')
	chrom,pos,regionid = row
	pos = int(pos)
	# markedhash[(chrom,pos)] = 1
	# print 'regionid', regionid
	if not curreg: 
		curreg = regionid
		currmin = pos
		currmax = pos
	elif curreg == regionid:
		currmax = pos
	else:
		# print 'adding', curreg
		try: mbc[chrom] += [[chrom,currmin,currmax,curreg]]
		except KeyError: mbc[chrom] = [[chrom,currmin,currmax,curreg]]
	
		curreg = regionid
		currmin = pos
		currmax = pos
# fallen out - print last one
try: mbc[chrom] += [[chrom,currmin,currmax,curreg]]
except KeyError: mbc[chrom] = [[chrom,currmin,currmax,curreg]]


# plot profile of probability values and true bound regions
if SPATIALPROFILE:
	print 'working on spatial profile...'
	lw = None#21000
	hg = None#24000
	
	X = vslice(nucprofile[chrom],0)#[lw:hg]
	# Y = ut.vslice(nucprofile[chrom],1)#[lw:hg]
	# Y = ut.kernelsmooth(Y, diameter)
	# nucprofile = zip(X,Y)

	Y = vslice(ptmprofile[chrom],1)#[lw:hg]
	# Y = ut.kernelsmooth(Y, diameter)
	ptmprofile = zip(X,Y)

	# spatial correlation
	# for chrom in [chrom]:#sorted(mbc.keys()):
	truedat = [(i,0) for i in xrange(chromlen[chrom])]
	# truedat = [(i,0) for i in xrange(35000)]
	for chrom,mn,mx,GID in mbc[chrom]:
		for i in range(mn,mx): truedat[i] = [i,1.2]
	
	if PLOTINCLUDE:
		plot.scatter([truedat,ptmprofile], xlabel='Position', ylabel='Posterior probability', style='lines', file=outdir+'PP profile Q%s.pdf'%(Qthresh), showstats=0, legends=['True', 'Ptm'], lineWeights=[2,2], colors=[1,3], custom='set yrange [0:1.05]; set size ratio .25')
	
		for lw,hg in [(0,1000), (1000,2000), (5000,10000)]:
			plot.scatter([truedat[lw:hg],ptmprofile[lw:hg]], xlabel='Position', ylabel='Posterior probability', style='lines', file=outdir+'PP profile Q%s %s,%s.pdf'%(Qthresh,lw,hg), showstats=0, legends=['True', 'Ptm'], lineWeights=[2,2], colors=[1,3], custom='set yrange [0:1.25]; set size ratio .25')
	
		# plot.scatter([nucprofile,ptmprofile], xlabel='Position', ylabel='Posterior probability', style='lines', file=outdir+'PP profile.pdf', showstats=0, legends=['Nuc', 'Ptm'], lineWeights=[2,2], colors=[1,3], custom='set yrange [0:1.05]; set size ratio .25')
		# sys.exit()

		# plot.hist(metaQ, nbins=25, title='quality distribution')
		# plot.hist(metaD, nbins=25, title='depth distribution')

	# for chrom in Dct.keys():
		# Dct[chrom]['bounds'][1] = Dct[chrom]['bounds'][0]+len(Dct[chrom]['P'])



# compute distances

# print '- Loading data...'
d,r,c = readTable(bedfile, rownames=0)

# if method == 'macs':
# 	d = filter(lambda x: float(x[8])<=FDRCUTOFF*100., d)
# 	print 'there are %s peaks above %s FDR'%(len(d), FDRCUTOFF)


# size distribution
sizes = map(lambda x: abs(int(x[1])-int(x[2])), d)
hist = histogram(sizes, bins=bins)
hist += [('Total', len(sizes))]
# io.printTable(hist, file=sizedir+'%s peak size distribution Q%s.txt'%(name,Qthresh))

if PLOTINCLUDE:
	plot.hist(sizes, bins=[0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=sizedir+'%s peak size distribution freq Q%s.pdf'%(name, Qthresh), custom='set size ratio .3; set yrange [0:.4]', xlabel='Peak size (nt)', ylabel='Frequency', yfreq=1)#'; set yrange [0:%s]'%(rnge+.05*rnge))
	plot.hist(sizes, bins=[0,5,10,25,50,75,100,125,150,175,200,225,250,275,300,350,400,450,500,1000,2000,5000,10000], file=sizedir+'%s peak size distribution Q%s.pdf'%(name, Qthresh), custom='set size ratio .3; set yrange [0:*]', xlabel='Peak size (nt)', ylabel='Count', yfreq=0)#'; set yrange [0:%s]'%(rnge+.05*rnge))

	# rnge = ut.listRange(sizes)
	plot.hist(sizes, bins=bins, file=sizedir+'%s peak size distribution Q%s.png'%(name,Qthresh), custom='set size ratio .3', xlabel='Peak size (nt)', ylabel='Count')#'; set yrange [0:%s]'%(rnge+.05*rnge))
	# plot.hist(sizes, bins=[1,5,10,25,50,75,100,125,150,175,200,225,250,275,300], file=outdir+'Size distro %s.pdf'%(chrom))



NPREDICTIONS = len(d)

if DISTANCEHISTOGRAM:
	print '- Computing pairwise distance distribution between peaks to confirm gap joining'
	ds = []
	
	tups = map(lambda x: (x[0],int(x[1]),int(x[2])), d)
	dist = lambda x,y: min( map(lambda z: abs(z[0]-z[1]), [(x[0],y[0]), (x[0],y[1]), (x[1],y[0]), (x[1],y[1])]) )
	
	
	ct = 0
	for i in range(len(tups)):
		# if random.random() > 0.01: continue
		for j in range(len(tups)):
			# only compare peaks on the same scaffold
			if i != j and tups[i][0] == tups[j][0]:
				if ct % 100000 == 0: sys.stdout.write(',%s'%(ct)); sys.stdout.flush()
				
				ds += [dist(tups[i][1:],tups[j][1:])]
				ct += 1
	sys.stdout.write('Done.\n'); sys.stdout.flush()
	
	print '- Plotting histogram'
	bins = [0,5,10,20,30,40,50,75,100,150,200,250,300,400,500,1000,1500,2000,2500,5000,10000,100000,1000000]
	printList(ds,file=outdir+method+' distances.txt')
	# plot.hist(ds,bins=bins,file=outdir+method+' distance distro.pdf')
	if PLOTINCLUDE: plot.hist(ds,bins=bins,file=outdir+method+' distance distro.pdf', logscale='')

# set up d as a dictionary keyed by scaffold: d is the bedfile contents / the peaks
ddct = dict([(d[i][0],[]) for i in range(len(d))])
for i in range(len(d)): ddct[d[i][0]] += [d[i]]

# clear existing file
os.system('rm -f "%s"'%(outdir+'performance_report_Q%s.txt'%(Qthresh)))
# create a blank file?
# os.system('touch "%s"'%(outdir+'performance_report_Q%s.txt'%(Qthresh)))

if chrom not in ddct: sys.exit("no data. exiting.")

N = 0
M = []
NM = []
minds = []

within = 0
without = 0
total = 0

seen = {}

# iterate over the true binding regions and compare to predictions
for chromy in sorted(mbc.keys()):
	for chromz,mn,mx,GID in mbc[chromy]:
		BOUNDS = [mn,mx]
		TRUECENTER = min(*BOUNDS)+abs(mx-mn)/2.
		
		# all predictions on this chromosome
		ddce = ddct[chromz]
		
		P = map(lambda x: (int(x[1]),int(x[2])), ddce) # start,stop
		
		PREDCENTERS = None
		if method == 'glimmr':
			# improved modal center prediction
			PREDCENTERS = map(lambda x: int(x[3]), ddce)
		else:
			# standard median center prediction
			PREDCENTERS = map(lambda x: min(*x)+abs(x[1]-x[0])/2., P)
		
		# find the prediction that is closest to the current true region
		pairs = dict(map(lambda i: (abs(PREDCENTERS[i]-TRUECENTER),(PREDCENTERS[i],P[i])), xrange(len(PREDCENTERS))))
		mind = minna(map(lambda pred: abs(pred-TRUECENTER), PREDCENTERS))
		if mind == nastr: continue
		
		# pairs = dict(map(lambda x: (min(map(lambda B: abs(B-x), TRUECENTER)),x), PREDCENTERS))
		# mind = ut.minna(map(lambda x: min(map(lambda B: abs(B-x), TRUECENTER)), PREDCENTERS))
		
		# pairs = dict(map(lambda x: (min(map(lambda B: min(abs(B-x[0]),abs(B-x[1])), TRUECENTER)),x), PREDCENTERS))
		# mind = ut.minna(map(lambda x: min(map(lambda B: min(abs(B-x[0]),abs(B-x[1])), TRUECENTER)), PREDCENTERS))
		
		
		# the best distance value and region
		thep = pairs[mind][0]
		thepeak = pairs[mind][1]
		
		if valueWithinRange(thepeak,TRUECENTER):
			within += 1
			minds += [mind] # only count distance for true positives
		else: without += 1
		total += 1
		
		if VERBOSE: print GID, chromz, mn, mx, 'closest prediction', thep, 'dist = %s'%(mind)
		
		N += 1
		
		# compare true with predicted boundary regions - look for overlap
		MATCH = False
		m1 = map(lambda Pi: intervalWithinRange(Pi, BOUNDS), P)
		m1 = filter(lambda x: x!=0 and True, m1)
		
		# m1 = filter(lambda Pi: ut.valueWithinRange(Pi, BOUNDS[0]), P)
		# m2 = filter(lambda Pi: ut.valueWithinRange(Pi, BOUNDS[1]), P)
		
		if len(m1) > 0: MATCH = True
		
		# if len(m1) > 0 or len(m2) > 0: MATCH = True
		if MATCH: M += [GID]
		else: NM += [GID]
		
		# seen[B] = 1

# print 'Tested N=%s, M=%s, NM=%s'%(N,len(M), len(NM))
# print 'Within %s (%s), Without %s (%s), Total %s'%(within, within/float(total), without, without/float(total), total)

thebins = [0,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,175,200,300,400,500,600,700,800,900,1000,2500,5000,10000,20000]
printTable(histogram(minds, bins=thebins), file=outdir+method+' mindist hist.txt')
if PLOTINCLUDE: plot.hist(minds, bins=thebins, file=outdir+method+' mindist Q%s.pdf'%(Qthresh), custom='set size ratio .25; set yrange [*:*]', xlabel='Distance to true binding center')


fh = None
if method == 'glimmr':
	fh = open(outdir+'performance_report_Q%s.txt'%(Qthresh), 'w')
elif method == 'macs':
	fh = open(outdir+'macs_performance_report.txt', 'w')
for od in [sys.stdout, fh]:
	# report sensitivity and specificity
	print >>od, 'Performance report:'
	print >>od, 'Predicted peaks (Q>%s) = %s'%(Qthresh, NPREDICTIONS)
	print >>od, 'True bound regions = %s'%(N)
	
	TP = within
	TN = 0
	FP = max(0, NPREDICTIONS - TP)
	FN = N - TP
	TPR = TP/float(TP + FN)
	FDR = FP/float(FP + TP)
	# FPR = FP/float(NPREDICTIONS)

	print >>od, 'TP = %s, TPR = %s, FN = %s'%(TP, TPR, FN)
	print >>od, 'FP = %s, FDR = %s'%(FP, FDR)
	print >>od, 'Min distance: Average = %s, Median = %s, STD = %s'%(avg(minds),median(minds),stdev(minds))
fh.close()


if SPATIALPROFILE:
	lw = None#21000
	hg = None#24000
	
	# spatial correlation
	for chromy in sorted(mbc.keys()):
		truedat = [(i,0) for i in xrange(chromlen[chromy])]
		
		# truedat = [(i,0) for i in xrange(35000)]
		
		for chromz,mn,mx,GID in mbc[chromy]:
			for i in range(mn,mx): truedat[i] = [i,1]

		predat = [(i,0) for i in xrange(chromlen[chromy])]
		# predat = [(i,0) for i in xrange(35000)]
		for info in ddct[chromy]:
			mn = int(info[1]); mx = int(info[2])
			for i in range(mn,mx): predat[i] = [i,.8]
			
		# true = ut.vslice(truedat,1)
		# pred = ut.vslice(predat,1)
		# print 'Correlation of true vs prediction PTM maps: s=%s'%(ut.cor(true,pred, method='spearman'))
		# print 'Correlation of true vs prediction PTM maps: r=%s, s=%s'%(ut.cor(true,pred, method='pearson'), ut.cor(true,pred, method='spearman'))
		
		if PLOTINCLUDE:
			for lw,hg in [(0,1000), (1000,2000), (5000,10000)]:
				plot.scatter([truedat[lw:hg],predat[lw:hg]], style='lines', file=outdir+'True vs predicted plot Q%s %s,%s,%s.pdf'%(Qthresh,chromy,lw,hg), xlabel='Position', ylabel='Binding', custom='set yrange [0:1.05]; set size ratio .25', lw=2, showstats=0, legends=['True', 'Pred'])
		
			plot.scatter([truedat,predat], style='lines', file=outdir+'True vs predicted plot Q%s %s.pdf'%(Qthresh,chromy), xlabel='Position', ylabel='Binding', custom='set yrange [0:1.05]; set size ratio .25', lw=2, showstats=0, legends=['True', 'Pred'])
		