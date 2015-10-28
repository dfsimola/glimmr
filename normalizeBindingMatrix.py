#! /usr/bin/env python

import re, sys, os, stat, random, math
from glimmrAccessories import *

version = """\nnormalizeBindingMatrix.py, version %s
 
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
cfunc = lambda x: lg(addna(x,1))

# USER PREFS
# -------------------------------------------------------------------
indirfile = ''
outdir = 'debug/'
outfile = ''
SAVE = True
FORCECORRECTION = None
RAWSTATS = True
SDSCALE = False
DEBUG = False
statsfile = None

spikeInValue = None

# CLI ARGS
# -------------------------------------------------------------------
help = '\nHelp for normalizeBindingMatrix.py\nlog transform and mean center matrix of glimmr data.'
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
	elif re.match('out|^o$', arg): outfile = val
	elif arg == 'stats': statsfile = val
	elif arg == 'nosave': SAVE = False; ai-=1
	# elif arg == 'sdscale': SAVE = False; SDSCALE = True; ai-=1
	elif arg == 'force': FORCECORRECTION = map(float,val.split(','))
	elif arg == 'debug': DEBUG = True; outdir = val
	elif arg == 'raw': RAWSTATS = True; ai-=1
	elif arg == 'spikein': spikeInValue = float(val)
	
	# elif arg == 'logscale': RAWSTATS = False; ai-=1
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)


# DO WORK
# -------------------------------------------------------------------

if not len(outfile) and SAVE == True: sys.exit('Error: must specify outfile using -o <outfile.txt>')

f = indirfile
print >> sys.stderr, f.split("/")[-1]
if spikeInValue: print >> sys.stderr, 'Further scaling by spike in value:', spikeInValue

# scan and take average method
asum = 0
allvals = []
pallvals = []
alen = 0
cavg = None
tol = .01
nchange = 0

def process(line,USERAW=True):
	row = line[:-1].split('\t')
	
	def iproc(val):
		global asum, alen, allvals, pallvals
		try: cv = (USERAW and floatna(val)) or cfunc(val)
		except Exception,e: print 'error', val
		# only consider positive values - the 0s are gaps in data
		if cv != nastr and cv > 0:
			asum += cv
			alen += 1
			allvals += [cv]
	
	# def iprocraw(val):
	# 	global asum, alen, allvals, pallvals
	# 	try: cv = ut.floatna(val)
	# 	except Exception,e: print 'error', val
	# 	# only consider positive values - the 0s are gaps in data
	# 	if cv != ut.nastr and cv > 0:
	# 		asum += cv
	# 		alen += 1
	# 		allvals += [cv]
	#
	# if RAWSTATS: map(iprocraw,filter(lambda x: random.random() < .1, row[1:]))
	# else: map(iproc,filter(lambda x: random.random() < .1, row[1:]))
	
	map(iproc,filter(lambda x: random.random() < .1, row[1:]))


def outprocess(line,mean,sd,shift,pipe):
	row = line[:-1].split('\t')
	
	# this is a gain correction for the range of each sample
	# vals = map(lambda x: str(ut.roundna(3)(ut.divna(cfunc(x), val))), row[1:])
	
	vals = []
	for x in row[1:]:
		# log scale
		# v = ut.addna1(shift)(ut.powna(math.e)(ut.divna(ut.subna(ut.ln(cfunc(x)),mean), sd)))
		# v = ut.addna1(shift)(ut.divna(ut.subna(cfunc(x),mean), sd))
		
		# normal scale
		v = addna1(shift)(divna(subna(cfunc(x),mean), sd))
		vals += [str(roundna(3)(v > 0 and v or 0))]
	
	# vals = map(lambda x: str(ut.roundna(3)(ut.addna1(1.)(ut.divna(ut.subna(cfunc(x),mean), sd)))), row[1:])
	
	# no we want to subtract the mean from each value
	# vals = map(lambda x: str(ut.roundna(3)(ut.subna(cfunc(x), val))), row[1:])
	
	pipe.write('\t'.join([row[0]]+vals)+'\n')

#
def outprocessgain(line,gain,pipe):
	row = line[:-1].split('\t')
	
	vals = []
	for x in row[1:]:
		# normal scale
		v = cfunc(divna(x,gain))
		# v = cfunc(ut.divna(ut.subna(x,gain[0]),gain[1]))
		# v = cfunc(ut.divna(ut.subna(x,gain[0]),gain[1]))
		
		# normalize on the log scale
		# v = ut.divna(cfunc(x),gain)
		
		vals += [str(roundna(3)(v > 0 and v or 0))]
	
	pipe.write('\t'.join([row[0]]+vals)+'\n')
	


#
def outprocesssimple(line,pipe):
	row = line[:-1].split('\t')
	
	vals = []
	for x in row[1:]:
		# normal scale
		v = cfunc(x)
		vals += [str(roundna(3)(v > 0 and v or 0))]
	
	pipe.write('\t'.join([row[0]]+vals)+'\n')
	

# -------------------------
# process file data
if not FORCECORRECTION:
	print >> sys.stderr, ' - Loading...'
	fh = open(f)
	
	ct = 1
	for line in fh:
		ct += 1
		process(line,RAWSTATS)
		
		if alen > 0:
			curavg = asum/float(alen)
			try:
				# print line[0:line.index('\t')], curavg, cavg, 'nchange', nchange
				if cavg and abs(curavg-cavg) < tol:
					cavg = curavg
					if nchange == 20000: break
					nchange += 1
				else: nchange = 0
				cavg = curavg
				# print 'cavg', curavg
			except ValueError: pass
	fh.close()

	# print >> sys.stderr, ' - Computing average...'
	# scaffoldnames = ut.vslice(triples,0)
	# metasum = sum(ut.vslice(triples,1))
	# metalen = sum(ut.vslice(triples,2))
	
	# estimate location and scale parameters assuming normal distribution
	
	# stats for normalization
	metaavg = avg(allvals)
	metasd = stdev(allvals)
	metamed = median(allvals)
	meta75 = percentile(allvals,.75)
	
	# estimate location and scale parameters assuming log-normal distribution
	# allvals = filter(lambda x: x>0, allvals)
	# mu = sum([math.log(x+.5,math.e) for x in allvals])/float(len(allvals))
	# metaavg = math.e**mu
	# metasd = (sum([(math.log(x+.5,math.e) - mu)**2 for x in allvals])/float(len(allvals)))**.5
	
else:
	metaavg = FORCECORRECTION[0]
	metasd = FORCECORRECTION[1]
	
# -------------------------

if DEBUG:
	fname = getPath(outfile)+getLabel(outfile)+'.pre.sample.distribution.txt'
	printList(filter(lambda x: random.random() < .05, allvals),file=fname)

# print >> sys.stderr, 'Average=%s, SD=%s (n=%s)'%(metaavg, metasd, alen)
# printTable(ut.histogram(allvals,bins=ut.makeRange(0,11,interval=.5)))
# print >> sys.stderr, ' - Quantiles: .5=%.3f, .75=%.3f, 0.9=%.3f, 0.95=%.3f, 0.99=%.3f, 0.999=%.3f'%(qs[.5],qs[.75],qs[.9],qs[.95],qs[.99],qs[.999])
tagline = getLabel(outfile)
if not SAVE: tagline = getLabel(indirfile)

# save a stats file
fo = sys.stdout
if statsfile: 
	# always save stats on log scale
	if RAWSTATS: allvals = map(cfunc,allvals)
	meany = avg(allvals)
	standy = stdev(allvals)
	qs = quantiles(allvals, Q=[0.05,0.15,0.25,0.5,0.75,0.8,0.85,0.9,0.95,0.99,0.999])
	del allvals
	fo = open(statsfile,'a')
	print >>fo, tagline+'\t'+'Average=%s, SD=%s (n=%s)'%(meany, standy, alen)
	print >>fo, tagline+' Quantiles: .05=%.3f .15=%.3f .25=%.3f .5=%.3f .75=%.3f .8=%.3f .85=%.3f .9=%.3f 0.95=%.3f 0.99=%.3f 0.999=%.3f'%(qs[.05], qs[.15], qs[.25], qs[.5], qs[.75], qs[.8], qs[.85], qs[.9], qs[.95], qs[.99], qs[.999])
	fo.close()

# if poststatsfile:
# 	fo = open(poststatsfile,'a')
# 	print >>fo, tagline+'\t'+'Average=%s, SD=%s (n=%s)'%(metaavg, metasd, alen)
# 	print >>fo, tagline+' Quantiles: .05=%.3f .15=%.3f .25=%.3f .5=%.3f .75=%.3f .8=%.3f .85=%.3f .9=%.3f 0.95=%.3f 0.99=%.3f 0.999=%.3f'%(pqs[.05], pqs[.15], pqs[.25], pqs[.5], pqs[.75], pqs[.8], pqs[.85], pqs[.9], pqs[.95], pqs[.99], pqs[.999])
# 	fo.close()

# save data
if SAVE:
	print >> sys.stderr, ' - Saving...'
	fh = open(f)
	fho = open(outfile,'w')
	
	# gain-correction - equivalent to centering a log-normal distribution
	# assuming equal variance
	# [outprocess(line,0,max(metaavg,1.),0,fho) for line in fh]
	
	
	# gain correction of raw values
	print >> sys.stderr, 'Positive Mean=%s'%(metaavg)
	
	# scale by mean
	# [outprocessgain(line,metaavg,fho) for line in fh]
	
	
	K = meta75
	# if we have spike in control, divide by this value
	if spikeInValue: K = meta75*spikeInValue
	print >> sys.stderr, '%s, K=%s'%(outfile, K)
	
	# scale by 75th percentile (upper quartile)
	[outprocessgain(line,K,fho) for line in fh]
	
	
	# scale by median
	# [outprocessgain(line,metamed,fho) for line in fh]
	
	# [outprocessgain(line,(metasd,metaavg),fho) for line in fh]
	# [outprocess(line,0,max(metaavg,1.),0,fho) for line in fh]
	
	# the mean and sd of the log2+(x+1) transformation is so close to Z that there is no need to manipulate further
	# [outprocesssimple(line,fho) for line in fh]
	
	# Z-transformation, setting the mean to +1
	# [outprocess(line,metaavg,metasd,1.,fho) for line in fh]
	
	
	
	# what could be a better transformation?
	# need to understand the biases of glimmr better, but on first glance
	# there is some kind of gain bias due to crosslinking incubation time
	# so the distribution assumption should be exponential (log-normal on log scale)
	# regress the log(data) against a lognormal distribution and invert
	# this is similar just to a Z transformation, except that Z over-corrects
	
	fho.close()
	fh.close()






