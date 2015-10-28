#! /usr/bin/env python

import re, sys, os, stat, glob, random                    # python
from glimmrAccessories import *

version = """\ntransformGlimmrMatrices.py, version %s
 
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

integrator = lambda x: subna(*x)
operator = lambda x: len(x)==0 and nastr or integrator(x)
transform = lambda x: x

def transform2(x):
 	y = transform(x)
	return y != nastr and y or 0

na2naut = lambda x: x != nastr and x or 0

operatorstr = ''
transformstr = ''

# USER PREFS
# -------------------------------------------------------------------
file1 = ''
file2 = ''
indirfile = ''
outdir = 'subtractMatrices/'
outfile = None
PLOTLIST = None
PLOT = False
Fcov = 1 # coverage bias correction

upto = 0

# CLI ARGS
# -------------------------------------------------------------------
help = """
Usage: transformGlimmrMatrices.py ...

Given 2 glimmr matrix files, return the difference between them"""

nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('^in$|^i$', arg): file1 = val
	elif arg == 'indir': indirfile = val
	elif arg == 'j': file2 = val
	elif re.match('out|^o$', arg): outdir = slash(val)
	elif arg == 'fn': outfile = val
	elif arg == 'summary': PLOTLIST = val
	elif arg == 'plot': PLOT = True; ai-=1
	elif arg == 'op': operatorstr = val
	elif arg == 'transform': transformstr = val
	elif arg == 'cov': Fcov = float(val)
	elif arg == 'upto': upto = int(val)
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

if PLOT: import plot

# DO WORK
# -------------------------------------------------------------------
if operatorstr == 'subtract': integrator = lambda x: subna(*null2na(x))
elif operatorstr == 'divide': integrator = lambda x: divna(*x)
elif operatorstr == 'multiply': integrator = lambda x: multna(*x)
elif operatorstr == 'add': integrator = lambda x: addna(*x)
elif operatorstr == 'average': integrator = avg

if transformstr == 'log': transform = log
if transformstr == 'lg': transform = lg

if PLOTLIST or PLOT:
	of = None
	L = []
	if outfile:
		of = getLabel(outfile)
		odof = outdir+'%s_sampled_differences.txt'%(of)
	
	# first check if this is the sampled diffs file
	if not outfile and 'sampled_differences' in PLOTLIST:
		print >> sys.stderr, 'loading data', PLOTLIST
		if not upto:
			L = readList(PLOTLIST, dtype=floatna)
		else:
			L = readListPartial(PLOTLIST, dtype=floatna, upto=upto)
		
		of = getLabel(PLOTLIST)
	
	elif os.access(odof, os.F_OK):
		print >> sys.stderr, 'loading data', odof
		if not upto:
			L = readList(odof, dtype=floatna)
		else:
			L = readListPartial(odof, dtype=floatna, upto=upto)
	
	else:
		fh = open(PLOTLIST)
		for line in fh:
			row = line[:-1].split('\t')
			scaf = row.pop(0)
			# print 'working', scaf
			L += map(floatna, random.sample(row,int(.1*len(row))))
		fh.close()
		
		printList(L,file=odof)
		
	if PLOT:
		print 'Plotting...'
		legend = '.95=%s, .99=%s, .999=%s'%(percentile(L,.95),percentile(L,.99),percentile(L,.999))
		tL = random.sample(L,int(.25*len(L)))
		plot.hist(tL, xlabel='difference',file=outdir+'%s_sampled_diff_hist.pdf'%(of), yfreq=1, logscale='y', bins=makeRange(-6,6,n=100), legend=legend)
		
		aL = map(absna,L)
		legend = '.9=%s, .95=%s, .99=%s, .999=%s'%(percentile(aL,.9),percentile(aL,.95),percentile(aL,.99),percentile(aL,.999))
		tL = random.sample(aL,int(.25*len(aL)))
		
		plot.hist(tL, xlabel='abs difference',file=outdir+'%s_sampled_absdiff_hist.pdf'%(of), yfreq=1, bins=makeRange(0,6.1,interval=.1), legend=legend)
		
		sys.stdout.write('%s:%s\n'%(of,legend))
	else:
		aL = map(absna,L)
		legend = '.9=%s, .95=%s, .99=%s, .999=%s'%(percentile(aL,.9),percentile(aL,.95),percentile(aL,.99),percentile(aL,.999))
		sys.stdout.write('%s:%s\n'%(of,legend))
		
	sys.exit()



fo = sys.stdout
createdir(outdir)
if outfile: fo = open(outdir+outfile, 'w')


files = [file1,file2]
if len(indirfile):
	if isdir(indirfile) or os.access(indirfile,os.F_OK):
		files = getFiles(indirfile)
	else:
		files = glob.glob(indirfile)

print '\n%s INPUT FILES\n'%(len(files)), files
print
print 'OPERATOR:', operatorstr
print

fh = {}
for f in files: fh[getLabel(f)] = open(f)
sk = fh.keys()

L = []
nscafs = 0
try:
	while 1:
		scafs = []
		line = [[] for k in sk]
		plen = [0 for k in sk]
		sys.stderr.write('Processing row...'); sys.stderr.flush()
		ct = 0
		FLAG = 0
		for k in sk:
			info = fh[k].next()[:-1].split('\t')
			scaf = info.pop(0)
			# profile = map(float,info[1:])
			plen[ct] = len(info)
			line[ct] = info
			
			# print 'test', scaf, info
			
			# if scaf == 'scaffold1':
			# 	print k, len(map(ut.floatna, info))
			
			scafs += [scaf]
			ct += 1
		
		# if scaf == 'scaffold1': sys.exit()
		
		sys.stderr.write('%s (n=%s bp)'%(scafs[0], len(line[0]))); sys.stderr.flush()
			
		# ensure 
		if variance(plen) > 0.0:
			print >> sys.stderr, 'Error: different scaffold lengths', zip(scafs,plen)
			fo.write(scafs[0]+'\n')
		else:
			# this doesn't work for large genomes
			# linet = ut.transpose(line)
			# D = map(lambda x: str( transform(ut.multna(Fcov, operator(*x))) ), linet)
			
			if Fcov != 1:
				D = [str( na2naut(multna(Fcov, operator([ln[i] for ln in line]))) ) for i in range(len(line[0]))]
			else:
				# if scaf == 'scaffold101':
				# 	for i in range(len(line[0])):
				# 		print i, [ln[i] for ln in line], operator([ln[i] for ln in line])
				
				D = [str( na2naut(operator([ln[i] for ln in line])) ) for i in range(len(line[0]))]
				# print 'LEN D', len(D)
				# if operatorstr == 'average':
				# 	D = [str( na2naut(operator([ln[i] for ln in line])) ) for i in range(len(line[0]))]
				# else:
				# 	# try:
				# 	D = [str( na2naut(operator([ln[i] for ln in line])) ) for i in range(len(line[0]))]
				# 	# except ValueError:
				# 	# 	print 'Value', [[ln[i] for ln in line] for i in range(len(line[0]))]
				# 	# 	break
			
			fo.write(scafs[0]+'\t'+'\t'.join(D)+'\n')
			
			# sample differences to plot distribution
			L += map(floatna, random.sample(D,int(.1*len(D))))
			
		nscafs += 1
		sys.stderr.write('done row.\n'); sys.stderr.flush()
		
except StopIteration:
	print >> sys.stderr, 'DONE'

printList(L,file=outdir+'%s_sampled_differences.txt'%(getLabel(outfile)))

if outfile:
	fo.close()
