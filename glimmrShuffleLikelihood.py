#! /usr/bin/env python

import re, sys, os, stat, glob, random
from glimmrAccessories import *

version = """\nglimmrShuffleLikelihood.py, version %s
 
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


# USER PREFS
# -------------------------------------------------------------------
indirfile = ''
outdir = ''
DRYRUN = False
refGenomeFile = None
nthreads = 1

# CLI ARGS
# -------------------------------------------------------------------
help = '\nUsage: glimmrShuffleLikelihood.py -i <FILE_OR_DIR> -r <REF_GENOME.FA> [-o <OUTDIR=PWD>] [--dryrun]  \nshuffles chromosome and position information for a provided likelihood file'
nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

# if len(args)==1: sys.exit('\n'.join(help.split('\n')))

if argstr.count('-h') or argstr.count('-help'):
	print >> sys.stderr, help,
	for k in command.keys():
		print >> sys.stderr, '  %s:%s'%(k,command[k])
	print >> sys.stderr
	sys.exit()


while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('in|^i$', arg): indirfile = val; nhelps += 1
	elif re.match('out|^o$', arg): outdir = slash(val)
	elif arg == 'r' or arg == 'ref': refGenomeFile = val
	elif arg == 'x': nthreads = int(val)
	elif arg == 'dryrun': DRYRUN = True; ai-=1
	
	elif re.match('^help|h$', arg.lower()): sys.exit(help)
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)

files = []
if isdir(indirfile) or os.access(indirfile,os.F_OK): files = getFiles(indirfile)
else: files = glob.glob(indirfile)

if not files: sys.exit('No input files.')

createdir(outdir)


if not DRYRUN:
	print 'Loading reference genome...'
	G = readFasta(refGenomeFile)
	allscaffolds = {}
	for k in G.keys(): allscaffolds[k] = len(G[k]['seq'])
	# del G
	asd = allscaffolds.keys()

	# list of scaffold,weight pairs
	# total = float(ut.sumna(allscaffolds.values()))
	# urn = [(s,length/total) for s,length in allscaffolds.items()]


for f in files:
	
	lab = getLabel(f)
	newf = outdir+lab[:lab.index('.mother')]+'.shuffle.mother.unsort.'+getSuffix(f)
	newfsort = outdir+lab[:lab.index('.mother')]+'.shuffle.mother.'+getSuffix(f)
	
	if os.access(newfsort, os.F_OK):
		print >> sys.stderr, 'File exists, skipping: %s'%(newfsort)
		continue
	
	print 'File: %s'%(f)
	
	if DRYRUN: 
		print 'Output: %s'%(newfsort)
		continue
	
	fh = open(f)
	fho = open(newf,'w')
	
	for line in fh:
		row = line[:-1].split('\t')
		try:
			crand = row[0]# random.choice(asd) # fix chrom
			# right way to do it is to weight the sampling by size of each scaffold
			# but this is way too slow
			# crand = ut.weightedChoice(urn)
			
			try:
				cpos = random.randint(0,allscaffolds[crand]-1)
				# print >> sys.stderr, 'choice', crand, cpos, 'vs original'#, row[:3], row[3:]
				
				newrow = [crand,str(cpos),G[crand]['seq'][cpos]]+row[3:]
				print >> fho, '\t'.join(newrow)
			except KeyError:
				print >> sys.stderr, 'Error: corrupted row: \n%s'%(row)
				pass
		except ValueError:
			print >> sys.stderr, 'Error: could not parse line:', line[:-1]
	fh.close()
	fho.close()
	
	# sort it
	print 'Sorting shuffled likelihood file...'
	call = 'sort -k 1,1 -k 2,2n -S 50%'+' --parallel=%s "%s" > "%s"'%(nthreads, newf,newfsort)
	print 'Call:', call
	os.system(call)
	os.system('rm "%s"'%(newf))
	print 'Done'
