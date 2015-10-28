#! /usr/bin/env python

import re, sys, os, stat
from glimmrAccessories import *

version = """\nglimmrToMatrix.py, version %s
 
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


indirfile = ''
reffile = ''
reflengths = ''
value = 'Q'
outfile = None

help = """Usage:
glimmrToMatrix.py -i <glimmr/dir/or/file/> -r refgenome.fa | -lengths chromlengths.txt [-o outfile.txt -value Q]

Convert a glimmr peak_calls file to genome-wide chromosome x position matrix.
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
	elif arg == 'o' or arg == 'outfile': outfile = val
	elif arg == 'ref' or arg == 'r': reffile = val
	elif arg == 'lengths': reflengths = val
	elif arg == 'value': value = val
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)

if not indirfile or not (reffile or reflengths): sys.exit(help)

scafLengths = {}

sys.stderr.write('Loading reference...'); sys.stderr.flush()
if reffile:
	G = readFasta(reffile)
	for k in G.keys(): scafLengths[k] = len(G[k]['seq'])
	sys.stderr.write('n=%s scaffolds\n'%(len(G.keys()))); sys.stderr.flush()
elif reflengths:
	scafLengths = dict(map(lambda x: (x[1],int(x[0])), readTable(reflengths,header=0,rownames=0)[0]))
	sys.stderr.write('n=%s scaffolds\n'%(len(scafLengths.keys()))); sys.stderr.flush()

allscaffolds = {}
for k in scafLengths.keys(): allscaffolds[k] = {}


scafcount = 1
fh = open(indirfile)

# special treatment of first line
line = fh.next()
row = line[:-1].split('\t')
chrom, pos, ref, model, Q, D, rat1, rat2 = row[:8]
currScaf = chrom
currPileup = {}

for line in fh:
	# just add info to current scaffold until we reach the next scaffold
	
	row = line[:-1].split('\t')
	# C3619755        21      A       nuc     36      4       1       3.0     nuc:2.298e-03,3,0,0,0.0 inp:4.596e-07,-,-,-,-
	
	# assumes a glimmr posterior file format
	chrom, pos, ref, model, Q, D, rat1, rat2 = row[:8]
	try: pos = int(pos)
	except ValueError: continue
	
	theval = None
	if value == 'Q': theval = float(Q)
	elif value == 'nuc.inp': theval = float(rat1)
	elif value == 'ptm.nuc': theval = float(rat2)
	elif value == 'lik': theval = -float(model) # - log likelihood
	else: sys.exit('no matching value: %s'%(value))
	
	# print 'test', chrom,pos,theval
	
	# dump the current scaffold's data whenever we change chromosomes
	if chrom != currScaf:
		if chrom not in scafLengths: 
			print >> sys.stderr, 'Error unrecognized name', chrom
			continue
		
		# create a list from 0 to scaffold length
		tmp = [0 for i in xrange(scafLengths[currScaf])]
		
		for cpos in currPileup.keys():
			# no addition
			try: 
				tmp[cpos] += currPileup[cpos]
				# if tmp[pos] > 160: print 'issue!', tmp[pos]
			except IndexError: 
				print >> sys.stderr, 'Storing error: chrom=%s, currScaf=%s, cpos=%s, len=%s'%(chrom,currScaf,cpos,len(tmp))
		allscaffolds[currScaf] = '\t'.join([currScaf, '\t'.join(map(str, tmp))])
		# done with the information in this scaffold
		
		print >> sys.stderr, 'Saving %s (%s scaffolds completed)'%(currScaf, scafcount)
		
		# update given new chrom
		currScaf = chrom
		currPileup = {pos:theval}
		scafcount += 1
	else:
		try: currPileup[pos] = theval
		except ValueError: print >> sys.stderr, 'error', line[:-1]

# dump last row... only if scafcount is 1
# ------------------------------------------
print >> sys.stderr, 'DUMPING', currScaf, scafcount
# if scafcount == 1:
# create a list from 0 to scaffold length
tmp = [0 for i in xrange(scafLengths[currScaf])]

for cpos in currPileup.keys():
	# no addition
	try: 
		tmp[cpos] += currPileup[cpos]
		# if tmp[pos] > 160: print 'issue!', tmp[pos]
	except IndexError: 
		# print >> sys.stderr, 'TESTING', chrom, currScaf, scafLengths[currScaf], len(tmp)
		print >> sys.stderr, 'Storing error: chrom=%s, currScaf=%s, cpos=%s, len=%s'%(chrom,currScaf,cpos,len(tmp))
			
allscaffolds[currScaf] = '\t'.join([currScaf, '\t'.join(map(str, tmp))])
# done with the information in this scaffold
	
# update given new chrom
currScaf = chrom
currPileup = {pos:theval}
scafcount += 1
print >> sys.stderr, 'Finished processing %s scaffolds'%(scafcount)
# ------------------------------------------

fh.close()

# save entire matrix
# do we need to tab out each row to same length for R?
scafcount = 1
print >> sys.stderr, 'Saving matrix...'
outfh = sys.stdout
if outfile: outfh = open(outfile,'w')
for scaf in sorted(allscaffolds.keys()):
	if allscaffolds[scaf] == {}: print >> outfh, scaf
	else: print >> outfh, allscaffolds[scaf]
	scafcount += 1
outfh.close()
if outfile: outfh.close()
