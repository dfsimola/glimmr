#! /usr/bin/env python

import re, sys, os, stat, glob
from glimmrAccessories import *

version = """\ncleanGlimmrLikFile.py, version %s
 
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

# USER PREFS
# -------------------------------------------------------------------
indirfile = ''
headerfile = ''
reffile = ''
outdir = ''
PRESORT = False
DRYRUN = False

# CLI ARGS
# -------------------------------------------------------------------
help = """\nUsage: cleanGlimmrLikFile.py -i <FILE_OR_DIR> [-f <HEADER_LIST>] [-r <REFFILE] [-o <OUTDIR>]

Helper script used to remove duplicate rows from a glimmr likelihood/mother file if needed. 
File must be presorted by chrom and pos: sort -k 1,1 -k 2,2n <FILE>
"""

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
	elif arg == 'f': headerfile = val
	elif arg == 'r': reffile = val
	
	elif arg == 'sort': PRESORT = True; ai-=1
	elif arg == 'dryrun': DRYRUN = True; ai-=1
	
	elif re.match('^help|h$', arg.lower()): sys.exit(help)
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)


# # USER PREFS
# # -------------------------------------------------------------------
# assert len(sys.argv) == 3, 'Usage: cleanGlimmrLikFile.py <FILENAME> <HEADERS> > <STDOUT>'
# indirfile = sys.argv[1]
# headerfile = sys.argv[2]

curscaf = None
seen = {} # postitions

validheaders = {}
if reffile:
	G = readFasta(reffile)
	validheaders = makeDict(G.keys())
	del G
else:
	validheaders = makeDict(readTable(headerfile,header=0)[1])

print 'Valid headers', len(validheaders), validheaders.keys()[:10]

files = []
if isdir(indirfile) or os.access(indirfile,os.F_OK):
	files = getFiles(indirfile)
else:
	files = glob.glob(indirfile)

createdir(outdir)

for f in files:
	
	if PRESORT:
		sys.exit('not yet implemented')
	
	if os.access(outdir+getFilename(f), os.F_OK):
		print >> sys.stderr, 'File exists, skipping: %s'%(outdir+getFilename(f))
		continue
	
	if DRYRUN: 
		print ' IN: %s'%(f)
		print 'OUT: %s'%(outdir+getFilename(f))
		continue
	
	fh = open(f)
	fho = open(outdir+getFilename(f),'w')
	
	
	for line in fh:
		row = line[:-1].split('\t')
		try:
			chrm,pos = row[:2]
			if curscaf != chrm: 
				curscaf = chrm
				seen = {} # only need to maintain for current scaffold
				print >> sys.stderr, 'Current:', chrm,
				if chrm in validheaders: print >> sys.stderr, 'VALID'
				else: print >> sys.stderr, 'BAD'
			
			if pos not in seen and chrm in validheaders: print >> fho, line,
			# else: print >> sys.stderr, 'SEEN:', line[:-1]
		
			seen[pos] = None
		except ValueError:
			print >> sys.stderr, 'Error: could not parse line:', line[:-1]
	fh.close()
	fho.close()
