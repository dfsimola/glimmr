#! /usr/bin/env python

import re, sys, os, stat, glob
from glimmrAccessories import *

version = """\nglimmrChecklist.py, version %s
 
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
outfile = ''

BATCH = False

# CLI ARGS
# -------------------------------------------------------------------
help = '\nUsage: glimmrChecklist.py ...\nScan a glimmr log file directory for completed files'
nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

# if len(args)==1: sys.exit('\n'.join(help.split('\n')))

if argstr.count('-h') or argstr.count('-help'):
	print >> sys.stderr, help
	sys.exit()


while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('in|^i$', arg): indirfile = val; nhelps += 1
	elif re.match('out|^o$', arg): outfile = val
	elif arg == 'batch': BATCH = True; ai-=1
	elif re.match('^help|h$', arg.lower()): sys.exit(help)
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)

# I/O WORK
# -------------------------------------------------------------------

indirs = [indirfile]
if BATCH: indirs = getDirectories(indirfile)

for indirfile in indirs:
	files = []
	if isdir(indirfile) or os.access(indirfile,os.F_OK):
		files = getFiles(indirfile)
	else:
		files = glob.glob(indirfile)
		
	# DO WORK
	# -------------------------------------------------------------------
	
	print >> sys.stderr, indirfile
	
	total = 0
	done = []
	notdone = []
	
	outfile = '/var/tmp/grepout.txt'
	for f in files:
		lab = getLabel(f).split('.')[-1]
		call = "grep 'Done proc' '%s' > '%s'"%(f,outfile)
		# call = "grep 'Results for' '%s' > '%s'"%(f,outfile)
		os.system(call)
		lst = readList(outfile)
		# print 'test', lab, lst
		if len(lst):
			done += [lab]
		else:
			notdone += [lab]
			
		total += 1
	
	print >> sys.stderr, 'COMPLETE\t%s'%(len(done))
	print >> sys.stderr, 'INCOMPLETE\t%s'%(len(notdone))
	print >> sys.stderr, 'TOTAL\t%s'%(total)
	
	printList(notdone)
	if outfile: 
		if getSuffix(outfile) != 'txt': createdir(outfile)
		printList(notdone,file=outfile+'.%s.txt'%(getLabel(outfile)))















