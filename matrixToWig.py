#! /usr/bin/env python

import re, sys, os, stat, glob
from glimmrAccessories import *

version = """\nmatrixToWig.py, version %s
 
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
outdirfile = ''

format = 'fixedStep'
step = 100 # resolution of track intervals (nt)
tfunc = lambda x: x

# header information
title = ''
description = 'none'
vis = 'full'
color = '128,50,128'
scale = 'on'
wfunc = 'maximum'
smooth = 'off'


# CLI ARGS
# -------------------------------------------------------------------
help = """
Usage: matrixToWig.py ...

Convert glimmr matrix file to wig format
  name              trackLabel           # default is "User Track"
  description       centerlabel          # default is "User Supplied Track"
  visibility        full|dense|hide      # default is hide (will also take numeric values 2|1|0)
  color             RRR,GGG,BBB          # default is 255,255,255
  altColor          RRR,GGG,BBB          # default is 128,128,128
  priority          N                    # default is 100
  autoScale         on|off               # default is on
  alwaysZero        on|off               # default is off
  gridDefault       on|off               # default is off
  maxHeightPixels   max:default:min      # default is 128:128:11
  graphType         bar|points           # default is bar
  viewLimits        lower:upper          # default is range found in data
  yLineMark         real-value           # default is 0.0
  yLineOnOff        on|off               # default is off
  windowingFunction maximum|mean|minimum # default is maximum
  smoothingWindow   off|[2-16]           # default is off
"""

nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

# if len(args)==1: sys.exit('\n'.join(help.split('\n')))

while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('in|^i$', arg): indirfile = val; nhelps += 1
	elif re.match('out|^o$', arg): outdirfile = val
	elif arg == 'format': format = val
	elif arg == 'step': step = int(val)
	elif re.match('^help|h$', arg.lower()): sys.exit(help)
	
	elif arg == 'title': title = val
	elif arg == 'description' or arg == 'desc': description = val
	elif arg == 'color': color = val
	
	elif arg == 'transform': tfunc = eval(val)
	
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)

# I/O WORK
# -------------------------------------------------------------------

# DO WORK
# -------------------------------------------------------------------

span = step


files = []
if indirfile.count('*'): files = glob.glob(indirfile)
elif os.access(indirfile,os.F_OK): files = getFiles(indirfile)
else: sys.exit('cannot access %s'%(indirfile))

if len(files) > 1:
	outdirfile = slash(outdirfile)
	createdir(outdirfile)

for infile in files:
	# header information
	if not title: title = getLabel(infile)
	header = 'track type=wiggle_0 name="%s" description="%s" visibility=%s color=%s autoScale=%s windowingFunction=%s smoothingWindow=%s'%(title,description,vis,color,scale,wfunc,smooth)
	
	sys.stderr.write('Working: %s\n'%(infile)); sys.stderr.flush()
	fho = None
	if len(files) == 1 and outdirfile[-1]!='/': fho = open(outdirfile,'w')
	else: fho = open(outdirfile+getLabel(infile)+'.wig', 'w')
	
	fho.write(header+'\n'); fho.flush()
	
	# stream process the matrix file line by line
	fh = open(infile)
	for line in fh:
		sys.stderr.write('Processing...'); sys.stderr.flush()
		row = line[:-1].split('\t')
		chrom = row.pop(0)
		sys.stderr.write('%s. '%(chrom)); sys.stderr.flush()
		
		# process profile for this chromosome into bins of step size
		start = 0
		if len(row) == 1: continue
		profile = [r3(avg(map(tfunc,row[i:i+step]),ret=0)) for i in xrange(start,len(row),step)]
		if not len(profile): continue
		sys.stderr.write('Saving wig track...'); sys.stderr.flush()
		
		# declaration (note positions are 1-offset)
		stuff = '%s chrom=%s start=%s step=%s span=%s\n'%(format,chrom,start+1,step,span)
		for i in range(len(profile)): stuff += '%s\n'%(profile[i])
		fho.write(stuff)
		
		sys.stderr.write('Done.\n'); sys.stderr.flush()
		
	fh.close()
	fho.close()
	
