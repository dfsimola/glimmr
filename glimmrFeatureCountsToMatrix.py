#! /usr/bin/env python

import os, re, sys, math, time, stat, time, random, subprocess, decimal

sys.path.insert(0, "/home/dsimola/glimmr.v1.8/")
import glimmrAccessories as ga

help = """Usage: %s [options]

"""%(sys.argv[0])

"""%s"""%(help)

def pipeit(s, nl=0, pipe=sys.stderr):
	pipe.write(s+(nl>0 and ''.join(['\n' for z in xrange(nl)]) or '')); pipe.flush()

def slash(astr):
	if not len(astr): return '/'
	elif astr[-1] !='/': astr+='/'
	elif astr[len(astr)-2:len(astr)] == '//': astr = astr[:-1]
	return astr

DRYRUN = False

indirfile = ''
outfile = ''

# CLI ARGS
# -------------------------------------------------------------------
help = '\nUsage: %s --in/-i <FILE> [--out/-o <OUTFILE>]\n'%(sys.argv[0])
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
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)


# ga.createdir(outdir)

files = ga.getFiles(indirfile, exclude='.summary')

# thread vectors together - assume same ordering of rows
dt = []
rownames = []
columns = []

for f in files:
	print >> sys.stderr, '- Loading %s'%(f)
	d,r,c = ga.readTable(f, rownames=0)
	
	# sid = ga.getLabel(f).)
	
	# Geneid	Chr	Start	End	Strand	Length	DS68871_AG16844.unique.map
	# chunk0001_1_1	chr1	16140	16200	+	61	1
	# chunk0001_2_1	chr1	51868	52040	+	173	0
	# chunk0001_3_1	chr1	57280	57354	+	75	1
	
	columns += [c[6].strip('.unique.map')]
	
	if not len(rownames):
		rownames = map(lambda row: '%s:%s-%s'%(row[1],row[2],row[3]), d)
		
		
		# rownames = r
	
	Y = ga.vslice(d,6)
	dt += [Y]


if len(outfile):
	ga.printTable(ga.transpose(dt),['Count']+columns,rownames, file=outfile)
else:
	ga.printTable(ga.transpose(dt),['Count']+columns,rownames)



