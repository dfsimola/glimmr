#! /usr/bin/env python

# brief methods for allele specific mapping
# given reference genome with assumed heterozygous characters (iupac), glimmr will first map reads that do not overlap the iupac sites. then glimmr will generate 2 reference genomes disambiguating the heterozygous sites into one of two files randomly, and map remaining reads to each of these haploid references allowing missing characters. the three map files will then be merged, assuring that, in the case of multiple valid alignments of the same read to the same exact position in the genome, the best alignment is chosen using bowtie2's alignment score (AS:i:<N>).


GLIMMR_VERSION = 1.8

# Version info
# - added ability to map to a heterozygous genome natively or produce allele specific maps for SE or PE data
# - fixed between sample coverage normalization, correcting effective coverage of translikelihoods
# changed default buffersize to reduce disk access
# - added new runtime flags:
#	- --standard [new default]: --lambda 1 --d 1
#	- --complete: --lambda 0.67 --d 10
# - hugly "/path/to/ref/chroms/*.fa": calls --rest chr1,chr2,...,chrN while specifying only the appropriate chrom.fa file to load into memory (only works with --d 1)

"""GlimmrHet: genomic locus identification using multiply maped reads and a diploid (heterozygous) reference genome
 
 Version: %s
 
 Daniel F. Simola, PhD (simola@upenn.edu)
 Laboratory of Shelley L. Berger, PhD
 University of Pennsylvania
 October 2015
 
 Copyright (c) 2015, Daniel F. Simola and Shelley L. Berger, University of
 Pennsylvania.  All Rights Reserved.
 
 You may not use this file except in compliance with the terms of our
 License. You may obtain a copy of the License at XXX
 
 Unless required by applicable law or agreed to in writing, this
 software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 CONDITIONS OF ANY KIND, either express or implied.  See the License
 for the specific language governing permissions and limitations
 under the License.
"""%(GLIMMR_VERSION)

import os, re, sys, math, time, stat, time, random, subprocess, decimal
from operator import add,mul
from string import *

# GLOBAL VARIABLES AND FUNCTIONS
# -----------------------------------------------------------------
# -----------------------------------------------------------------
multichromtoken = '.'         # separator for --rest <chr1.chr2.....chrn> separation
nastr = '"NA"'                # represents missing value
infy = decimal.Decimal('inf') # 
defaultprec = 25              # 28 is the default for the class
bufferlimit = 10000          # number of lines to buffer in memory (10MB)
PCR_HASH_WIDTH = 5
PCR_HASH_NMAX = 5-1

getLabel = lambda astr: '.'.join(astr.split('/')[-1].split('.')[:-1])
getFile = lambda astr: astr.split('/')[-1]
getPath = lambda astr: '/'.join(astr.split('/')[:-1])+'/'
getSuffix = lambda astr: astr.split('/')[-1].split('.')[-1]

all_jobs = []

def launchProcesses(jobs, nsubprocs=1, sleep=1):
	# global all_jobs, running_jobs, nlaunced, ncomplete, nkilled
	
	global all_jobs
	all_jobs += jobs
	
	# multiprocessing queue
	# ---------------------
	njobs = len(all_jobs) # initial job count
	nlaunched = 0
	running_jobs = []
	ncomplete = 0
	nkilled = 0
	
	while len(all_jobs):
		# as long as there are more jobs to be processed
		if nlaunched - ncomplete - nkilled < nsubprocs:
			print >> sys.stderr, ' - Scheduling job', nlaunched+1, '(%s completed, %s killed | x%s processors)'%(ncomplete, nkilled, nsubprocs)
			# and we can launch another job
			# fetch the next job
			call = all_jobs.pop(0)
			# and execute it
			job = subprocess.Popen(call, shell=True)
			running_jobs += [job]
			nlaunched += 1
			time.sleep(sleep)
		else:
			# pipeit('waiting...')
			# check if existing job has been killed
			nkilled = 0
			marked = []
			for ji in xrange(len(running_jobs)): 
				if running_jobs[ji].poll() != None: 
					# print 'poll killed', rj[ji].poll(), len(running_jobs)
					ncomplete += 1
					marked += [running_jobs[ji]]
			for m in marked:
				running_jobs.pop(running_jobs.index(m))
					
			# if nkilled > 0: print '%s jobs killed'%(nkilled)
			# check to see if an existing job has finished
			# ncomplete = len(getFiles(completedir, include='likelihood'))
			time.sleep(sleep)
		# print 'tet', nlaunched
	for job in running_jobs: job.wait() # halt code until all jobs are finished.
	return 1

def mergeSamFiles(lst,outfile,header=False):
	for x in lst:
		assert os.access(x, os.F_OK), 'Cannot access %s'%(x)
	
	f1 = lst.pop(0)
	# assume first file has desired sam header
	if header==True:
		os.system('cp "%s" "%s"'%(f1,outfile))
	else:
		fho = open(outfile,'a')
		fh = open(f1, 'r')
		for line in fh:
			if line[0] != '@': print >> fho, line,
		fh.close()
		fho.close()
		
	# read subsequent lst files and concatenate non-header lines to outfile
	fho = open(outfile,'a')
	for l in lst:
		fh = open(l, 'r')
		for line in fh:
			if line[0] != '@': print >> fho, line,
		fh.close()
	fho.close()
	

def getFilename(astr):
	return astr.split('/')[-1]

def indexFile(fn, outfile):
	"Records number of bytes into file for first row of each chromosome"
	
	if not os.access(fn, os.F_OK): sys.exit('Error: cannot index file %s. Check file access.'%(fn))
	
	dct = {} # keyed by chromosome
	
	curchrom = None
	
	count = 1; nbytes = 0
	
	fh = open(fn)
	for line in fh:
		ll = len(line)
		if line[0] == '@': nbytes += ll; continue
		
		try:
			lst = line[:-1].split('\t')
			# chrom = lst[2]
			# pos = int(lst[3])-1
			
			if lst[2] not in dct:
				dct[lst[2]] = nbytes
			nbytes += len(line)
		# exception if a line beginning with @
		except IndexError: nbytes += ll; pass
		except ValueError: nbytes += ll; pass
	fh.close()
	
	printTable(sorted(dct.items()), file=outfile)
	return 1


def indexDegFile(infile, outbase, filterids={}):
	"Creates a multi-alignment index and a metaindex for this file. Format allows stream-processing"
	
	if not os.access(infile, os.F_OK): sys.exit('Error: cannot index file %s. Check file access.'%(infile))
	pipeit('indexing...')
	
	# helper function for fast printing
	printme = lambda pipe, astr: lambda tup: pipe.write('%s\t%s\t%s\n'%(tup[0], tup[1], astr))
	
	# output file
	fidx = open('%s.index.unsort'%(infile), 'w')
	
	# temporary storage to summarize info for all deg alignments for a given hitid
	curHitid = None; curPese = None; curInfo = []; curByte = []
	count = 1; nbytes = 0 # update line by line
	
	fh = open(infile)
	for line in fh:
		if count % 1000000 == 0: pipeit('%.1e,'%(count))
		count+=1; ll = len(line)
		
		try:
			lst = line[:-1].split('\t')
			chrom = lst[2]
			if chrom == '*': nbytes += ll; continue
			
			hitid = hash(lst[0]); pos=int(lst[3])-1; seq=lst[9]; qual=lst[10]; rl = len(seq)
			
			# must record all instances of multiply-mapped reads (relevant to bounds)
			# but for each entry store first line for a new readid and number of aligments
			# this allows us to avoid the otherhits dictionary
			
			if curHitid == None:
				# So this is the first read we encounter in the file
				curHitid = hitid # name of the multiply-mapped read
				curPese = pese # store in PE or SE part of dictionary
				curInfo = [(chrom,pos,rl)] # relevant positions for this read
				curByte = [(nbytes,ll)] # minimal info for this read's alignments
				
			elif hitid == curHitid:
				# continuing with multiple alignments of the same read
				curInfo += [(chrom,pos,rl)] # add this alignment
				curByte += [(nbytes,ll)]
				
			else:
				# Then we moved on to a new read. Save previous information
				
				# don't save PCR duplicates
				if curHitid in filterids: continue
				
				# Write a single line with all of this information to each chromosome index file
				# iterate over alignments for this read
				# single str representation of all alignments for this read, wrt this chr,pos
				# ccb = ','.join(map(lambda x: '%s:%s'%(x[0],x[1]), curByte[:maxaln]))
				ccb = '%s\t%s'%(curByte[0][0], len(curByte))
				map(printme(fidx,ccb), curInfo)
				
				# for cc,cp,clen in curInfo: # store each alignment, pointing them to the block of bytes
				# 	fidx.write('%s\t%s\t%s\n'%(cc, cp, ccb))
				
				# update
				curHitid = hitid
				curPese = pese
				curInfo = [(chrom,pos,rl)]
				curByte = [(nbytes,ll)]
			
		# exceptions if line begins with @
		except IndexError: pass
		except ValueError: pass
		nbytes += ll
	fh.close()
	fidx.close()
	
	# sort the index files spatially
	pipeit('sorting...')
	os.system('sort -k 1,1 -k 2,2n -S 50%'+' --parallel=%s "%s" > "%s"'%(nsubprocs, '%s.index.unsort'%(infile), '%s.index'%(infile)))
	os.system('rm -f "%s"'%('%s.index.unsort'%(infile)))
	
	# create a meta-index file that stores the starting bytes for each chromosome
	fh = open('%s.index'%(infile))
	dct = {}
	nbytes = 0
	for line in fh:
		chrom,pos = line[:-1].split('\t')[:2]
		if chrom not in dct: dct[chrom] = nbytes
		nbytes += len(line)
	fh.close()
	printTable(sorted(dct.items()), file='%s.metaindex'%(infile))
	return 1


def indexDegFileParts(infile, outbase, filterids={}, bufsize=5000000):
	"Creates a multi-alignment index and a metaindex for this file. Format allows stream-processing"
	
	if not os.access(infile, os.F_OK): sys.exit('Error: cannot index file %s. Check file access.'%(infile))
	
	# output file
	# fidx = open('%s.index.unsort'%(infile), 'w')
	finalindex = '%s.index'%(infile)
	fmetaindex = '%s.metaindex'%(infile)
	
	# -------------------------------------------
	# store temp file names, but clear first if exists
	outdir = infile+'.indices/'
	if os.access(outdir,os.F_OK): os.system('rm -rf "%s"'%(outdir))
	createdir(outdir)
	
	Buff = {}; cbuf = 0
	IH = {}
	fh = open(infile)
	try:
		line = fh.next()
		while line[0] == '@':
			# @SQ	SN:chr1	LN:230208
			if line[0:3] == '@SQ':
				chrom = line[:-1].split('\t')[1].split(':')[1]
				IH[chrom] = outdir+'%s.%s.index'%(outbase,chrom)
				Buff[chrom] = []
			line = fh.next()
	except StopIteration: fh.close(); os.system('rmdir "%s"'%(outdir)); return 1
	fh.close()
	# -------------------------------------------
	
	# temporary storage to summarize info for all deg alignments for a given hitid
	curHitid = None; curPese = None; curInfo = []; curByte = []
	count = 1; nbytes = 0 # update line by line
	
	fh = open(infile)
	for line in fh:
		if count % 1000000 == 0: pipeit('%.1e,'%(count))
		count+=1; ll = len(line)
		
		try:
			lst = line[:-1].split('\t')
			chrom = lst[2]
			if chrom == '*': nbytes += ll; continue
			
			hitid = hash(lst[0]); pos=int(lst[3])-1; seq=lst[9]; qual=lst[10]; rl = len(seq)
			
			# must record all instances of multiply-mapped reads (relevant to bounds)
			# but for each entry store first line for a new readid and number of aligments
			# this allows us to avoid the otherhits dictionary
			
			if curHitid == None:
				# So this is the first read we encounter in the file
				curHitid = hitid # name of the multiply-mapped read
				curPese = pese # store in PE or SE part of dictionary
				curInfo = [(chrom,pos,rl)] # relevant positions for this read
				curByte = [(nbytes,ll)] # minimal info for this read's alignments
				
			elif hitid == curHitid:
				# continuing with multiple alignments of the same read
				curInfo += [(chrom,pos,rl)] # add this alignment
				curByte += [(nbytes,ll)]
				
			else:
				# Then we moved on to a new read. Save previous information
				
				# don't save PCR duplicates
				if cfPCRduplicates(filterids,curHitid): continue
				
				# Write a single line with all of this information to each chromosome index file
				# iterate over alignments for this read
				# single str representation of all alignments for this read, wrt this chr,pos
				ccb = map(lambda x: '%s:%s'%(x[0],x[1]), curByte)
				for cc,cp,clen in curInfo: # store each alignment, pointing them to the block of bytes
					
					Buff[cc] += ['%s\t%s\t%s'%(cc,cp,','.join(ccb))]
					cbuf += 1
					# print >> fidx, '%s\t%s\t%s'%(cc,cp,','.join(ccb))
				
				# update
				curHitid = hitid
				curPese = pese
				curInfo = [(chrom,pos,rl)]
				curByte = [(nbytes,ll)]
			
		# exceptions if line begins with @
		except IndexError: pass
		except ValueError: pass
		nbytes += ll
		
		# flush the buffer
		if cbuf >= bufsize:
			for cc in Buff.keys(): 
				if len(Buff[cc]): 
					ftmp = open(IH[cc], 'a')
					for line in Buff[cc]: ftmp.write('%s\n'%(line))
					ftmp.close()
					Buff[cc] = []
			cbuf = 0
		
	fh.close()
	# fidx.close()
	
	# flush once more
	for cc in Buff.keys(): 
		if len(Buff[cc]): 
			ftmp = open(IH[cc], 'a')
			for line in Buff[cc]: ftmp.write('%s\n'%(line))
			ftmp.close()
			Buff[cc] = []
	
	
	# # sort the index files spatially
	# pipeit('sorting...')
	# os.system('sort -k 1,1 -k 2,2n "%s" > "%s"'%('%s.index.unsort'%(infile), '%s.index'%(infile)))
	# os.system('rm -f "%s"'%('%s.index.unsort'%(infile)))
	
	# clear final index
	ftmp = open(finalindex,'w'); ftmp.close()
	
	# sort the index files by position and concatenate to final index file in sorted order
	pipeit('sorting...')
	for chrom in sorted(IH.keys()):
		if os.access(IH[chrom], os.F_OK): os.system('sort -k 1,1 -k 2,2n -S 50% --parallel=%s "%s" >> "%s"'%(nsubprocs, IH[chrom], finalindex))
	os.system('rm -rf "%s"'%(outdir))
	
	# create a meta-index file that stores the starting bytes for each chromosome
	fh = open(finalindex)
	dct = {}
	nbytes = 0
	for line in fh:
		# try:
		# chrom,pos,alns = line[:-1].split('\t')
		chrom = line[0:line.index('\t')]
		if chrom not in dct: dct[chrom] = nbytes
		# except ValueError: pass
		nbytes += len(line)
	fh.close()
	
	ftmp = open(fmetaindex, 'w')
	for chrom,line in sorted(dct.items()): ftmp.write('%s\t%s\n'%(chrom,line))
	ftmp.close()
	
	# printTable(sorted(dct.items()), file=fmetaindex)
	return 1


def loadPCRduplicatesOLD(filename):
	# load the read ids to ignore (due to PCR duplicates)
	badIDs = {}
	try:
		lh = open(filename)
		for line in lh:
			count, badid = line[:-1].split('\t')[:2]
			badIDs[hash(badid)] = None
		lh.close()
	except IOError: pipeit('Warning: Cannot access PCR dup file %s'%(filename),1)
	return badIDs


def intoPieces(x,width,nmax):
	out = []
	i = 0
	npieces = 0
	while i < len(x) and npieces < nmax:
		out += [x[i:i+width]]
		npieces += 1
		i += width
	out += [x[i:]]
	while len(out) < nmax: out += ['']
	return out


def cfPCRduplicates(D,astr):
	global PCR_HASH_WIDTH, PCR_HASH_NMAX
	p = intoPieces(astr, PCR_HASH_WIDTH, PCR_HASH_NMAX)
	# expectation is exception, so make faster using if else
	
	if p[0] in D:
		if p[1] in D[p[0]]:
			if p[2] in D[p[0]][p[1]]:
				if p[3] in D[p[0]][p[1]][p[2]]:
					if p[4] in D[p[0]][p[1]][p[2]][p[3]]:
						print 'filtering', astr, p
						return True
					else: return False
				else: return False
			else: return False
		else: return False
	else: return False
	
	# try: 
	# 	D[p[0]][p[1]][p[2]][p[3]][p[4]]
	# 	return True
	# except KeyError: return False


def loadPCRduplicates(filename):
	# load the read ids to ignore (due to PCR duplicates)
	global PCR_HASH_WIDTH, PCR_HASH_NMAX
	D = {}
	
	def loadIntoHash(count, astr):
		p = intoPieces(astr, PCR_HASH_WIDTH, PCR_HASH_NMAX)
		try: D[p[0]][p[1]][p[2]][p[3]][p[4]] = count
		except KeyError:
			try: D[p[0]][p[1]][p[2]][p[3]] = {p[4]: count}
			except KeyError:
				try: D[p[0]][p[1]][p[2]] = {p[3]: {p[4]: count}}
				except KeyError:
					try: D[p[0]][p[1]] = {p[2]: {p[3]: {p[4]: count}}}
					except KeyError:
						D[p[0]] = {p[1]: {p[2]: {p[3]: {p[4]: count}}}}
						
		# control on accurate storage
		# if ''.join(pieces) != astr: 
		# 	print "HEY!"
		# 	print 'input ', astr, ''.join(pieces) == astr
		# 	print 'output', ''.join(pieces)
		# 	if ''.join(pieces) != astr: print pieces
		# 	print
		
		return
	
	
	try:
		lh = open(filename)
		[loadIntoHash(*line[:-1].split('\t')[:2]) for line in lh]
		lh.close()
	except ImportError: pipeit('Warning: Cannot access PCR dup file %s'%(filename),1)
	return D



def pipeit(s, nl=0, pipe=sys.stdout):
	pipe.write(s+(nl>0 and ''.join(['\n' for z in xrange(nl)]) or '')); sys.stdout.flush()
	
	

def isnan(x):
    return str(x) == str(1e400*0)


def maxlenval(tup):
	mlv = 0
	for i in tup:
		x = len(str(i))
		if x>mlv:mlv = x
	return mlv

def fixspacing(item,nsp):
	v = str(item)
	if len(v) == nsp: return v
	delt = nsp-len(v)
	if delt < 0:
		v = v[0:nsp]
	else:
		for i in range(delt): v+=' '
	return v

def transpose(mat):
	"""Given a 2d array of arrays, will invert the dimensions."""
	newm = []
	for c in range(len(mat[0])):
		newrow = []
		for r in range(len(mat)):
			newrow.append(mat[r][c])
		newm.append(newrow)
	return newm


fixNull = lambda nastr: lambda x: x=='' and nastr or x
def flatten(sequence):
	def rflat(seq2):
		seq = []
		for entry in seq2:
			if seqin([entry]):
				seq.extend([i for i in entry])
			else:
				seq.append(entry)
		return seq
	
	def seqin(sequence):
		for i in sequence:
			## all sequences have '__contains__' in their dir()
			## parentheses present to aid commenting mid-condition
			if ('__contains__' in dir(i) and type(i) != str and type(i) != dict):
				return True
		return False
	
	seq = [sequence][:] # in case parameter isn't already a sequence
	while seqin(seq): seq = rflat(seq)
	return seq


def createdir(dirname):
	if not dirname: return 0
	if not os.access(dirname, os.F_OK):
		try: os.mkdir(dirname); return 1
		except OSError: pass
	return 0

def record(astr, filename, ioMethod='a'):
	if ioMethod == 'a':
		assert os.access(filename, os.F_OK), 'Cannot access '+filename
	
	pmrfh = open(filename, ioMethod)
	print >>pmrfh, astr
	pmrfh.close()
	
	return 1


def trimComments(dat):
	# parse comments
	ndat = []
	for line in dat:
		# only care about non commented lines
		if not re.match('(^\s*#.*$)|(^\s$)', line):
			# parse any commentary on this line
			cleanme = re.split('^(.*)#.*$', line)
			if len(cleanme) > 1: cleanme.pop(0)
			ndat.append(cleanme[0].strip('\n'))
	return ndat

def getFiles(indirfile, include=[], exclude=[]):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	return getContents(indirfile, typ='files', include=include, exclude=exclude)

def getDirectories(indirfile, include=[], exclude=[]):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	files = getContents(indirfile, typ='dirs', include=include, exclude=exclude)
	return map(slash, files)

def getContents(indirfile, typ='all', include=[], exclude=[]):
	if type(include)!=type(list()): include = [include]
	if type(exclude)!=type(list()): exclude = [exclude]
	keepfiles = []
	if len(include) == 0: include = ['']
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
						# fix oddities in item
						item = item.replace('+','\+')
						
						# print 'TEST', item,'vs',f
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
						# else: print 'no match', item,'vs',f,re.match('.*%s.*'%(item), f)
				# else:
					# print 'ODDITY', f
					
			for f in files:
				inexclude = False
				for item in exclude:
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

# def readList(fname, header=False):
# 	crap,thelist,crap = readTable(fname, header=header)
# 	return flatten(thelist)
def readList(fname, dtype=None, delim='\n'):
	fh = open(fname)
	lst = None
	if dtype: lst = map(lambda x: dtype(x[:-1]), fh)
	else: lst = map(lambda x: x[:-1], fh)
	fh.close()
	if delim != '\n': return lst[0].split(delim)
	return lst

def printList(lst, title="", delim="\n", pipe=sys.stdout, file='', sort=False, newline='\n', ioMethod='w'):
	if file: pipe = open(file,ioMethod)
	if sort: lst.sort()
	if title: 
		lst = [title]+lst#print >> pipe, title
	for i in range(len(lst)):
		if i == len(lst)-1:
			print >> pipe, str(lst[i]),
		else:
			print >> pipe, str(lst[i])+delim,
	if newline: print >> pipe, newline,
	if file: pipe.close()

def readTable(filename, header=True, rownames=True, delim="\t", comments=True, keepKey=False):
	fh = open(filename)
	table = fh.readlines()
	fh.close()
	if comments: table = trimComments(table)
	
	tabHeader = []; rowIDs = []; data = []
	if header: 
		tabHeader = table.pop(0).strip('\n').split(delim)
		if not keepKey: tabHeader.pop(0)
		# else: tabHeader = table.pop(0).strip('\n').split(delim)
	for line in table:
		sline = line.strip('\n').split(delim)
		if rownames:
			rowIDs.append(sline.pop(0))
		data += [sline]
	return data, rowIDs, tabHeader

def printTable(tab, header=[], rows=[], delim="\t", newline='\n', file='', pipe=sys.stdout, ioMethod='w'):
	
	if file: pipe = open(file,ioMethod)
	
	if header:
		printList(header, delim=delim, pipe=pipe, newline=newline)
	if len(rows):
		assert len(rows) == len(tab), "Rows and tab must have same length."
		for ri in range(len(tab)):
			r = [rows[ri]]+tab[ri]
			printList(r, delim=delim, pipe=pipe, newline=newline)
	else:
		for r in tab:
			printList(r, delim=delim, pipe=pipe, newline=newline)
	
	if file: pipe.close()

def printFormattedTable(tab, header=[], rows=[], delim=" ", newline='\n', file='', pipe=sys.stdout, nastr='"NA"', ioMethod='w', colSpacing=[]):
	"""not good for tables > say 100 MB"""
	if file: pipe = open(file,ioMethod)
	if tab == [] and file: pipe.close(); return
	tab2 = tab[:]
	if rows:
		for i in range(len(tab)): tab2[i] = [rows[i]]+tab[i]
	
	# need to determine for each column the max length value and space to that
	collens = colSpacing[:]
	if not len(collens):
		# print 'printformattedtable: fixing col spacing'
		temp = transpose(tab2)
		if len(header):
			for i in range(len(temp)):
				r = temp[i]
				r.append(header[i])
				collens.append(maxlenval(r))
		else:
			for r in temp: collens.append(maxlenval(r))
	
	# now create a new tuple with proper spacing
	temp = []; sheadings = []
	if len(header):
		for i in range(len(header)):
			sheadings.append(fixspacing(header[i],collens[i]))
	
	for row in tab2:
		newrow = []
		for i in range(len(row)):
			val = row[i]==nastr and '.' or row[i]
			#if val==nastr: val = '.' # SAS missing value
			newrow.append(fixspacing(val,collens[i]))
		temp.append(newrow)
		
	# pipe = open(file,ioMethod)
	printTable(temp,header=sheadings,pipe=pipe,delim=delim,newline=newline)
	if file: pipe.close()

def nts2iupac(nts):
	X = {'A':'A', 'T':'T', 'C':'C', 'G':'G', 'N':'N', 'U':'U', '-':'-'}
	
	X['GT'] = 'K'
	X['AC'] = 'M'
	X['CGT'] = 'B'
	X['ACG'] = 'V'
	X['CG'] = 'S'
	X['AT'] = 'W'
	X['AGT'] = 'D'
	X['CT'] = 'Y'
	X['AG'] = 'R'
	X['ACT'] = 'H'
	
	X['AA'] = 'A'
	X['TT'] = 'T'
	X['CC'] = 'C'
	X['GG'] = 'G'
	
	nts = ''.join(sorted(map(lambda x: x.upper(), nts)))
	return X[nts]
	
	

def readFasta(faname, usealtid=False, split='^>', idPat='', addins=[], TOUPPER=False, VERBOSE=False):
	"""Parses a multiple fasta file into a dictionary object keyed by fasta identifier. 
	Option to key the dictionary by an alternative id found in the fasta header: altid=True|False.
	If directory of several fasta files is provided as faname, loads them into single dictionary.
	"""
	fadict = {} # return dictionary of fasta entries
	faj = '' # load all information into single data string, faj
	files = getFiles(faname)
	
	for f in files:
		if VERBOSE: print '- Loading %s'%(f)
		fh = open(f)
		fa = fh.readlines()
		fh.close()
		# remove comment lines
		faj += ''.join(filter(lambda x: not re.match(' *#.*', x) and True, fa))
	
	
	# parse by entry (>)
	getentries = re.compile(split, re.MULTILINE)
	faentries = re.split(getentries, faj)
	
	# for some reason the first entry is the null string
	faentries.pop(0)
	
	# parse individual entries
	for entry in faentries:
		# first  has format >name other-header-info
		(header, seq) = entry.split("\n", 1)
		# trim off the '>' character
		theid = ""
		altid = ""
		info = ""
		# further split up the info - for use of alternative identifier
		pattern = '^(\S+) (.*)$'
		if re.match(pattern, header):
			(crap, theid, info, crap) = re.split(pattern, header)
			info = info.strip(' ').strip('\t')
			if len(info) > 0:
				#print "\""+info+"\""
				(crap, altid, info, crap) = re.split('(\S+)(.*)', info)
				info = info.strip(' ').strip('\t')
		else:
			theid = header
			altid = ""
		
		# remove newlines so the sequence data is contiguous
		seq = re.sub("\n", "", seq)
		# remove terminal whitespace
		seq = seq.strip(' ')
		
		if idPat != '':
			crap, theid, crap = re.split(idPat, theid)
			#print 'theid', theid
		
		for pre,suf,pat in addins:
			if re.match(pat, theid):
				# print 'theid', theid, 'pat', pat
				theid = pre+theid+suf
				# print 'new', theid
		
		if TOUPPER:
			seq = seq.upper()
		
		# build the entry
		# add a little check to see if there are repeat entries
		if fadict.has_key(theid):
			print theid, ": This key has already been used to store another fasta entry!"
		else:
			# build the entry for this id
			# allow choice of primary dict key
			if theid and altid and usealtid:
				#print "OK", theid, altid
				fadict[altid] = {'id':theid, 'altid':altid, 'info': info, 'seq': seq}
			else:
				fadict[theid] = {'id':theid, 'altid':altid, 'info': info, 'seq': seq}
	return fadict

def list2fasta(seqs=[], ids=[]):
	dct = dict()
	uid = 0
	for i in range(len(seqs)):
		theid = str(uid)
		try: theid = ids[i]
		except IndexError: pass
		dct[theid] = {'id':theid, 'seq':seqs[i], 'info':''}
		uid += 1
	return dct

def printFasta(fadict, file='', pipe=sys.stdout, key='seq', usealtid=False):
	"""Print a dictionary of fasta entries to a file in mfasta format"""
	keys = fadict.keys()
	keys.sort()
	
	if file: pipe = open(file,'w')
	
	for theid in keys:
		if 'altid' in fadict[theid] and not usealtid:
			header = '>'+theid+' '+fadict[theid]['altid']+' '+fadict[theid]['info']
		elif 'id' in fadict[theid]:
			header = '>'+theid+' '+fadict[theid]['id']+' '+fadict[theid]['info']
		else:
			header = '>'+theid+' '+fadict[theid]['info']
		print >> pipe, header
		print >> pipe, fadict[theid][key]
		print >> pipe
	
	if file: pipe.close()

def fasta2Phylip(dct={}, file='phylipi.txt', infile='', verbose=False):
	"""Convert a fasta file into Phylip interleaved format."""
	# print infile
	d = dct
	
	if infile: d = readFasta(infile)
	
	names = d.keys()
	# fix length
	maxlen = 10+1
	shorts = ['' for i in range(len(names))]
	for i in range(len(names)): 
		shorts[i] = fixspacing(names[i], maxlen)
	blank = fixspacing('',maxlen)
	taxa = len(names)
	# print d
	# print 'test', d[d.keys()[0]]
	alilen = len(d[d.keys()[0]]['seq'])
	seqs = ['' for i in range(taxa)]
	
	if verbose: print shorts
	if verbose: print taxa, alilen
	if verbose: print map(len, seqs)
	enc = 'ascii'
	for i in range(taxa): seqs[i] = d[names[i]]['seq'].encode(enc)
	
	# print out an interleaved phylip formatted file
	pos = 0
	fh = open(file, 'w')
	# alilen = 715
	fh.write(str(taxa)+' '+str(alilen)+'\n'.encode(enc))
	
	for i in range(taxa):
		fh.write(shorts[i]+seqs[i][pos:pos+maxlen]+' '+seqs[i][pos+maxlen:pos+2*maxlen]+' '+seqs[i][pos+2*maxlen:pos+3*maxlen]+' '+seqs[i][pos+3*maxlen:pos+4*maxlen]+'\n'.encode(enc))
	pos = 2*maxlen
	fh.write('\n'.encode(enc))
	
	while pos < alilen:
		for i in range(taxa):
			fh.write(blank+seqs[i][pos:pos+maxlen]+' '+seqs[i][pos+maxlen:pos+2*maxlen]+' '+seqs[i][pos+2*maxlen:pos+3*maxlen]+' '+seqs[i][pos+3*maxlen:pos+4*maxlen]+'\n'.encode(enc))
		pos += 4*maxlen
		fh.write('\n'.encode(enc))
	
	fh.close()
	
	return 1


def unique(lst):
	"""Returns a subset of lst containing the unique elements."""
	d = {}
	for i in lst:
		if i not in d: d[i] = 0
		else: d[i] += 1
	e = d.keys(); e.sort()
	return e


def unzip(thread):
	"""Reverse of the built-in method zip."""
	try:
		if len(thread):
			return [[thread[i][j] for i in range(len(thread))] for j in range(len(thread[0]))]
		else: return []
	except ValueError: return []


eq = lambda s: lambda i: str(s)!=str(i)
filterstr = lambda s: lambda v: filter(eq(s), v)

def argmax(v):
	"""Return argmax of a list"""
	mi = 0
	mv = v[0]
	for i in range(1,len(v)):
		if float(v[i]) > float(mv):
			mv = float(v[i])
			mi = i
	return mi

def argmin(v):
	"Return argmin of a list"
	mi = 0
	mv = v[0]
	for i in range(1,len(v)):
		if v[i] < mv:
			mv = v[i]
			mi = i
	return mi


def slash(astr):
	if not len(astr): return '/'
	elif astr[-1] !='/': astr+='/'
	elif astr[len(astr)-2:len(astr)] == '//': astr = astr[:-1]
	return astr

def makeContigs(intervals):
	intervals.sort()
	# print intervals[0:100]
	# combine overlaps
	coverage = 0 # total bases mapped
	contigs = []
	curr = 1
	while curr < len(intervals):
		bak = curr
		start = intervals[curr-1][0]; stop = intervals[curr-1][1]
		
		while curr < len(intervals) and intervals[curr][0] <= stop:
			if intervals[curr][1] > stop: stop = intervals[curr][1]
			curr += 1
		
		contigs += [(start,stop)]
		coverage += abs(start-stop)
		curr += 1
	# print 'contigs', contigs
	return {'contigs':contigs, 'coverage':coverage}


roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))
floatna = lambda x: x == nastr and nastr or float(x)
intna = lambda x: (x == nastr or x == '' or x == 'NA' and nastr) or int(x)
sumna = lambda lst: sum(filterstr(nastr)(map(floatna, lst)))

def factit(x):
	"""Iterative implementation of factorial function."""
	if x <= 1: return 1
	a = [0 for i in range(x+1)]
	a[0] = 1; a[1] = 1
	for i in range(2,x+1): a[i] = i*a[i-1]
	return a[x]

nchoosek = lambda n,k: factit(n)/(factit(k)*factit(n-k)) # binomial coefficient
prod = lambda lst: reduce(mul, lst) # product function

def avg(lst):
	try: return reduce(add, lst)/float(len(lst))
	except TypeError: return nastr

def avg0(lst):
	try: return reduce(add, lst)/float(len(lst))
	except TypeError: return 0

def wavg(v, w):
	if len(v) == 0: return ret
	if len(w) > 0 and len(w) != len(v): return ret
	return sum([v[i]*float(w[i]) for i in range(len(v))])


def harmonicMean(v, minval=1e-20):
	if len(v) == 0: return 0
	# prevent 0 division
	v = map(lambda x: x==0 and minval or x, v)
	return len(v) / sum(map(lambda x: 1./x, v))

def wHarmonicMean(v, w, minval=1e-20):
	if len(v) == 0: return 0
	if len(w) > 0 and len(w) != len(v): return 0
	# prevent 0 division
	v = map(lambda x: x==0 and minval or x, v)
	# return len(v) / sum(map(lambda x: 1./x, v))
	return sum(w) / sum(map(lambda x,y: y/x, v,w))

def geometricMean(lst):
	try: return reduce(mul, lst)**1/float(len(lst))
	except TypeError: return 0


ssd = lambda lst, mu: reduce(add, map(lambda x: (x-mu)*(x-mu), lst))
def variance(lst):
	try: return ssd(lst,avg(lst))/float((len(lst)-1))
	except TypeError: return nastr

def sd(lst):
	try: return math.sqrt(variance(lst))
	except TypeError: return nastr

def percentile(v=[], k=0.5):
	"""
	Return the value of array that corresponds to the kth percentile of the data.
	"""
	if len(v) <= 0: return 0.0
	temp = map(float, v)
	temp.sort()
	n = len(temp)
	idx = float(k)*(n+1)
	lidx = math.floor(idx)
	hidx = math.ceil(idx)
	if idx < 1.0:
		return float(temp[0])
	elif idx > n-1:
		return float(temp[n-1])
	elif not (idx == lidx or idx == hidx):
		return float(temp[int(lidx)] + temp[int(hidx)])/2.0
	else:
		return float(temp[int(idx)])

median = lambda v: percentile(v,.5)
getQuants = lambda dat: ','.join([str(q)+':'+str(round(percentile(dat, q),4)) for q in [0.05, 0.25, 0.5, 0.75, 0.95]])

# list operations
head = lambda x: x[0]
tail = lambda x: x[len(x)-1]
car = lambda x: x[0]
cdr = lambda x: x[len(x)-1]

# arithmetic operations with missing values
addna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) + float(v)
subna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) - float(v)
multna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) * float(v)
divna = lambda m, v: ((m==nastr or v==nastr or float(v)==0.0) and nastr) or float(m) / float(v)
roundna = lambda e: lambda x: ((x==nastr or e==nastr) and nastr) or round(float(x),int(e))

def makeRange(xmin,xmax,interval=None,n=None):
	if n == None: n = int(abs(xmax-xmin)/float(interval))
	else: 
		interval = (xmax-xmin)/float(n)
		interval += interval/float(n-1)
	return [xmin+interval*i for i in range(n)]


def histogram(dat, nbins=100, bins=[], categorical=False, yfreq=False, logscale=False, logbase=10, shift=0, dtype=None):
	import math
	if logscale: dat = map(lambda x: math.log(x,logbase), dat)
	
	if len(bins) > 0: 
		# print 'problem?'
		nbins = len(bins)
		givenbins = True
	else: givenbins = False
	
	# need to create a histogram
	filt = filterstr(nastr)(dat)
	if dtype: filt = map(dtype, filterstr(nastr)(dat))
	dat = None
	if not len(filt): return [[0,0]]
	
	counts = [0 for i in range(nbins)]
	if not givenbins and not categorical: bins = makeRange(min(filt), max(filt), n=nbins)
	elif not givenbins: bins = unique(filt)
	
	# minor bin correction for plotting
	if shift != 0: bins = map(lambda x: x+shift, bins)
	
	# this is like O(N2) too slow
	if not categorical:
		for j in range(len(filt)):
			k = 0
			# print j, len(filt)
			while k < nbins and float(filt[j]) > bins[k]: 
				k +=1
			if k == nbins: k -= 1
			counts[k] += 1
	else:
		bindct = {}
		for b in bins: bindct[b] = 0
		bk = bindct.keys()
		for j in range(len(filt)): bindct[filt[j]] += 1
		# convert back to a histogram
		bins = []; counts = []
		for key in bindct.keys(): bins += [key]; counts += [bindct[key]]
	
	tmp = []
	intnbins = int(nbins)
	# counts or frequency?
	if yfreq:
		tot = float(sum(counts))
		if logscale: tmp = [[logbase**bins[k],counts[k]/tot] for k in range(intnbins)]
		else: tmp = [[bins[k],counts[k]/tot] for k in range(intnbins)]
	if logscale: tmp = [[logbase**bins[k],counts[k]] for k in range(intnbins)]
	else: tmp = zip(bins,counts)
	
	if categorical: tmp = sorted(tmp)
	
	return tmp


def reHist(hist,bins, freq=False):
	# aggregates an existing histogram in terms of new bins
	
	newhist = dict( map(lambda x: (x,0), bins) )
	nbins = len(bins)
	# which bins corresponds to this key?
	binmap = {}
	
	for k in transpose(hist)[0]:
		newbini = 0
		while newbini < nbins and bins[newbini] < k:
			newbini += 1
		if newbini >= nbins: newbini = nbins-1
		binmap[k] = bins[newbini]
	
	for k,v in hist:
		newhist[binmap[k]] += v
	
	ret = sorted(newhist.items())
	
	#
	if freq:
		tot = float(sum(transpose(ret)[1]))
		# if logscale: tmp = [[logbase**bins[k],counts[k]/tot] for k in range(intnbins)]
		# else: tmp = [[bins[k],counts[k]/tot] for k in range(intnbins)]
		ret = [[bins[k],ret[k][1]/tot] for k in range(len(bins))]
	
	return ret



log = lambda x: math.log(x, 10)
ln = lambda x: math.log(x, math.e)

# dictionary keyed by ascii character returning phred q-score
# ascii 33 to 126 - full range
fastqchars = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '\{', '\|', '\}', '~']
A2Q = dict(zip(fastqchars, range(len(fastqchars))))

# string conversion to integer (potentially a long)
encodeNt = lambda x: x=='A'and'1' or x=='T'and'2' or x=='C'and'3' or x=='G'and'4' or x=='N'and'5'
encodeSeq = lambda seq: int(''.join(map(encodeNt, seq)))
encodeSeqStr = lambda seq: ''.join(map(encodeNt, seq))
decodeNt = lambda x: x=='1'and'A' or x=='2'and'T' or x=='3'and'C' or x=='4'and'G' or x=='5'and'N'
decodeSeq = lambda val: ''.join(map(decodeNt, str(val)))

encodeQ = dict(zip(fastqchars, map(str,range(10,len(fastqchars)+10)))) # every character is 2 integer digits
decodeQ = dict(zip(map(str,range(10,len(fastqchars)+10)), fastqchars))
encodeQual = lambda qual: int(''.join(map(lambda x: encodeQ[x], qual)))
pairUP = lambda qual: map(lambda x: qual[x-1]+qual[x] ,range(1,len(qual),2))
decodeQual = lambda qual: ''.join(map( lambda x: decodeQ[x], pairUP(str(qual)) ))

isNT = lambda x: (x=='A' or x=='T' or x=='C' or x=='G' and True) or False

def rc(seq=''):
	"Translates a nucleotide sequence (including degenerate symbols) into its reverse complement. Always returns uppercase."
	
	# m=(a/c) => k=(t/g)
	# r=(a/g) => y=(t/c)
	# s=(c/g)
	# v=(a/g/c) => b=(t/g/c)
	# w=(a/t)
	# h=(a/t/c) => d=(a/t/g)
	# n=(a/t/g/c)
	srev = lambda c:\
		(c=='A'and'T') or (c=='T'and'A') or (c=='C'and'G') or (c=='G'and'C') or (c=='U'and'A') or\
		(c=='a'and't') or (c=='t'and'a') or (c=='c'and'g') or (c=='g'and'c') or (c=='u'and'a') or\
		(c=='-'and'-') or (c=='M'and'K') or (c=='K'and'M') or (c=='R'and'Y') or (c=='Y'and'R') or\
		(c=='-'and'-') or (c=='m'and'k') or (c=='k'and'm') or (c=='r'and'y') or (c=='y'and'r') or\
		(c=='S'and'S') or (c=='V'and'B') or (c=='B'and'V') or (c=='W'and'W') or \
		(c=='s'and's') or (c=='v'and'b') or (c=='b'and'v') or (c=='w'and'w') or \
		(c=='H'and'D') or (c=='D'and'H') or 'N' or \
		(c=='h'and'd') or (c=='d'and'h') or 'n'
		
	rcseq = map(srev, seq)
	rcseq.reverse()
	rcseq = ''.join(rcseq)
	return rcseq


def getChromsFromHeader(mapfile, fileformat='sam'):
	chrm = {}
	
	if fileformat == 'sam':
		fh = open(mapfile)
		SEEN = False
		for line in fh:
			if line[0:3]=='@SQ':
				crap,scaf,ct = line[:-1].split('\t')
				scaf = scaf[3:]
				ct = int(ct[3:])
				chrm[scaf] = ct
				SEEN = True
			elif SEEN: break
		fh.close()
	elif fileformat == 'fasta':
		G = readFasta(mapfile)
		for k in G.keys():
			chrm[k] = len(G[k]['seq'])
		del G
	return chrm

# matches based on shared chrom,pos,strand triples
def filterReplicateReadsPos(indirfile, outdir, patterns=[], verbose=False, tol=3):
	"Read a SAM formatted map file and record a list of read ids that share the same chrom,pos,strand triple."
	
	fastqchars = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '\{', '\|', '\}', '~']
	A2Q = dict(zip(fastqchars, range(len(fastqchars))))
	aytocue = lambda x: A2Q[x]
	# print 'indirfile', indirfile
	for pattern in patterns:
		# load the reads files corresponding to a single sample id
		files = getFiles(indirfile, include=pattern)
		files = filter(lambda x: re.match('.*(SE|PE).*', x), files)
		print 'files', files
		H = {}
		badchrpos = {}
		name = 'nemo'
		
		for f in files:
			print '- Processing %s'%(f.split('/')[-1])
			
			name = f.split('/')[-1].split('.')[0]
			
			# need to do 2 passes through the read map
			# first determine which positions have multiple bad reads
			# then get the hit ids at these bad positions and choose 1 to keep
			fh = open(f)
			count = 1
			
			for line in fh:
				if count % 1000000 == 0: 
					nbadids = 0
					for rlc in badchrpos.iterkeys():
						for stnd in badchrpos[rlc].iterkeys():
							nbadids += len(badchrpos[rlc][stnd].keys())
					print '   - reached %s | %s bad reads, and counting...'%(count, nbadids)
				# if count == 100000: break
				try:
					lst = line[:-1].split('\t')
					chrom = lst[2]
					if chrom == '*': continue
					
					hitid = lst[0]; flag=lst[1]; pos=int(lst[3])-1
					
					# retrive run,lane parts from hitid
					idparts = hitid.split('-')
					rl = idparts[0]
					
					# retrieve strand
					strand = '+'
					if int(flag) & 16: strand = '-'
					rlchr = '%s:%s'%(rl,chrom)
					
					# initialize if need be
					if rlchr not in badchrpos: badchrpos[rlchr] = {'+':{}, '-':{}}
					if rlchr not in H: H[rlchr] = {'+':{}, '-':{}}
					
					# hash the run,lane,chr,pos,strand info to find replicates
					try: H[rlchr][strand][pos] += 1
					except KeyError: H[rlchr][strand][pos] = 1
					
					count+=1 # we have officially seen this read
					
					if H[rlchr][strand][pos] > tol:
						# this position has a false sampling issue
						try: badchrpos[rlchr][strand][pos] = H[rlchr][strand][pos]
						except KeyError: badchrpos[rlchr][strand] = {pos:H[rlchr][strand][pos]}
					
				# exceptions if line begins with @
				except IndexError: pass
				except ValueError: pass
			fh.close()
			print '   - evaluated %s reads'%(count-1)
			
		histdat = [] # list of the shared start counts across positions
		for rlchr in H.iterkeys():
			for strand in H[rlchr].keys():
				# this is a dictionary of positions
				histdat += H[rlchr][strand].values()
		del H
		
		lhd = len(histdat) # this is the number of loci with overlapping reads
		havg = 0; hmedian = 0; hist = None
		bins = [1,2,3,4,5,6,7,8,9,10]+range(15,105,5)
		if lhd > 0:
			havg = sum(histdat)/float(lhd)
			hmedian = median(histdat)
			hist = histogram(histdat,bins=bins)
			del histdat
			
			sys.stdout.write('   Saving redundant genomic positions...'); sys.stdout.flush()
			outfile = outdir+'%s_redundantLoci.txt'%(name)
			of = open(outfile,'w')
			for rlchr in badchrpos.iterkeys(): 
				#                   x = strand, {pos:int, pos:int, ...}
				for (strand,posdct) in badchrpos[rlchr].iteritems():
					# pos dict is a dictionary of pos:int items
					for pos,ct in posdct.iteritems():
						print >>of, '\t'.join(map(str, [rlchr, pos, strand, ct]))
			of.close()
			sys.stdout.write('Sorting...'); sys.stdout.flush()
			os.system('sort -n -k 4 -r "%s" -o "%s"'%(outfile, outfile))
			sys.stdout.write('Done\n'); sys.stdout.flush()
		# ----------------------------------------------------
		
		
		# second pass over the map files
		print '   Saving redundant reads...'
		seen = {}
		nmappedreads = 0 # number of mapped reads
		nbadids = 0 # reset counter
		
		outfile = outdir+'%s_redundantIDs.txt'%(name)
		of = open(outfile,'w')
		
		for f in files:
			fh = open(f)
			count = 1
			for line in fh:
				if count % 1000000 == 0: print '   - reached %s | %s bad reads, and counting...'%(count, nbadids)
				
				try:
					lst = line[:-1].split('\t')
					chrom = lst[2]
					if chrom == '*': continue
					
					hitid = lst[0]; pos=int(lst[3])-1; flag=lst[1] 
					qual=lst[10]
					# summarize this into single value
					Q = sum(map(aytocue, qual))
					
					# retrive run,lane parts from hitid
					idparts = hitid.split('-')
					rl = idparts[0]
					
					strand = '+'
					if int(flag) & 16: strand = '-'
					
					rlchr = '%s:%s'%(rl,chrom)
					
					count+=1
					nmappedreads += 1
					
					if rlchr in badchrpos and pos in badchrpos[rlchr][strand]:
						# this is a position with error
						keystr = '%s,%s,%s'%(rlchr,pos,strand)
						# but we want to keep 1 of these reads, so track the number of reads we store
						# this is our memory bottleneck, potentially
						try: seen[keystr] += [(Q,hitid)]
						except KeyError: seen[keystr] = [(Q,hitid)]
						
						if keystr in seen and len(seen[keystr]) == badchrpos[rlchr][strand][pos]:
							# this is the last read to be seen that is part of the current bad position
							bads = seen[keystr]
							# now choose 1 read not to flag, and flag the rest
							# choose the one with highest quality
							bads.sort()
							bads.pop()
							for Q,anid in bads: 
								# we should also save chrom and start location for degenerate reads
								print >> of, '%s\t%s\t%s\t%s'%(badchrpos[rlchr][strand][pos], anid, chrom, lst[3])
								nbadids += 1
								
							# which means we can delete this position from seen now that we have recorded it
							del seen[keystr]
							# and we can delete it from the bad positions
							del badchrpos[rlchr][strand][pos]
							
				# exceptions if line begins with @
				except IndexError: continue
				except ValueError: continue
			fh.close()
			
		del seen
		del badchrpos
		
		sys.stdout.write('   - Sorting output...'); sys.stdout.flush()
		os.system('sort -n -k 1 -r "%s" -o "%s"'%(outfile, outfile))
		sys.stdout.write('\n'); sys.stdout.flush()
		print '   - Final: %s bad reads'%(nbadids)
		
		if lhd > 0:
			failing = nbadids/float(nmappedreads)
			passing = (nmappedreads - nbadids)/float(nmappedreads)
			leg = '%s: tol=%s; %s pass (%.3f), %s fail (%.3f); Avg=%.2f, Median=%.2f (n=%s)'%(name,tolerance,nmappedreads-nbadids,passing,nbadids,failing,havg,hmedian,lhd)
			printTable(hist,file=outdir+'%s_%sread_error_distribution.txt'%(name,tolerance))
			
			# how to replot this if needed
			try: 
				import plot
				# bins = [1,2,3,4,5,6,7,8,9,10]+range(15,105,5)
				# hist = readTable(outdir+'%s_%sread_error_distribution.txt'%(name,tolerance),header=0,rownames=0)
				# hist = hist[0]
				plot.hist(histdat=hist, xlabel='Number of reads sharing chrom,pos,strand', ylabel='Frequency of distinct genomic loci', bins=bins, legend=leg, file=outdir+'%s_%sread_error_distribution.pdf'%(name,tolerance), rotateLabels=-45, custom='set size ratio .33; set yrange [0:1.1]; set grid', logscale='', yfreq=1, showstats=0)
				plot.hist(histdat=hist, xlabel='Number of reads sharing chrom,pos,strand', ylabel='Frequency of distinct genomic loci', bins=bins, legend=leg, file=outdir+'%s_%sread_error_distribution_log.pdf'%(name,tolerance), rotateLabels=-45, custom='set size ratio .33; set yrange [1e-8:5]; set grid', logscale='y', yfreq=1, showstats=0)
			except ImportError: pass
	return 1


# matches based on shared sequence - exact sequence
def filterReplicateReads(indirfile, outdir, name, ud, mm, T=2, seed=[0,36], verbose=False, redo=False, replot=False):
	"Read a SAM formatted map file and record a list of read ids that share the same sequence exactly (on either strand). Keep tol of the instances."
	
	mmfrac = {1:1.0, 2:2/3., 3:1/3.}
	
	import random
	
	pf = outdir+'%s_%s_T=%s_read_error_distribution.txt'%(name,ud,T)
	outfile = outdir+'%s_%s_redundantIDs.txt'%(name,ud)
	
	mbins = [1,2,3,4,5,6,7,8,9,10]+range(15,105,5)+range(150,550,50)
	
	if replot:
		try: 
			import plot
			
			plotdir = outdir+'plots/'
			createdir(plotdir)
			
			tolerance = []; nbadids = []; nmappedreads = []
			nmmb = []; passing = []; failing = []
			
			lhd = []
			hist = [[x,0] for x in mbins]
			
			for ud1 in ud:
				pf = outdir+'%s_%s_T=%s_read_error_distribution.txt'%(name,ud1,T)
				try:
					d,r,c = readTable(pf,header=0)
					nrow = r.index('Name')
					
					thist = zip(map(int,r[:nrow]),[int(i[0]) for i in d[:nrow]])
					for i in range(len(thist)): hist[i] = [mbins[i],hist[i][1]+thist[i][1]]
					
					name_,nbadids_,nmappedreads_,nmmb_,passing_,failing_,havg_,hmedian_,lhd_ = [i[0] for i in d[nrow:nrow+9]]
					
					nmmb_ = int(nmmb_); passing_ = float(passing_)
					nbadids_ = int(nbadids_); failing_ = float(failing_); havg = float(havg_)
					hmedian_ = float(hmedian_); lhd_ = float(lhd_); nmappedreads_ = int(nmappedreads_)
					
					nbadids += [nbadids_]; nmappedreads += [nmappedreads_]
					nmmb += [nmmb_]; passing += [passing_]; failing += [failing_]
					
					lhd += [lhd_]
					
					# passpct = sum(nmmb)/float(sum(nmappedreads))
					
					# passtotal += sum(nmmb)
					# passdenom += sum(nmappedreads)
					
					# havgtotal += sum(nmappedreads)
					# havgdenom += sum(lhd)
					
				except IOError: pass
			
			passpct = divna(sum(nmmb), sum(nmappedreads))
			failpct = subna(1,passpct)
			passpct = roundna(3)(passpct)
			failpct = roundna(3)(failpct)
			
			havg = roundna(2)(divna(sum(nmappedreads),sum(lhd)))
			
			sg = sum(nmmb); sb = sum(nbadids)
			
			print name.replace('_','-'),sg,passpct,sb,failpct,havg,sg+sb
			
			
			leg = '%s: %s pass (%s), %s fail (%s); Avg=%s (n=%d)'%(name.replace('_','-'),sg,passpct,sb,failpct,havg,sg+sb)
			
			plot.hist(histdat=hist, xlabel='Fragment multiplicity (Number of identical reads per fragment)', ylabel='Frequency of fragments', bins=mbins, legend=leg, file=plotdir+'%s_T=%s_read_error_distribution.pdf'%(name,T), rotateLabels=0, custom='set size ratio .33; set yrange [0:1.1]; set grid', logscale='', yfreq=1, showstats=0,xticrounding=int)
			plot.hist(histdat=hist, xlabel='Fragment multiplicity (Number of identical reads per fragment)', ylabel='Frequency of fragments', bins=mbins, legend=leg, file=plotdir+'%s_T=%s_read_error_distribution_log.pdf'%(name,T), rotateLabels=0, custom='set size ratio .33; set yrange [1e-8:5]; set grid', logscale='y', yfreq=1, showstats=0,xticrounding=int)
			
		except ImportError: pass
		except IOError: 
			print >> sys.stderr, '  Cannot access %s'%(pf)
			pass
		return 1
	
	
	elif os.access(outfile, os.F_OK) and not redo: return 1
	
	# touch early to ensure other jobs won't start working on this
	of = open(outfile, 'w')
	
	# global fastqchars;
	# A2Q = dict(zip(fastqchars, range(len(fastqchars))))
	# aytocue = lambda x: A2Q[x]
	
	# load the SE AND PE maps for this sampleID
	if not PESECOMBO:
		files = getFiles(indirfile, include='^%s\.[S|P]E\.%s'%(name,ud))
		files = filter(lambda x: re.match('.*(SE|PE).*', x), files)
	else:
		# print 'ALL FILES', name, getFiles(indirfile)
		files = getFiles(indirfile, include='%s\.%s\.map$'%(name,ud))
		# print 'FILTERED', files
	
	H = {} # keyed by read sequence; returns dictionary of (hitid,Q) items having same sequence
	# needed in the case that we have multiple lanes for the same sampleID
	premap = {}; map2pre = {}; precount = 1 
	
	if not len(files): pipeit('Warning: no map files found for sample %s.%s'%(name,ud),1)
	
	for f in files:
		pipeit('- Processing %s...'%(f.split('/')[-1]))
		
		# First determine which reads are shared, then choose which to keep
		fh = open(f)
		count = 1
		for line in fh:
			if count % 1000000 == 0: pipeit('%s,'%(count))
			
			# if count == 5000000: break
			
			try:
				lst = line[:-1].split('\t')
				if lst[2] == '*': continue # if no chrom then read is not mapped
				
				hitid = lst[0]
				seq = hash(lst[9][seed[0]:seed[1]]) # compressed seed portion of sequence
				# seq = encodeSeq(lst[9][:seed])
				
				# if we are saving Q values
				# Q = sum(map(aytocue, lst[10])) # summarized quality (INT) for ranking
				
				# retrive run,lane parts from hitid
				# only keep small part of hitid that is variable within run,lane
				idparts = hitid.split('-')
				prepart = idparts[0]+'-'+idparts[1]+'-'
				if prepart not in premap: 
					premap[prepart] = precount
					map2pre[precount] = prepart
					precount += 1
				hitid = '%s|'%(premap[prepart])+'-'.join(idparts[2:])
				
				# if degenerate we should only include ids in top tier?
				# if ud == 'degenerate':
				# 	diffs = int(lst[-1].split(':')[-1])
				# 	if diffs > 0: 
				# 
				
				# store the sequence in dictionary to find replicates
				if seq in H: 
					H[seq] += [hitid]
					# H[seq] += [(hitid,Q)] # might use less memory than dictionary
					# H[seq][hitid] = None # need dictionary to deal with READS not alignments
				else:
					H[seq] = [hitid]
					
					# # this doesn't make sense...
					# # try the reverse complement; this ensures that
					# # duplicates of the same fragment are only recorded once
					# # regardless of mapping orientation
					# rseq = hash(rc(lst[9])[seed[0]:seed[1]])
					# # rseq = encodeSeq(rc(lst[9])[:seed])
					# if rseq in H: H[rseq] += [hitid+'---']
					# else: H[seq] = [hitid+'+++']

					# if rseq in H: H[rseq][hitid] = None
					# else: H[seq] = {hitid:None}
					# if rseq in H: H[rseq] += [(hitid,Q)]
					# else: H[seq] = [(hitid,Q)]
				
				# now it is possible that this read has already occurred (multiply-mapped file)
				# this situation is taken care of since the sequence will hash
				
				count+=1 # we have officially seen this read
				
			# exceptions if line begins with @
			except IndexError: pass
			except ValueError: pass
			
		fh.close()
		pipeit('evaluated %s reads'%(count-1), 1)
		
	
	pipeit('  - Saving redundant reads...',1)
	nbadids = 0 # reset counter
	# number of reads with same sequence over all unique sequences in map file
	histdat = []
	# iterate over the sequences
	count = 1
	if ud == 'unique':
		for seq in H.iterkeys():
			if count % 1000000 == 0: pipeit('   - reached %s | %s bad reads, and counting...'%(count, nbadids),1)
			count += 1
			
			duphitids = H[seq]
			nrss = len(duphitids) # number of reads sharing same sequence
			histdat += [nrss]
			
			if nrss > T:
				# for unique be more conservative and keep only 1 read (pop this to remove from the filter)
				# duphitids.pop(random.randint(0,nrss-1))
				
				# keep the expected number T
				for foo in range(T): duphitids.pop(random.randint(0,len(duphitids)-1))
				
				# all remaining reads are bad - add to list of badids
				# save them to be masked during likelihood calculation
				for hitid in duphitids:
					# reconstruct original hitid
					idx = hitid.index('|')
					ohitid = map2pre[int(hitid[:idx])]+hitid[idx+1:]
					print >> of, '%s\t%s'%(nrss, ohitid)
					nbadids += 1
	else:
		
		for seq in H.iterkeys():
			if count % 1000000 == 0: pipeit('   - reached %s | %s bad reads, and counting...'%(count, nbadids),1)
			count += 1
			
			# from the list of alignments with same sequence, 
			# count up the number of alignments per hitid
			alnForID = {}
			for hid in H[seq]:
				try: alnForID[hid] += 1
				except KeyError: alnForID[hid] = 1
				
			duphitids = alnForID.keys() # this is unique list of reads with same sequence
			nrss = len(duphitids) # number of reads sharing same sequence
			histdat += [nrss]
			
			# Need to determine if this sequence is due to PCR duplication
			# for unique reads, every hitid occurs once and maps to the same location
			# for non-unique reads, each hitid maps to multiple locations
			# however these hitids were identified as having the same sequence
			# so they all map to the same set of locations
			
			# so just choose one hitid to get the number of alignments
			# actually a few rare cases have a mixture e.g. 50, 50, 50, 3, so take max
			Tri = T*max(alnForID.values())
			# T = mmfrac[mm]*max(alnForID.values())
			
			# now we can compare the number of shared sequences nrss (observed potential duplicates)
			# to the number of alignments (expected true reads of same sequence)
			# assuming uniform genomic coverage
			# if observed > expected, then this sequence has duplicates
			# keep the expected number of reads, and record as false the remainder
		
			# the caveat here is that number of alignments depends on k mismatches
			# thus the expected true reads is a function f(n,k) = {k=1 and xn, k=2 and yn, k=3 and zn}
			# where x,y, and z, which indicate the fraction of alignments with 0 mismatches
			# we can estimate x,y, and z as the median over all multiply-mapped reads 
			# forget about this for now and just fix the values as x=1, y=.67, z=0.33
			# given this expectation, we can compute a P-value using Poisson
			# that said, the more conservative thing to do is simply use 
			# the observed number of alignments as the threshold
			# and when k > 1 we aren't being too lenient since the information from this
			# read is proportionally diluted anyway 
			# (the effect of 1 extra read is negligible for n >> 1 alignments)
			
			if nrss > Tri:
				# keep the expected number T
				for foo in range(Tri): duphitids.pop(random.randint(0,len(duphitids)-1))
					
				# now we have badids
				# save them to be masked during likelihood calculation
				flushcount = 0
				for hitid in duphitids:
					# reconstruct original hitid
					idx = hitid.index('|')
					ohitid = map2pre[int(hitid[:idx])]+hitid[idx+1:]
					print >> of, '%s\t%s'%(nrss, ohitid)
					nbadids += 1
					
					flushcount += 1
					if flushcount % 100 == 0: of.flush(); sys.stdout.flush(); flushcount = 0
					
	
	pipeit('  - Final count of %s bad reads'%(nbadids),1)
	
	# OLD
	# if we are saving Q values
	# if lhs > udTol[ud]:
	# 	# then this sequence has duplicates
	# 	bads = []
	# 	for hitid,Q in badhitids:
	# 		idx = hitid.index('|')
	# 		#                  reconstructed hitid
	# 		bads += [(Q,map2pre[int(hitid[:idx])]+hitid[idx:])]
	# 		
	# 	# keep 1 of the reads (the one with highest quality) and toss the rest
	# 	bads.sort(); bads.pop()
	# 	for Q,anid in bads: 
	# 		print >> of, '%s\t%s'%(lhs, anid)
	# 		nbadids += 1
	
	of.close()
	
	del H # no longer need this
	
	pipeit('  - Summarizing sequence multiplicity in a histogram...')
	lhd = len(histdat) # number of unique sequences in map file
	nmappedreads = sum(histdat) # total number of reads in map file
	havg = 0; hmedian = 0; hist = None
	if lhd > 0:
		havg = nmappedreads/float(lhd) # average reads per sequence
		hmedian = median(histdat)      # median reads per sequence
		hist = histogram(histdat,bins=mbins) # histogram of reads per sequence
	del histdat
	pipeit('done',1)
	
	
	if lhd > 0:
		failing = nbadids/float(nmappedreads)
		passing = (nmappedreads - nbadids)/float(nmappedreads)
		
		printTable(hist,file=pf)
		
		tab = [['']]
		tab += [['Name',name]]
		tab += [['Num bad reads',nbadids]]
		tab += [['Total reads',nmappedreads]]
		tab += [['Passing reads',nmappedreads-nbadids]]
		tab += [['Passing pct',passing]]
		tab += [['Failing pct',failing]]
		tab += [['Avg reads per seq',havg]]
		tab += [['Median reads per seq',hmedian]]
		tab += [['Num unique sequences',lhd]]
		printTable(tab,file=pf, ioMethod='a')
	
	# don't bother
	# pipeit('- Sorting output...')
	# os.system('sort -n -k 1 -r "%s" -o "%s"'%(outfile, outfile))
	# pipeit('Done.',1)
	
	return 1




#
def threadListPairs(A,B,inner=True,lena=1):
	pairs = []
	genes = []
	Adct = {}
	for i in range(len(A)): 
		# print 'storing', A[i][1], 'into key', A[i][0]
		Adct[A[i][0]] = A[i][1]
	
	for i in range(len(B)):
		name = B[i][0]
		if name in Adct:
			pairs += [[Adct[name], B[i][1]]]
			# print 'adding', pairs[-1]
			genes += [name]
		elif not inner:
			if lena > 1: pairs += [([nastr for x in xrange(lena)], B[i][1])]
			else: pairs += [(nastr, B[i][1])]
			genes += [name]
	# print 'returning', pairs
	return genes, pairs
#
def threadSetPairs(lst,inner=True,flat=True):
	if len(lst) ==1: 
		n,d = unzip(lst[0])
		d = [[x] for x in d]
		return n,d
	if not len(lst): return [],[]
	
	# if list contents are themselves lists or tuples, make sure they are listed
	if type(lst[0][1][1]) == type(()):
		# print 'we have tuples', lst[0][1]
		# print lst[0]
		lst[0] = [[lst[0][i][0], list(lst[0][i][1])] for i in xrange(len(lst[0]))]
		lst[1] = [[lst[1][i][0], list(lst[1][i][1])] for i in xrange(len(lst[1]))]
	
	ans = threadListPairs(lst[0], lst[1], inner)
	# print 'answer', ans
	
	lena = 2
	for i in range(2,len(lst)):
		ans = threadListPairs(zip(*ans), lst[i], inner, lena=lena)
		lena += 1
	# final flattening of data
	if flat:
		# print 'flattening', ans[1]
		# return [ans[0],[flatten(row) for row in ans[1]]]
		return [ans[0],map(flatten,ans[1])]
	else:
		return [ans[0],ans[1]]
#
def joinSet(lst,inner=True,flat=True):
	return threadSetPairs(lst,inner,flat=flat)



# =========================================================
# PROGRAM INITIALIZATION
# =========================================================

# CLI VARIABLES
# -------------------------------------------------------------------
USESGE = False
SGEMEM = '-l mem_free=8G -l h_vmem=8G'
xflag = 'bychrom'
KEEPALIVE = False
JOINPIECES = False
CLEANJOINPIECES = False
FORCEJOIN = False
SGEPYTHONPATH = '/ifs/h/simola/core/simolacode/build/lib.linux-x86_64-2.4'
USESWALLOW = False; swallowpath = 'swallow'

python = 'python'
glimmrpath = ''
glpath = slash(getPath(sys.argv[0]))+'glimmrLikelihoodHet.v%s.py'%(GLIMMR_VERSION)
outdir = slash(os.getcwd())

# btdirs = ['/data15/simola/solseq/bowtie-0.11.3/', '/home/simola/core/bowtie/', '/Users/simola/bin/bowtie-0.12.5/', '/ifs/h/simola/core/bowtie/']
btdirs = []
thdirs = []
cufflinksdir = ''

btmapdir = ''; SETMAPDIR = False
mapfile = ''
zkey = 'glimmr_results' # default folder name for processed samples
bfaref = ''

MAPCUT = False
PID = None

########################### GLIMMR PARAMETERS ###########################
##########################################################################

# PART I: ALIGNMENT PARAMETERS
# ----------------------------
RUNALIGN = False                # align reads to reference genome using bowtie
REALIGN = False                 # realign reads
ALIGNTOOL = 'bowtie2'           # bowtie, tophat, bowtie2
STRATEGY = 'best'                # broad mapping strategy
                                # one of: all, unique, best, bestno [default: all]
                                # unique: genotype using uniqely aligned reads only (--bestm -m 1)
                                # best: genotype using best-guess aligned reads only (--bestk -m 1)
                                # if UNIQ = BEST = False, align using --bestk/m -m d
                                # best-noguess: exclude best reads if multiply best alignments exist

alignmentPolicy = 'bestk'       # bowtie1: basic, best, all, bestk, bestm, maq-bestm, bestm-stratum, best-stratum
                                # bowtie2: all, bestk, bestm
                                
mismatches = {'SE':1, 'PE':1}   # number of mismatches for bowtie alignment (-v mismatches[x])
maxnumalign = 1                 # maximum number of alignments a read could have, 
                                # subject to -k mismathces (ONLY FOR policies bestm and bestk)
                                # set this to the expected read coverage

alignmentMode = '--end-to-end --very-sensitive' # or --local --very-sensitive-local - bowtie2 only

readfileformat = 'fastq'        # or fasta
# qualstr = '--phred64-quals'     # default is 64-base phred scores for Illumina GA 1.3+
qualstr = '--phred33-quals'     # default is 33-base phred scores for Illumina GA 1.3+
trim5 = ''
trim3 = ''

ITER_TRIM = []# [(15,15)]           # default is no iterative re-mapping of trimmed reads
ALIGNTOTRANSCRIPTOME = False    # align unmapped reads to a transcriptome via tophat
gfffile = None                  # gff annotation file for ALIGNTOTRANSCRIPTOME

BOWTIEFLAGS = '' # user specified flags
TOPHATFLAGS = ''# '--segment-mismatches 3 --segment-length 20'


# bowtie1 settings
MAPWITHQUALITY = False          # quality-blind or quality-aware mapping - bowtie1 only
COLORSPACE = False              # use bowtie colorspace
UNIQ = False                    # genotype using uniqely aligned reads only (--bestm -m 1)
BEST = False                    # genotype using best-guess aligned reads only (--bestk -m 1)
                                # if UNIQ = BEST = False, align using --bestk/m -m d
BESTNO = False                  # best-noguess; exclude best reads if multiply best alignments exist


minins = 0; maxins = 300        # need to recalculate - get this from pileup or sam file...
GETINSERTSIZE = False           # estimate PE insert size distributions
insstat = '.99'                 # take as insert bounds the inner insstat percentiles

userAlignArgs = ''              # any additional arguments fed directly to bowtie
COUNTREADS = False              # count number of mapped reads from bowtie output

# PART II: READ FILE ORGANIZATION
# -------------------------------
ORGANIZEREADS = False           # organize bowtie ouput into 4 read files per sampleID: unique/deg/PE/SE
REORGANIZEREADS = False         # redo this read organization
PESECOMBO = True                # if false, save PE and SE alignments to single map file
UDPESECOMBO = False             # if true, save all alignments to a single map file
ALLELESPECIFIC = False          # prepare 2 separate haploid read maps for a diploid reference
FORCEHAPDIP = None              # user can override Diploid setting in map file
EXORGANIZE = False              # organize a preexisting map file
EXOSORT = True                  # by default sort the ex read map
EXPATH = None

KEEPBOWTIE = False              # retain the bowtie output in addition to organized map files...
CLEARBOWTIE = False             # remove bowtie map files after read organization
INDEXMAPFILES = False           # sort unique map files by chrom,pos
INDEXUNI = False                # only index unique reads
INDEXDEG = False                # only index deg reads
RESORT = False                  # resort unique index
INDEXONLY = True                # assume map is already sorted (this is done by -O/-rO)
LINKUNIQUE = False              # if you only have a unique map, this will link file name to what sniper wants
DEGREADSTREAM = True            # stream deg read index

FILTERPCRDUPLICATES = False     # remove putative PCR replicates from unique reads (becomes True if -O is specified)
pcrfilterlength = [0,75]        # how many bases to consider a read as PCR duplicate by matching sequence
                                # default avoids the first several cycles which on average are lower quality
                                # this is a tradeoff between sensitivity and error rate (TP vs FP)
pcrfiltermm = 1                 # number of allowed mismatches within range of pcrfilterlength
pcrfilterT = 10                 # 2 used in simola et al 2013; max num of replicates before we start truncating

replotFD = False                # generate plots for this
redoFD = False                  # redo filtering, replacing existing files
andplotFD = False               # plot after filtering
REDUNDANCYFILTER = True         # mask PCR duplicates during deg file indexing

headerformat = 'sam'            # by default read list of reference scaffolds/chromosomes from header in sam file
boundsfile = ''
outfile = '' # save the posterior probs here
lowestQ = 2000 # based on decimal precision of 200

BYID = False                    # run user initiated glimmrHet.py call for each id in the map file separately
                                # the benefit of this is to decrease storage requirements for large map files

# CUFFLINKS analysis
# ------------------
samheaderfile = ''              # '/Volumes/Mater/genomes/valid.header.cflo.v3.0.sam'
USERPROVIDEDSAM = False         # True/False
RUNCUFFLINKS = False            # perform cufflinks analysis
CUFFLINKSFLAGS = ''             # user provided flags for cufflinks
REDOCUFFLINKS = False           # reprocess cufflinks files even if they already exist in destination folder
CUFFLABEL = ''                  # suffix for cufflinks/ directory

# PART III: posterior computation
# -------------------------------
LIKCALLING = False              # whether or not to perform genotyping
REDOLIKCALLING = False          # whether to redo SNP calling on a sampleID which has already been typed

PEAKCALLING = False
REDOPEAKCALLING = False
REDOFINALPEAKCALLING = False    # correct for multiple testing to return final significant snps

USESINGLETONS = True            # genotype using both PE and SE reads
boundsfile = ''                 # file containing chrom\tlow\thigh rows for restricting genotyping
singlesampleid = []             # name of single sampleID to genotype

LOWSID = None; HIGHSID = None   # if user wants to bound which samples are genotyped in a map file
LIPieces = 10                   # Break loading of unique reads into this number of pieces
                                # greatly reduces memory requirements
nsubprocs = 1                   # number of processors for likelihood/posterior computation
stagger_launch = 1              # seconds between launches

maxreadlen = 150                # maximum read length (for precomputing probabilities)

fragDiameter = 150              # mean of the bioanalyzer size distribution
fragUpperBound = 25             # extra ~1SD from the bioanalyzer size distribution

Lambda = 1.0                    # weighting factor; larger values up-weight global binomial pr.
                                # previous default was 0.67 for dual global + site-specific binomial model
expPhred = 35.                  # global expected per-base error rate

# Qsample = .5                  # fraction of biological sample containing PTM mark

# Pab = 33.7                      # this is the fold enrichment of antibody compared to input
# Perr = 1/float(Pab)             # probability that a read counts as evidence for the wrong hypothesis
Perr = 0.05                     # estimated from ant data (see map file)
Pxab = 1.                       # this is the relative efficiency of the ptm antibody over the nuc antibody
                                # used when computing P(ptm | data) - weight the P(nuc | D) x Pxab
priordef = 0.05                 # default sample probability

DOBNORM = False
betweenSID = []                 # list of other sampleids to consider for read normalization
swapBase = None
betweenBase = []                # for simplicity, just provide core sample name differences

# if one channel has drastically excess coverage may want to subsample randomly to avoid excessive bias
subsampleFold = 1.0 # 1.0 is no change

mapsuffix = '.sam'

# Note you must assume that Perr <= Perr; that is the antibody is better than random

# if blank, only return mother file
PvalueRange = [] #[20, 40, 80]

# OTHER PARAMETERS
# ----------------
READLENBYUSER = False
DIAGNOSTICS = False           # generate summary statistics for genotype results
REDOALL = False               # redo everything
DRYRUN = False                # dryrun the pipeline
VERBOSE = False
ESTIMATEGLOBALERROR = False
MERGESIDS = None              # list of sampleIDs to merge parts into single mother file per sampleID
REMOVEPREMERGED = False       # whether to delete the pre-merged parts after the merge
MERGESORT = False

SAVEALLPOSTERIOR = False      # trigger to report all results regardless of posterior variable state
# DEGCACHESIZE = None           # number of degenerate positions cached before clearing (defaults to glimmrLikelihood value)

HUGLY = None # provide the directory of fasta files for more memory-efficient threading
bounds2Ref = {}

# CLI ARGS
# -------------------------------------------------------------------
help = """
glimmrHet.py v%s

Full help menu: glimmrHet.py --help/-h
Description and requirements: glimmrHet.py --version/-v
Web page: XXX
Latest version: XXX

Usage: glimmrHet.py -r ref.fa (-m mapfile.txt|-id samplename) [options]

ChIPseq analysis with an existing read map file:
   The read map file should be named <mapdir>/samplename.map, where samplename
   is provided in SampleID field of the map file (see below).
   
   If the read map file contains unique reads only (if already sorted, 
   include --nosort to avoid resorting):
   > glimmrHet.py -r ref.fa -i <mapdir> --link --index -id samplename -s
   
   If the read map file contains both unique and non-unique reads, use 
   --exO to partition the file into unique and non-unique portions:
   > glimmrHet.py -r ref.fa -i <mapdir> --exO --index -id samplename -s

Read mapping and genotyping:
   > glimmrHet.py -r ref.fa -m mapfile.txt -z
   
   This will align raw NGS reads, organize and index reads, perform 
   genotyping, and provide diagnostics and new assemblies for all 
   samples listed in the map file. Results will be stored in a new 
   directory 'glimmr_results' in your current working directory.
   
   Note, -z is equivalent to writing 
   "-a -O --index -s -d --saveass -Q 40 --fwer 0.05,50"

Use glimmrHet.py -h for a complete description of arguments.


Standard arguments
==================
   -m <STR>:       map file (see below for format)
   -r <STR>:       reference genome (required, fasta format)
   -i <STR>:       directory to find pre-existing map file(s)
   -o <STR>:       directory in which to save results 
                   (default: current working directory ./)
  --name <STR>:    name of the output/results directory 
                   (default: sniper_results)
   -z:             perform/redo all procedures (equivalent to 
                   -a -O --index -s -d --saveass -Q 40 --fwer 0.05,50)
   -a/-ra:         (re)align reads (generate sam-format read maps with bowtie)
   -O/-rO:         (re)organize read maps into 2 files: 
                     sampleid.unique.map, sampleid.degenerate.map
                   then build an index for each file
  --splitmap       Separate single-end and paired-end reads into different map files
  --link           if you only have a unique map, this links to that file 
                   directly rather than organizing new files
   -exO:           organize a pre-existing read map into the 4 files above 
                   (Default: sampleID.sam, where sampleID is in mapfile)
   -s/-rs:         perform/redo SNP calling (and significance filtering)
   -so/-rso <Q>    perform/redo SNP calling, reporting only SNPs with stringency >= Q
  --sig:           perform significance filtering on existing genotypes
  --resig:         redo significance filtering on existing genotypes
  --diag:          generate diagnostics file summarizing SNP calls
  --saveass:       create a new fasta file of the reference genome modified 
                   by significant SNP calls (IUPAC format). NB, by default
                   -Q 40 is used, but you can return multiple assemblies: 
                   Ex: --saveassembly -Q 20,30,40,100

Key parameters
==============
EX: glimmrHet.py -m <map_file.txt> -r <ref_genome.fasta> -z -k 2 
    -d 10 -e 30 -l 0.5 --prior resequence --theta 0.01

Arguments:
  -k <1,2,3>:      maximum number of alignment mismatches (default: 2)
  -d <INT>:        maximum number of alignments per read (default: 50)
  -e <INT>:        global sequencing error rate (specified as phred value; 
                   E.g. Q=30 == P<0.001) (default: 35)
  -l <FLOAT>:      reweight the likelihood model towards read-specific or 
                   global binomial probability (default: 0.67 towards binomial)
  --prior <STR>:   prior probability distribution (default: resequence)
                   Options: 
                   resequence:     h0:1-theta-theta^2, t1:theta/2, h2:theta/2, t2:theta^2
                   resequence-het: h0:1-theta-theta^2, t1:theta, h2:theta^2/2, t2:theta^2/2
                   uniform:        h0:1-3*theta, t1:theta, h2:theta, t2:theta
                   maq:            h0:(1-theta-theta^2)/2, t1:theta, h2:(1-theta-theta^2)/2, t2:theta^2
                   naive:          h0:0.25, t1:0.25, h2:0.25, t2:0.25
                   haploid:        h0:1-theta, t1:theta
  --theta <FLOAT>: expected divergence (heterozygous and homozygous) from 
                   reference genome (default: 0.001)

Read map indexing
=================
EX: glimmrHet.py -m <map_file.txt> -r <reference_genome.fasta> --index

Arguments:
  --index:         sorts and indexes the unique SE and PE map files
  --reindex:       replaces an existing index
  --(re)indexuni:  index unique maps only
  --(re)indexdeg:  index degenerate maps only
  --nosort:        index only, assuming map files are already sorted

MEMORY / RUNTIME / SCHEDULING options
=====================================
EX: glimmrHet.py -m <map_file.txt> -r <reference_genome.fasta> -s -x 5
For each sampleID, break the genome into 5 pieces, then launch and maintain 
up to 5 parallel processes on your local machine. The total number of jobs 
equals 5 times the number of sampleIDs in the map file. If --sge is added 
schedule jobs across CPUs using Sun Grid Engine.

EX: glimmrHet.py -m <map_file.txt> -r <reference_genome.fasta> -s -x 15 
    --rest <restriction_file.txt>

This will launch and maintain up to 15 parallel processes, where the total 
number of jobs is equal to the number of sampleIDs in the map file times the 
rows in <restriction_file.txt>.
   -x <INT>:    distribute jobs over <INT> processes on 1 machine (multi-core)
  --groupchrom: group multiple chromosomes into a single job (default for -x 1)
  --subchrom:   split each chromosome into -x jobs
  --bychrom:    process each chromosome as a separate file, regardless of -x
                (allows you to add more processors mid-run or recover from
                 errors more quickly.)
  --ids <STR[,STR,...,STR>: restrict genotyping to a subset of SampleIDs
 
  --stream:     stream the degenerate index (Default: True)
  --buffer <INT>: sets streaming buffer to approx. <INT> MB (Default: 10)
 
  --sge:        distribute jobs using Sun Grid Engine. If -x > 1 launch 
                <INT> jobs per sampleID, otherwise launch one job for each 
                different sample to be processed
  --lip <INT>:  partition map into <INT> chunks to lower memory requirements
                (NB, requires <INT> scans of the map file, increasing runtime)
  --sgemem:     Required memory to launch process on SGE (default: 16G)
  --rest:       Specify a region or file containing a list of regions over 
                which to perform genotyping (default: none)
                EX: "--rest chr1"; genotype the entirety of chr1 only
                EX: "--rest chr1.chr5"; genotype chrom1 and chr5 only
                EX: "--rest chr1,0,1000"; genotype positions 0 to 1000 
                    of chr1, inclusive
                EX: "--rest chr1,0,1000,-"; genotype all of chr1 excluding 
                    positions 0 to 1000
                A restriction file should be a tab-delimited list of one or 
                more regions similarly described, replacing commas with tabs.
                Execute this shell command to generate one:
                echo -e "chrom1\\nchrom2\\nchrom3" > myrest.txt, or
                echo -e "chrom1\\t0\\t1000\\nchrom5\\t500\\t50000" > myrest.txt


Other notable arguments
=======================
  --mindepth <INT>: minimum number of reads to evaluate a nucleotide locus
  --maxdepth <INT>: maximum number of reads to evaluate a nucleotide locus.
                    If > <INT> reads only the first <INT> will be used.
  --fasta:       read files do not contain quality scores
  --qual33:      quality scores follow a phred 33 scale instead of default 
                 phred 64 scale
  --peonly:      only use paired end reads for SNP identification
  --all:         use all reads for mapping (based on --policy bestm with 
                 maximum of d alignments)
  --uniq:        only use uniquely mapping reads for SNP identification. Also
                 if provided with -a for mapping
                 only return an alignment if it is unique in the reference 
                 genome given specified maximum mismatches
                 bowtie arguments: '-v mismatches -a -m 1 --best
  --best:        generate a best-guess read map, returning an alignment if it
                 is the one with the fewest mismatches
                 given specified maximum mismatches
                 bowtie arguments: '-v mismatches -k 1 --best'
  --bestno:      similar to --best but discard reads which have multiple 
                 equally best alignments
  --mapqual      use read quality values during mapping (bowtie -n mode with
                 default parameters)
  --colorspace:  input files are in ABI SOLiD colorspace format 
                 (suffix .csfasta and .qual)
                 NB, currently only works with single-end reads
  --ins <FLOAT>: estimate insert size range for paired-end reads as the 
                 <FLOAT> percentile of the empirical distribution 
                 (default: 0.99) (see below for more details)
  --globalerror: estimate expected per-nucleotide sequencing error rate from 
                 your read map
  --cstheta <FLOAT>: Haploid proportion of divergence between sample and 
                     reference genomes (default: 3/2 theta)
  --qual-int:    sets the --integer-quals flag in bowtie
  -Q <INT[,INT,...,INT]>: Comma-delimited list of phred-quality scores for
                          which to compute significance (Default: 40)
  -fwer <FLOAT,INT>: Filter significant SNPs using a family-wise error rate 
                     of <FLOAT> and an average read length of <INT>. This 
                     corrects the posterior probability for each genotyped 
                     locus for multiple testing of all loci that utilize the 
                     same set of NGS alignments. The Q stringency required for
                     significance is -10 log( <FLOAT> / (2 <INT> dbar) ), 
                     where <FLOAT> is the error rate (e.g. 0.01), <INT> is the
                     read length (e.g. 50), and dbar is the estimated mean 
                     number of total alignments for all reads that overlap the
                     current locus of interest (Default: 0.05,50).
  --nominal:     Instead of waiting until the whole genome is processed before
                 returning significant SNPs, this updates a file 
                 sampleID.nominal.txt with potential significant SNPs
                 concurrent with processing.
  --suffix <STR>:    Change default map file suffix to FILENAME.<STR> [Default: sam]

Bowtie mapping options
======================
  --bowtie/-btd: Specify the base directory of your bowtie installation 
                 (e.g. --btd /usr/bin/bowtie-0.12.5) (default: '')
  --keepbowtie/--kb: Do not delete the original bowtie SAM file following 
                     read organization
  --btx <INT>:   number of processors used for mapping
  --policy:      Specify Bowtie alignment policy (default: bestk)
                 Options:
                 * 'basic': 'bowtie -n <mismatches> -m 1 (Maq-like)
                 * 'all': 'bowtie -v <mismatches> -a (end-to-k all alignments)
                 * 'bestk': -v <mismatches> -k <max num alignments> --best
                 * 'bestm': -v <mismatches> -a -m <max num alignments> --best
  --kpe: Permit a maximum number of alignment mismatches for paired-end reads
  --kse: Permit a maximum number of alignment mismatches for single-end reads
  --trimreads n5,n3: For alignment, consider read substring from 1+n5...|r|-n3
  --realign/--ra: perform new Bowtie alignment, replacing existing alignments
  --integer-quals: invokes Bowtie flag of same name (default: NA)
  --btcustom:      arbitrary custom string fed directly to bowtie (default: NA)

Map file format
===============
The map file is a standard tab-delimited text file indicating how each raw 
sequence reads file should be processed by sniper. The format is organized as 
one row per sequenced lane and requires 7 columns of information row, as 
follows:

#!Run Lane SampleID       Alias RefGenome           PairedEnd Path
  1   1    myfirstgenome  foo   /path/to/ref_genome 0         /path/to/run1/
  2   5    mysecondgenome bar   /path/to/ref_genome 1         /path/to/run2/

(If you copy/paste this, make sure to convert spaces to tabs.)

In this example, the first sample was run on lane 1 and contains single-end 
reads. The second sample was sequenced on lane 5 and contains paired-end reads.
The standard mapper used with glimmrHet.py is bowtie, so the specified reference 
is the name of the respective index file, excluding file suffix.

The first line shows the headers but is treated as a comment; 
comments are discarded by Sniper.

Raw sequence read files should contain as a substring the label indicated 
in the Lane column and placed in the /path/to/reads/directory/ directory. 

For example the first sample file may be named 's_1_sequence.txt' or 
'testing1.txt'. The second sample should have 2 files since it is 
paired-end, named 's_5_1_sequence.txt' and 's_5_2_sequence.txt'. 
Colorspace files should take the format 's_1_sequence.csfasta' 
and 's_1_sequence.qual' for the fasta and quality files, respectively. 

Additional fields may be added to the map file after the required columns 
(e.g. concentration, date sequenced, author, etc.).

To specify the prior model and theta parameter value in the map file, 
include headings "Prior" and "Theta" after "Path" and include the 
particular prior and theta values for each sample. Missing values are 
tolerated, in which case those samples will be processed with the 
default values.

Estimating insert size distribution for paired-end data
=======================================================
Two-step execution for paired-end read data:
1. glimmrHet.py -m <map_file.txt> -r <reference_genome.fasta> --ins
2. glimmrHet.py -m <map_file.txt> -r <reference_genome.fasta> -z

By default paired-end reads are mapped using an insert size range of 
0 to 350 nt. Users may specify custom min and max bounds using

  --minins: Specify the minimium insert size for all samples
  --maxins: Specify the maximum insert size for all samples

Alternatively use '--ins <FLOAT>' to estimate the empirical insert size 
distribution for each sample specified in the map file:

  --ins <FLOAT>: Estimate insert size distributions for all samples specified 
                 in the map file using the percentile <FLOAT>[0,1] of the 
                 insert size distribution for alignment (default: 0.99)

Note that alignments must first be generated using '-a'. 
A file 'PE_insert_size.txt' will be created in the specified output directory 
that Sniper will use by default if samples are mapped again 
(using '--realign/--ra'). 
"""%(GLIMMR_VERSION)


version = """\nglimmrHet.py, version %s
January 2012
Daniel F. Simola, Ph.D.
Laboratory of Shelley L. Berger, Ph.D.
University of Pennsylvania
http://www.med.upenn.edu/berger/

About Glimmr: XXX

Requires Unix/Mac OS X/CYGWIN with Python 2.5 or 2.6 
(preferably 64-bit for memory allocation > 4GB)

Copyright (c) 2011, Daniel F. Simola and Shelley L . Berger, University of
Pennsylvania. All Rights Reserved.

You may not use this file except in compliance with the terms of our License. 
You may obtain a copy of the License at 
http://kim.bio.upenn.edu/software/LICENSE

Unless required by applicable law or agreed to in writing, this software is 
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, 
either express or implied.  See the License for the specific language governing
permissions and limitations under the License.

"""%(GLIMMR_VERSION)

nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

if len(args)==1: sys.exit('\n'.join(help.split('\n')[1:36]))

while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('out|^o$', arg): outdir = slash(val)
	elif arg == 'python' or arg == 'python-version': python = val
	elif arg == 'glimmrPath' or arg == 'glimmrpath' or arg == 'glpath': glpath = val
	elif re.match('^help|h$', arg.lower()): sys.exit(help)
	elif re.match('^version|v$', arg.lower()): sys.exit(version)
    
	elif re.match('^map$|^m$', arg): mapfile = val
	elif arg == 'hugly' or arg == 'refpieces' or arg == 'rp': HUGLY = val; xflag = 'subchrom'

	elif re.match('^name$', arg): 
		zkey = val
		if zkey[len(zkey)-1] == '/': zkey = zkey[0:len(zkey)-1]
	elif arg == 'fasta': readfileformat = 'fasta'; ai-=1
	elif re.match('^procs$|^x$', arg): nsubprocs = int(val)
	
	elif re.match('noqual|^nq$', arg): Lambda = -1
	elif arg == 'standard': Lambda = 1.0; maxnumalign = 1; STRATEGY = 'best'; alignmentPolicy = ''; ai-=1
	elif arg == 'complete': Lambda = 0.67; maxnumalign = 10; STRATEGY = 'all'; alignmentPolicy = 'bestk'; ai-=1

	elif re.match('redoall|^z$', arg): REDOALL = True; ai-=1
	elif re.match('dryrun', arg): DRYRUN = True; ai-=1
	elif re.match('verbose|^v$', arg): VERBOSE = True; ai-=1
	elif arg == 'keepalive': KEEPALIVE = True; ai -= 1
	
	elif re.match('^subchrom$', arg): xflag = 'subchrom'; ai-=1
	elif re.match('^bychrom$', arg): xflag = 'bychrom'; ai-=1
	elif re.match('^crosschrom$', arg): xflag = 'crosschrom'; ai-=1
	elif re.match('^groupchrom$', arg): xflag = 'groupchrom'; ai-=1
	elif re.match('^bysampleid$|^bysample$|^byid$', arg): BYID = True; ai-=1
	
	elif arg == 'sge': USESGE = True; ai-=1
	elif arg == 'memsge': SGEMEM = '-l mem_free=%s -l h_vmem=%s'%(val.upper(),val.upper())
	elif arg == 'token': multichromtoken = val
	
	elif re.match('^btdir$', arg): btdirs = [slash(val)]
	elif arg == 'thdir': thdirs = [slash(val)]
	elif arg == 'cuffdir': cufflinksdir = slash(val)
	elif re.match('^readmapdir$|^btmapdir$|^i$', arg): SETMAPDIR = True; btmapdir = val
	
	elif re.match('^mapcut$|^mc$', arg): MAPCUT = True; ai-=1
	
	elif arg == 'staggerlaunch': stagger_launch = int(val)
	elif arg == 'subsample': subsampleFold = float(val)
	elif re.match('^loadinpieces$|^lip$', arg): LIPieces = int(val)
	
	#elif arg == 'merge': MERGESIDS = val.split(','); REMOVEPREMERGED = False
	#elif arg == 'merges': MERGESIDS = val.split(','); REMOVEPREMERGED = False; MERGESORT = True
	#elif arg == 'merge+': MERGESIDS = val.split(','); REMOVEPREMERGED = True
	#elif arg == 'merges+': MERGESIDS = val.split(','); REMOVEPREMERGED = True; MERGESORT = True
	
	elif arg == 'merge': MERGESIDS = True; REMOVEPREMERGED = False; ai-=1
	elif arg == 'merges': MERGESIDS = True; REMOVEPREMERGED = False; MERGESORT = True; ai-=1
	elif arg == 'merge+': MERGESIDS = True; REMOVEPREMERGED = True; ai-=1
	elif arg == 'merges+': MERGESIDS = True; REMOVEPREMERGED = True; MERGESORT = True; ai-=1
	
	elif arg == 'join': JOINPIECES = True; CLEANJOINPIECES = False; KEEPALIVE = True; ai-=1
	elif arg == 'joinclean': JOINPIECES = True; CLEANJOINPIECES = True; KEEPALIVE = True; ai-=1
	elif arg == 'forcejoin': FORCEJOIN = True; ai-=1
	
	elif arg == 'k': mismatches = {'SE':int(val), 'PE':int(val)}
	elif arg == 'kpe': mismatches['PE'] = int(val)
	elif arg == 'kse': mismatches['SE'] = int(val)
	
	elif re.match('^degeneracy$|^d$', arg): maxnumalign = int(val)
	elif re.match('^efficiency$|^per$', arg): Perr = float(val)
	elif re.match('^xefficiency$|^px$', arg): Pxab = float(val)
	elif re.match('^lambda$|^l$', arg): Lambda = float(val)
	elif re.match('^expphred$|^e$', arg): expPhred = float(val)
	
	elif re.match('^ref$|^r$', arg): bfaref = val
	elif re.match('^trimreads$|^trim$', arg): trim5,trim3 = map(int, val.split(','))
	elif re.match('^itertrim$|^iter$', arg): ITER_TRIM = map(lambda x: map(int, x.split('-')), val.split(','))
	
	elif re.match('^align$|^a$', arg): RUNALIGN = True; ai-=1
	elif re.match('^aligntool$',arg): ALIGNTOOL = val
	elif re.match('^realign$|^ra$', arg): REALIGN = True; ai-=1
	
	elif re.match('^rna$', arg): ALIGNTOTRANSCRIPTOME = True; ai-=1
	elif arg.upper() == 'GTF' or arg.upper() == 'GFF': gfffile = val
	
	elif arg == 'tophat': ALIGNTOOL = 'tophat'; TOPHATFLAGS = val
	elif arg == 'tophat2': ALIGNTOOL = 'tophat2'; TOPHATFLAGS = val
	elif arg == 'bowtie2': ALIGNTOOL = 'bowtie2'; BOWTIEFLAGS = val
	
	elif re.match('countreads$', arg): COUNTREADS = True; ai-=1
	
	elif re.match('usesingletons|pese|sepe|^g$', arg): USESINGLETONS = True; ai-=1
	elif re.match('peonly|^pe$', arg): USESINGLETONS = False; ai-=1
	
	elif re.match('^qual33$', arg): qualstr = '--phred33-quals'; ai-=1
	elif re.match('^qual64$', arg): qualstr = '--phred64-quals'; ai-=1
	
	# mapping strategy selection
	elif re.match('^policy$|^pol$', arg): alignmentPolicy = val
	elif re.match('^all$|^total$', arg): STRATEGY = 'all'; ai-=1
	elif re.match('^uniq$|^unique$', arg): 
		# only return an alignment if it is unique in the reference genome
		# given specified maximum mismatches
		# STRATEGY = 'unique'; alignmentPolicy = 'bestm'; maxnumalign = 1; ai-=1

		STRATEGY = 'unique'; alignmentPolicy = 'bestk'; maxnumalign = 2; ai-=1
		# here, we find up to 2 alignments, run -O and then use the X.unique.map file
		# which contains the reads that only reported a single alignment.

	elif re.match('^best$|^bestguess$', arg): 
		# return an alignment if it is the one with the fewest mismatches
		# given specified maximum mismatches; quasi-randomly break ties
		STRATEGY = 'best'; alignmentPolicy = ''; maxnumalign = 1; ai-=1
	elif re.match('^bestno$|^bestnoguess$', arg):
		# return an alignment if it is the one with the fewest mismatches
		# but discard read if multiple equally best alignments exist
		STRATEGY = 'bestno'; alignmentPolicy = 'bestk'; maxnumalign = 2; ai-=1
		# here, we find up to 2 alignments, run -O and discard the read if both alignments
		# have same number of mismatches; this is more conservative than the standard bestguess

	elif arg == 'mapqual': MAPWITHQUALITY = True; ai-=1
	
	elif arg == 'end-to-end': alignmentMode = '--end-to-end --very-sensitive'; ai-=1
	elif arg == 'local': alignmentMode = '--local --very-sensitive-local'; ai-=1
	elif arg == 'more' or arg == 'bt-more' or arg == 'btmore' or arg == 'btcustom': userAlignArgs = val
	
	elif arg == 'ins': GETINSERTSIZE = True; insstat = str(val); RUNALIGN = True; minins=0; maxins=1000
	elif re.match('minins', arg): minins = int(val)
	elif re.match('maxins', arg): maxins = int(val)
	
	elif re.match('organize|^O$', arg): ORGANIZEREADS = True; FILTERPCRDUPLICATES = True; ai-=1
	elif re.match('reorganize|^rO$', arg): 
		ORGANIZEREADS = True; REORGANIZEREADS = True; KEEPBOWTIE = True; RESORT = True;
		FILTERPCRDUPLICATES = True; redoFD = True; ai-=1
	elif arg == 'as' or arg == 'alleles': ALLELESPECIFIC = True; ai-=1
	elif arg == 'haploid': FORCEHAPDIP = 0; ai-=1
	elif arg == 'diploid': FORCEHAPDIP = 1; ai-=1
	elif arg == 'combo': PESECOMBO = True; ai-=1
	elif arg == 'singlemap' or arg == 'mergemap': UDPESECOMBO = True; ai-=1
	elif arg == 'splitpese' or arg == 'splitmap': PESECOMBO = False; ai-=1
	elif re.match('exorganize|^exO$', arg): ORGANIZEREADS = True; EXORGANIZE = True; FILTERPCRDUPLICATES = True; ai-=1
	elif arg == 'exOno': ORGANIZEREADS = True; EXORGANIZE = True; FILTERPCRDUPLICATES = True; EXOSORT = False; ai-=1
	elif arg == 'exOpath': EXPATH = val
	elif arg == 'noexosort' or arg == 'exonosort': EXOSORT = False; ai-=1
		
	elif re.match('^keepbowtie$|^kb$', arg): KEEPBOWTIE = True; ai-=1
	elif re.match('^clearbowtie$|^cb$', arg): CLEARBOWTIE = True; ORGANIZEREADS = True; ai-=1
	elif arg == 'index': INDEXMAPFILES = True; ai-=1
	elif arg == 'reindex': INDEXMAPFILES = True; RESORT = True; ai-=1
	elif arg == 'indexuni': INDEXUNI = True; ai-=1
	elif arg == 'reindexuni': INDEXUNI = True; RESORT = True;ai-=1
	elif arg == 'indexdeg': INDEXDEG = True; ai-=1
	elif arg == 'reindexdeg': INDEXDEG = True; RESORT = True;ai-=1
	elif arg == 'nosort': INDEXONLY = True; ai-=1
	elif arg == 'link': LINKUNIQUE = True; ai-= 1
	
	elif arg == 'stream': DEGREADSTREAM = True; ai -= 1
	elif arg == 'nostream': DEGREADSTREAM = False; ai -= 1
	elif arg == 'buffer': bufferlimit = int(val)*100
	
	elif arg == 'keepdups' or arg == 'nopcr': REDUNDANCYFILTER = False; ai -= 1
	
	elif re.match('filterduplicates|^fd$', arg): FILTERPCRDUPLICATES = True; ai-=1
	elif re.match('nofilterduplicates|^nfd$', arg): FILTERPCRDUPLICATES = False; ai-=1
	elif re.match('plotduplicates|^pfd$', arg): FILTERPCRDUPLICATES = True; replotFD = True; ai-=1
	elif re.match('refilterduplicates|^rfd$', arg): FILTERPCRDUPLICATES = True; redoFD = True; ai-=1
	elif re.match('filterplotduplicates|^fd\+$', arg): FILTERPCRDUPLICATES = True; andplotFD = True; ai-=1
	elif re.match('refilterplotduplicates|^rfd\+$', arg): FILTERPCRDUPLICATES = True; redoFD = True; andplotFD = True; ai-=1
	elif arg == 'pcr-length': pcrfilterlength = map(int, val.split(','))
	elif arg == 'pcr-t': pcrfilterT = int(val)
	elif arg == 'pcr-mm': pcrfiltermm = int(val)
	
	# cufflinks
	elif arg in ['cuff','runcuff','runcufflinks','cufflinks']: RUNCUFFLINKS = True; ai-=1
	elif arg == 'cufflinksflags' or arg == 'cl-flags': CUFFLINKSFLAGS = val
	elif arg == 'redocuff' or arg == 'recuff' or arg == 'redocufflinks': REDOCUFFLINKS = True; ai-=1
	elif arg == 'cufflabel': CUFFLABEL = val
	elif arg == 'samheader': samheaderfile = val; USERPROVIDEDSAM = True
	
	elif re.match('^lik$|^k$', arg): LIKCALLING = True; ai-=1
	elif re.match('^relik$|^rk$', arg): LIKCALLING = True; REDOLIKCALLING = True; ai-=1
	elif re.match('posterior|^post$|^p$', arg): PEAKCALLING = True; ai-=1
	elif re.match('^redoposterior$|^repost$|^rp$', arg): PEAKCALLING = True; REDOPEAKCALLING = True; REDOFINALPEAKCALLING = True; ai-=1
	elif arg == 'prior': priordef = float(val)
	elif re.match('^minquality$|^Q$', arg): PvalueRange = [int(val)]
	elif re.match('^refheader$', arg): headerformat = 'fasta'; ai -= 1
	elif re.match('^restrict$|^rest$', arg): boundsfile = val
	elif re.match('^sampleid$|^id$', arg): singlesampleid = [val]
	elif re.match('^sampleids$|^ids$', arg): singlesampleid = map(lambda x: x.strip(), val.split(','))
	elif re.match('^subset|^sub$', arg): LOWSID,HIGHSID = map(int,val.split(','))
	
	elif re.match('^betsampleids$|^bet$|^bnormids$', arg): betweenSID = map(lambda x: x.strip(), val.split(',')); DOBNORM = True
	elif re.match('^betbase$|^betb$', arg): 
		DOBNORM = True
		tmp = val.split(':'); swapBase = tmp[0].strip(); betweenBase = map(lambda x: x.strip(), tmp[1].split(','))
	elif arg == 'bnorm' or arg == 'bnormpat': betweenBase = val.split(','); DOBNORM = True
		
	elif re.match('^readlen$|^rl$', arg): maxreadlen = int(val); READLENBYUSER = True
	elif re.match('^diameter$|^fraglen$|^fl$', arg): fragDiameter = int(val)
	elif re.match('^delta$', arg): fragUpperBound = int(val)
	
	elif re.match('^globalerror$|^ge$', arg): ESTIMATEGLOBALERROR = True; ai-=1
	
	elif re.match('^diagnostics$', arg): DIAGNOSTICS = True; ai-=1
	# elif re.match('degcache|^dc$', arg): DEGCACHESIZE = int(val)
	
	elif arg == 'saveall': SAVEALLPOSTERIOR = True; ai-=1
	
	elif arg == 'pid': PID = val
	
	elif arg == 'swallow': USESWALLOW = True; swallowpath = val
	# elif arg == 'swallowpath': swallowpath = val
	
	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)
# -------------------------------------------------------------------

if maxnumalign > 1:
	STRATEGY = 'all'; alignmentPolicy = 'bestk';

if HUGLY:
	if maxnumalign > 1: sys.exit('--hugly flag only compatible with --standard or -d 1.')
	import glob
	# want to queue up subchrom calls for each chrom in ref genome, loading only the single chrom at a time

	# by default if i use -subchrom - x 20 without -rest what happens?

	# set boundstr = chr1.chr2.chr3.,,,
	# and set ref to dir/chr1.fa etc, 
	reffiles = glob.glob(HUGLY)
	reffiles = filter(lambda x: '_' not in x, reffiles)
	# reffiles = getFiles(HUGLY,include='fa', exclude='_')
	# print 'REFFILES', reffiles
	
	reflabs = map(getLabel,reffiles)
	boundsfile = '.'.join(sorted(reflabs))
	# print 'boundsfile', boundsfile

	bounds2Ref = dict(zip(reflabs,reffiles))
	# print 'bounds2ref', bounds2Ref

glimmrpath = args[0]
pipeit('\n------> %s'%(getFilename(glimmrpath)),1)
pipeit('------> Path: %s\n'%(glimmrpath),1)

if PID: pipeit('\n'+' '.join(args)+'\n',1)

btdir = ''
for x in btdirs: 
	if os.access(x, os.F_OK): btdir = x; break

thdir = ''
for x in thdirs: 
	if os.access(x, os.F_OK): thdir = x; break


# contingency of flags
# --------------------
if REDOALL:
	RUNALIGN = REALIGN = True
	ORGANIZEREADS = REORGANIZEREADS = True
	FILTERPCRDUPLICATES = True
	INDEXMAPFILES = True; RESORT = True
	LIKCALLING = True; REDOLIKCALLING = True


# if ORGANIZEREADS and not bfaref: 
	# sys.exit('Error: reference fasta file must be provided with -O using --ref/-r <REF.FASTA>.')


# declaration of directories
# --------------------------
zdir = outdir+zkey+'/'

if not SETMAPDIR: btmapdir = zdir+'mapped_reads_all'
unmapdir = zdir+'unmapped_reads_all'
pmadir = zdir+'alignment_mismatches'
tmpstr = 'ALL'

if STRATEGY=='unique':
	if not SETMAPDIR: btmapdir = zdir+'mapped_reads_uni'
	unmapdir = zdir+'unmapped_reads_uni'
	pmadir = zdir+'alignment_mismatches_uni'
if STRATEGY=='best':
	if not SETMAPDIR: btmapdir = zdir+'mapped_reads_best'
	unmapdir = zdir+'unmapped_reads_best'
	pmadir = zdir+'alignment_mismatches_best'
if STRATEGY=='bestno':
	if not SETMAPDIR: btmapdir = zdir+'mapped_reads_bestno'
	unmapdir = zdir+'unmapped_reads_bestno'
	pmadir = zdir+'alignment_mismatches_bestno'

if MAPWITHQUALITY: 
	btmapdir += 'q'
	unmapdir += 'q'
	pmadir += 'q'

btmapdir += '/'
unmapdir += '/'
pmadir += '/'

zinsdir = zdir+'insertsize/'
badReadDir = zdir+'libraryErrors.%s-%snt.%smm.T%s/'%(pcrfilterlength[0],pcrfilterlength[1],pcrfiltermm,pcrfilterT)

tmpstr = 'ALL'
if STRATEGY=='unique': tmpstr = 'UNI'
elif STRATEGY=='best': tmpstr = 'BEST'
elif STRATEGY=='bestno': tmpstr = 'BESTNO'
if MAPWITHQUALITY: tmpstr += 'Q'
resdir = zdir+'results-L%s,E%s,%s,%s/'%(Lambda, expPhred, USESINGLETONS and 'PESE' or 'PE', tmpstr)
tmpdir = resdir+'tmp/'

# file suffix depending on map type
ubamap = '_%s'%(STRATEGY)
if MAPWITHQUALITY: ubamap += 'q'

logdir = zdir+'log_files/'
postlogdir = zdir+'log_files_posterior/'

# create new directories
# ----------------------
createdir(zdir); createdir(logdir); createdir(postlogdir)
createdir(outdir);   createdir(btmapdir)
createdir(zinsdir)
createdir(unmapdir); createdir(pmadir)

# decimal precision for PP computation
# decimal.getcontext().prec = defaultprec

# save a log file with the user's input arguments
logfile = zdir+'log_file.txt'
loglikfile = zdir+'log_lik_file.txt'
iomethod = 'w'
if os.access(logfile, os.F_OK): iomethod = 'a'
fh = open(logfile, iomethod)
print >> fh, ' '.join(args)
fh.close()

# parse map file and group replicates
# ----------------------------------------
if not mapfile:
	sys.exit('Map file must be provided with --map/-m <MAPFILE>. Try glimmrHet.py --help for more details.')

fh = open(mapfile)
runs = []; mapd = []; head = []; postcombos = []
INPOSTCOMBOS = False
ct = 0
for line in fh:
	# store header
	if line[0:2] == '#!': head = map(lambda x: x.strip(), line[2:-1].split('\t'))
	# enter posterior part
	if line[0:2] == '#@': INPOSTCOMBOS = True; continue
	
	if line[0] != '#':
		# trim past any place a # hash is present
		cline = line[:-1]
		# cline = cline.strip('\x00')
		try: cline = cline[:cline.index('#')]
		except ValueError: pass
		tmp = cline.split('\t')
		
		if not INPOSTCOMBOS:
			if len(tmp) > 1:
				runs += [tmp.pop(0)]
				mapd += [tmp]
		else:
			# this parses the enrichment compared to input
			if len(tmp) == 2 and tmp[0] == 'Perr':
				Perr = float(tmp[1])
			# this parses the cross antibody relative efficiency
			elif len(tmp) and tmp[0] != '':
				postcombos += [tmp]
	
	ct += 1

# works if you dont have the post combos section
# mapd,runs,head = readTable(mapfile,keepKey=True, header=False)
prioridx = None; rlidx = None; diamidx = None; subfidx = None
idxoffset = 1
if 'Prior' in head: prioridx = head.index('Prior') - idxoffset
if 'ReadLength' in head: rlidx = head.index('ReadLength') - idxoffset
if 'Diameter' in head: diamidx = head.index('Diameter') - idxoffset
if 'Subsample' in head: subfidx = head.index('Subsample') - idxoffset
# sizeidx = head.index('Size') - 1

seqRuns = {} # keyed by sampleID
rawfiles = {} # keyed by UID
allSeqRuns = {} # backup
fragDiameter4Sample = {} # fragDiameter by sample ID
sampleReadLengths = {}
sampleSubsample = {}
ct = 0
totalsamples = 0
rlnotissued = True

for i in range(len(runs)):
	ct += 1
	
	totalsamples += 1
	run,lane,sampleID,alias,refGenome,diploid,paired,path = runs[i], mapd[i][0], mapd[i][1], mapd[i][2], mapd[i][3], mapd[i][4], int(mapd[i][5]), mapd[i][6]
	if FORCEHAPDIP != None: diploid = FORCEHAPDIP
		
	if path == '': path = '/'
	else: path = slash(path)
	
	# make absolute path if local directory specified
	if path[0] != '/' or path == '/':
		abspath = slash(os.getcwd())
		# if this does not end at the current local path, add path
		if slash(abspath.split('/')[-2]) != path: abspath += path
		path = abspath
		
		path = slash(path)
	
	# print 'testing prior', sampleID, prioridx, mapd[i]
	
	# sample-specific prior probability
	try:
		if prioridx: prior = float(mapd[i][prioridx])
		else: prior = priordef
	except IndexError:
		prior = priordef
		pipeit('Warning: Prior not specified for sample %s... defaulting to %s'%(sampleID,prior),1)
	except ValueError:
		prior = priordef
		pipeit('Warning: Prior not specified for sample %s... defaulting to %s'%(sampleID,prior),1)
	
	
	# fetch sample-specific read length
	if rlidx:
		# print 'WTF', mapd[i], 'OK', mapd[i][rlidx]
		sampleReadLengths[sampleID] = int(mapd[i][rlidx].strip())
		# print sampleID, 'readlength', sampleReadLengths[sampleID]
	elif READLENBYUSER: 
		sampleReadLengths[sampleID] = int(maxreadlen)
	elif rlnotissued and LIKCALLING: 
		pipeit('Warning: read length not specified in map file or by command line (--rl INT). Required for posterior computation.',1)
		rlnotissued = False
	
	# sample-specific subcoverage
	if subfidx: 
		try: sampleSubsample[sampleID] = float(mapd[i][subfidx].strip())
		except IndexError: sampleSubsample[sampleID] = subsampleFold
		except ValueError: sampleSubsample[sampleID] = subsampleFold
	else: sampleSubsample[sampleID] = subsampleFold
	
	# check for sample-specific bandwidth parameter; take average of multiple samples
	try:
		if diamidx: 
			if sampleID in fragDiameter4Sample: fragDiameter4Sample[sampleID] += [int(mapd[i][diamidx].strip())]
			else: fragDiameter4Sample[sampleID] = [int(mapd[i][diamidx].strip())]
		else: fragDiameter4Sample[sampleID] = [fragDiameter]
	except IndexError:
		fragDiameter4Sample[sampleID] = [fragDiameter]
		pipeit('Warning: Diameter not specified for sample %s... defaulting to %s'%(sampleID,fragDiameter),1)
	except ValueError:
		fragDiameter4Sample[sampleID] = [fragDiameter]
		pipeit('Warning: Diameter not specified for sample %s... defaulting to %s'%(sampleID,fragDiameter),1)
	
	if MAPCUT and (ct < LOWSID or ct > HIGHSID): continue
	
	info = [refGenome, int(diploid), paired, path, prior, alias]
	try: seqRuns[sampleID][(run,lane)] = info
	except KeyError: seqRuns[sampleID] = {(run,lane):info}
	
	try: allSeqRuns[sampleID][(run,lane)] = info
	except KeyError: allSeqRuns[sampleID] = {(run,lane):info}
	
pipeit('- Read %s samples from map file.'%(totalsamples),1)

# average fragdiam
for k in fragDiameter4Sample.keys(): fragDiameter4Sample[k] = int(sum(fragDiameter4Sample[k])/len(fragDiameter4Sample[k]))

# ---------------------------------------------------------
sortseqruns = seqRuns.keys()
sortseqruns.sort()

if not MAPCUT:
	keeps = dict(map(lambda x: (x,None), sortseqruns[LOWSID:HIGHSID]))
	for sid in sortseqruns:
		if sid not in keeps:
			del seqRuns[sid]
			

if len(singlesampleid):
	sr2 = {}
	for sid in singlesampleid:
		if sid in seqRuns:
			bak = seqRuns[sid]
			sr2[sid] = bak
		else: pipeit('Warning: sample %s not found in map file.'%(sid),1)
	seqRuns = sr2

sortseqruns = seqRuns.keys()
sortseqruns.sort()

# =========================================================
# Part 0: Make sub calls if -byid was issued
# =========================================================
if BYID:
	jobs = []
	idcount = 1
	byidlogdir = zdir+'byid.runlog/'
	createdir(byidlogdir)
	
	NPROCS_BYID = 1
	NPROCS_BYCHR = 1
	
	if nsubprocs < len(sortseqruns):
		# we are processor limited, fastest is to distribute across ids / read maps
		NPROCS_BYID = nsubprocs
		NPROCS_BYCHR = 1
	else:
		# here we can launch 1 job per ID and potentially multiple jobs per ID
		NPROCS_BYID = len(sortseqruns)
		NPROCS_BYCHR = max(nsubprocs / int(len(sortseqruns)), 1)
	
	for anID in sortseqruns:
		# remove the by id flag !!!
		call = '%s '%(python)+' '.join(args)
		call = replace(call,'--byid', '')
		call = replace(call,'-byid', '')
		# also need to compute new CPU usage since we are distributing
		
		call += ' -x %s'%(NPROCS_BYCHR)
		call += ' --ids %s >> "%s"'%(anID,byidlogdir+'%s.runlog.txt'%(anID))
		pipeit('%s| SUBCALL>> %s'%(idcount, call),2)
		jobs += [call]
		idcount += 1
	
	# old - why?
	# launchProcesses(jobs,1)
	# new - 
	launchProcesses(jobs, NPROCS_BYID)
	sys.exit('DONE')


if LOWSID or HIGHSID: pipeit('- Processing samples %s'%(sortseqruns[LOWSID:HIGHSID]), 1)
else: pipeit('- Processing %s samples'%(len(sortseqruns)),1)
# ---------------------------------------------------------

# load bounds file
boundslabel = ''
seqbounds = {}
if boundsfile and os.access(boundsfile, os.F_OK):
	d,r,c = readTable(boundsfile, header=False)
	for i in range(len(d)):
		tmp = []
		try: tmp += [intna(d[i][0])]
		except IndexError: tmp += [nastr]
		except ValueError: tmp += [nastr]
		try: tmp += [intna(d[i][1])]
		except IndexError: tmp += [nastr]
		except ValueError: tmp += [nastr]
		try: tmp += [d[i][2]] # flag - + or -
		except IndexError: tmp += ['+']
		keyname = r[i].strip()
		try: seqbounds[keyname] += [tmp]
		except KeyError: seqbounds[keyname] = [tmp]
		
	parts = boundsfile.split('/')
	crap, boundslabel, crap = re.split('(.+)\..+$', parts[len(parts)-1])
	boundslabel = '.'+boundslabel
	
elif boundsfile:
	# this is a string formatted in a few possible ways
	# ex1: chr1[;chr2; ...;chrx]
	# ex2: chr1,500,1000,+
	# ex3: chr1,500,NA,+
	# ex4: chr1
	
	splt = boundsfile.split(multichromtoken)
	if len(splt) > 1:
		# then we have ex1
		for loc in splt:
			seqbounds[loc] = [(nastr, nastr, '+')]
		boundslabel = '.%s--%s'%(splt[0],splt[-1])
	else:
		# ex 2 or 3, or 4
		splt = boundsfile.split(',')
		if len(splt) == 1:
			loc = splt[0]
			seqbounds[loc] = [(nastr, nastr, '+')]
			# boundslabel = '.%s,%s,%s,%s'%(loc,nastr,nastr,'+')
			boundslabel = '.%s'%(loc)
		elif len(splt) == 3:
			loc, low, high = map(lambda x: x.strip(' '), splt)
			seqbounds[loc] = [(intna(low), intna(high), '+')]
			boundslabel = '.%s,%s,%s,%s'%(loc,low,high,'+')
		elif len(splt) == 4:
			loc, low, high,flag = map(lambda x: x.strip(' '), splt)
			seqbounds[loc] = [(intna(low), intna(high), flag)]
			boundslabel = '.%s,%s,%s,%s'%(loc,low,high,flag)
		else:
			sys.exit('Restriction flag not specified properly: '%(boundsfile))


# merge likelihood files if split
# -------------------------------
if MERGESIDS:
	if MERGESIDS == True: MERGESIDS = sortseqruns
	# print 'SIDS', MERGESIDS
	# sys.exit()
	if ALLELESPECIFIC:
		ms2 = []
		for sid in MERGESIDS:
			ms2 += [sid+'.h1']
			ms2 += [sid+'.h2']
		MERGESIDS = ms2
		
	for sid in MERGESIDS:
		filesToMerge = getFiles(resdir, include='%s.mother'%(sid), exclude=['%s.mother%s.txt'%(sid,boundslabel),'.bz','.zip','.bz2'])
		filesToMerge.sort()
		if not len(filesToMerge): pipeit('Warning: no files for %s found to merge. SKIP.'%(sid),1)
		mergedFile = resdir+'%s.mother%s.txt'%(sid,boundslabel)
		fh = open(mergedFile,'w'); fh.close()
		
		if MERGESORT:
			# first sort the contents of each partial file
			pipeit('- Sorting %s partial files...'%(len(filesToMerge)),1)
			for fx in xrange(len(filesToMerge)):
				pipeit('  %s) %s'%(fx+1, '/'.join(filesToMerge[fx].split('/')[-1])),1)
				os.system('sort -k 1,1 -k 2,2n -S 50% --parallel=%s "%s" > "%s"'%(nsubprocs, filesToMerge[fx], filesToMerge[fx]+'.sort'))
			
			pipeit('- Merging %s files (in above order) into %s'%(len(filesToMerge),mergedFile),1)
			for afile in filesToMerge:
				os.system('cat "%s" >> "%s"'%(afile+'.sort',mergedFile))
			
			# remove originals
			if REMOVEPREMERGED:
				for afile in filesToMerge: os.system('rm -f "%s"'%(afile)); os.system('rm -f "%s"'%(afile+'.sort'))
				
		else:
			pipeit('- Merging %s files (in order) into %s'%(len(filesToMerge),mergedFile),1)
			for fx in xrange(len(filesToMerge)): pipeit('  %s) %s'%(fx+1,filesToMerge[fx]),1)
			for afile in filesToMerge:
				os.system('cat "%s" >> "%s"'%(afile,mergedFile))
			
			# remove originals
			if REMOVEPREMERGED:
				for afile in filesToMerge: os.system('rm -f "%s"'%(afile))
	sys.exit('Done merging sampleIDs: %s'%(','.join(MERGESIDS)))

# DONE INITIALIZATION

if len(sortseqruns)  ==  0: sys.exit('- No samples to process. Done.')

# =========================================================
# Part I: Align reads to reference genome
# =========================================================
parts = mapfile.split('/')
peinsfile = outdir+'PE_insertsize_table-%s'%(parts[len(parts)-1])
peinsize = {} # dict specifying PE insert size statistics, keyed by sampleID
if os.access(peinsfile, os.F_OK):
	d,r,c = readTable(peinsfile, header=True)
	for i in range(len(d)):
		peinsize[r[i]] = {}
		for j in range(len(c)):
			peinsize[r[i]][c[j]] = float(d[i][j])


# save a sampleID to fastq files table
for sampleID in sortseqruns:
	for run,lane in seqRuns[sampleID].keys():
		refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
		UID = '%s-%s-%s'%(sampleID,run,lane)
		# find the file that is indicated in the map file
		fastqFiles = sorted(getFiles(path, include=lane, exclude=['pdf','zip']))
		
		if paired:
			# remove extraneous files that do not have the _1 or _2 in them
			fastqFilesTmp = filter(lambda x: ('_1' in x or '_2' in x) and True, fastqFiles)
			if not len(fastqFilesTmp):
				fastqFiles = filter(lambda x: ('R1' in x or 'R2' in x) and True, fastqFiles)
				del fastqFilesTmp
			else:
				fastqFiles = fastqFilesTmp
			
			# print 'test', fastqFiles
			if len(fastqFiles) > 2:
				print >> sys.stderr, 'Warning: Found more than 2 matching read files. Trying default s_lane_x_sequence.txt format.'
				fastqFiles = [path+'s_%s_1_sequence.txt'%(lane), path+'s_%s_2_sequence.txt'%(lane)]
			elif len(fastqFiles) < 2:
				print >> sys.stderr, 'Skipping: Check read file names or path. Did not find 2 files for sample %s in path "%s": %s'%(lane, path, fastqFiles)
				continue
			
			if not os.access(fastqFiles[0], os.F_OK) or not os.access(fastqFiles[1], os.F_OK): 
				print >> sys.stderr, 'Skipping: could not find some raw read files. Check paths.'
				continue
		else:
			if len(fastqFiles) > 2:
				print >> sys.stderr, 'Warning: Found more than 2 matching read files. Trying default s_lane_x_sequence.fastq format.'
				tmpfiles = ['s_%s_sequence'%(lane)]
				fastqFiles = sorted(getFiles(path, include=tmpfiles, exclude=['pdf','bz2','zip']))
				fastqFiles = filter(lambda x: '%s_1'%(lane) not in x and '%s_2'%(lane) not in x, fastqFiles)
		
		try: rawfiles[sampleID] += fastqFiles
		except KeyError: rawfiles[sampleID] = fastqFiles

# save mapping of raw reads to sampleIDs to file
printTable([[k]+rawfiles[k] for k in sorted(rawfiles.keys())],file=zdir+'readFilesToSampleIDs.txt',ioMethod='w')



if RUNALIGN or REALIGN:
	# --solexa1.3-quals: internal quality conversion from solexa to phred64
	# -t: print time elapsed
	# --pairtries <int>: try aligning a PE read more times than default (100)
	# --minins <int>: minimum insertion size (default=0).
	# --maxins <int>: maximum insertion size (default=250).
	# --n <int>: number of mismatches allowed in seed region (5' 28 bases)
	# -m <int>: toss alignment if it has > m allowable alignments to the reference
	# --sam: output in SAM format instead of bowtie reads
	
	readMapPipeFile = zdir+'read_mapping.log'
	readMapPipe = open(readMapPipeFile, 'a')
	
	pipeit('READ MAPPING (%s)...'%(ALIGNTOOL), 1)
	# call bowtie on all samples in mapfile
	for sampleID in sortseqruns:
		
		# check to see if we have completed maps
		MAPSARECOMPLETE = True
		if ALIGNTOTRANSCRIPTOME and not os.access(btmapdir+'%s.UD.sam'%(sampleID), os.F_OK):
			MAPSARECOMPLETE = False
		else:
			for ud in ['unique', 'degenerate']:
				if not os.access(btmapdir+'%s.%s.map'%(sampleID,ud), os.F_OK): 
					MAPSARECOMPLETE = False
		
		# set parm str each time here because insert size may be sample specific
		try: minins = peinsize[sampleID]['min'+insstat]; maxins = peinsize[sampleID]['max'+insstat]
		except KeyError: pass
		
		if not MAPSARECOMPLETE or REALIGN: pipeit('=> Aligning %s, min/max insert size = [%s,%s]'%(sampleID, minins, maxins),1)
		else: 
			pipeit('Read maps for %s are complete (%s)'%(sampleID, MAPSARECOMPLETE), 1)
			continue
		
		
		parmstr = '-t %s '%(qualstr); sinparmstr = '-t %s '%(qualstr)
		
		if ALIGNTOOL == 'bowtie2': 
			parmstr += '%s --minins %s --maxins %s '%(alignmentMode, minins, maxins)
			# parmstr += '--no-unal '; sinparmstr += '--no-unal '
			
			if alignmentPolicy == 'bestk':
				parmstr += '-k %s '%(maxnumalign)
				sinparmstr += '-k %s '%(maxnumalign)
			elif alignmentPolicy == 'bestm':
				parmstr += '-M %s '%(maxnumalign)
				sinparmstr += '-M %s '%(maxnumalign)
			elif alignmentPolicy == 'all':
				parmstr += '-a '
				sinparmstr += '-a '
			# elif alignmentPolicy == '':
			# this means default reporting of the best guess single alignment in bowtie2
			
		else:
			parmstr = '--pairtries 300 '
			sinparmstr = '%s '
			if readfileformat == 'fasta': parmstr += ' -f '; sinparmstr += ' -f '
			
			# output in sam format
			parmstr += ' --sam'; sinparmstr += ' --sam'
			
			# APPLY MISMATCHES
			# default mapping is v mode which is quality-blind
			mmstr = ' -v %s '%(mismatches['PE'])
			sinmmstr = ' -v %s '%(mismatches['SE'])
			if MAPWITHQUALITY:
				# use quality values for mapping
				# current is default seed length -l and max qual -e
				mmstr = ' -n %s '%(mismatches['PE'])
				sinmmstr = ' -n %s '%(mismatches['SE'])
		
			parmstr += mmstr
			sinparmstr += sinmmstr
			
			# APPLY STRATEGY
			if alignmentPolicy == 'basic':
				# for the maq-like policy
				parmstr += '--minins %s --maxins %s -m 1 --best'%(minins, maxins)
				sinparmstr += ' -m 1' # for singleton samples
		
			elif alignmentPolicy == 'all':
				# end-to-k policy (-v) with all valid alignments reported (-a)
				parmstr += '--minins %s --maxins %s -a'%(minins, maxins)
				sinparmstr += '-a' # for singleton samples
		
			elif alignmentPolicy == 'bestk':
				parmstr += '--minins %s --maxins %s -k %s --best'%(minins, maxins, maxnumalign)
				sinparmstr += '-k %s --best'%(maxnumalign) # for singleton samples
				
			elif alignmentPolicy == 'bestm':
				parmstr += '--minins %s --maxins %s -a -m %s --best'%(minins, maxins, maxnumalign)
				sinparmstr += '-a -m %s --best'%(maxnumalign) # for singleton samples
		
			elif alignmentPolicy == 'maq-bestm':
				parmstr += '--minins %s --maxins %s -a -m %s --best'%(minins, maxins, maxnumalign)
				sinparmstr += '-a -m %s --best'%(maxnumalign) # for singleton samples
		
			elif alignmentPolicy == 'bestm-stratum':
				parmstr += '--minins %s --maxins %s -a -m %s --best --strata'%(minins, maxins, maxnumalign)
				sinparmstr += '-a -m %s --best --strata'%(maxnumalign) # for singleton samples
		
			elif alignmentPolicy == 'best-stratum':
				parmstr += '--minins %s --maxins %s -a --best --strata'%(minins, maxins)
				sinparmstr += '-a --best --strata' # for singleton samples
		
		
		
		# trim reads
		if trim5: astr = ' --trim5 %s '%(trim5); parmstr += astr; sinparmstr += astr
		if trim3: astr = ' --trim3 %s '%(trim3); parmstr += astr; sinparmstr += astr
		
		# include arbitrary user specified arguments
		parmstr += userAlignArgs
		sinparmstr += userAlignArgs
		
		for run,lane in seqRuns[sampleID].keys():
			refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
			
			# command line overrides map file reference
			if bfaref: refGenome = bfaref
			
			# print 'entry', run, lane, 'het', hashet
			UID = '%s-%s-%s'%(sampleID,run,lane)
			if hashet and ALIGNTOOL != 'bowtie2': sys.exit('Heterozygous reference genome currently compatible with bowtie2 only. Use --tool bowtie2.')
			
			fastqFiles = rawfiles[sampleID]

			# determine if fasta format
			if len(fastqFiles) and ('.fasta' in fastqFiles[0]):
				parmstr += ' -f'
				sinparmstr += ' -f'
			
			
			if ALIGNTOOL == 'bowtie2':
				
				parts = refGenome.split('/')
				
				label = ''
				try: label = parts.pop()
				except IndexError: pass
				
				suffix = ''
				if '.fa' in label or '.fasta' in label:
					try: 
						more = label.split('.')
						# print 'more', more, 'label', label
					
						if len(more) > 1: 
							suffix = '.'+more.pop()
							label = '.'.join(more)
						else: label = '.'+more[0]
						# print 'A more', more, 'label', label
					except IndexError: pass
				
				refpath = ''
				if len(parts): refpath = '/'.join(parts)+'/'
				
				refGenomeBase = refpath+label
				# print 'base', refGenomeBase, refGenome.split('/')
				if hashet:
					refGenomeH = {}
					refGenomeH['h1'] = refpath+label+'.h1'
					refGenomeH['h2'] = refpath+label+'.h2'
				
				btidxsuffix = '.bt2'
				
				# generate index files if needed
				btindexfile = zdir+'bowtie.index.log'
				if not os.access(refGenomeBase+'.1'+btidxsuffix, os.F_OK):
					pipeit('- Generating reference genome index files...')
					refGenomePlus = refGenome
					if not hashet and suffix == '':
						if os.access(refGenome+'.fa', os.F_OK):
							refGenomePlus = refGenome+'.fa'
						elif os.access(refGenome+'.fasta', os.F_OK):
							refGenomePlus = refGenome+'.fasta'
						else:
							sys.exit('CANNOT ACCESS REFERENCE FILE.')
					
					call = '%sbowtie2-build "%s" "%s" >> "%s"'%(btdir,refGenomePlus,refGenomeBase,btindexfile)
					pipeit(call,2,readMapPipe)
					os.system(call)
					pipeit('done.',1)
				
				# heterozygous calls
				
				# partition reference into haploid components if desired
				if hashet:
					if not os.access(refGenomeH['h1']+suffix, os.F_OK):
						pipeit('- Partitioning diploid reference into 2 haploid files...')
						call = 'splitDiploidFasta.py -i "%s" -o "%s" >> "%s"'%(refGenome, refpath, readMapPipeFile)
						pipeit(call,1,readMapPipe)
						os.system(call)
						pipeit('done.',1)
						
					refGenome = refGenomeH
					
					if hashet and not (os.access(refGenome['h1']+'.1'+btidxsuffix, os.F_OK) or os.access(refGenome['h2']+'.1.ebwt', os.F_OK)):
						pipeit('- Generating haploid reference genome index files...')
						if not os.access(refGenome['h1']+'.1'+btidxsuffix, os.F_OK):
							call = '%sbowtie2-build "%s" "%s" >> "%s"'%(btdir,refGenome['h1']+suffix, refGenome['h1'], btindexfile)
							pipeit(call,2,readMapPipe)
							os.system(call)
							pipeit('done h1...')
						if not os.access(refGenome['h2']+'.1'+btidxsuffix, os.F_OK):
							call = '%sbowtie2-build "%s" "%s" >> "%s"'%(btdir,refGenome['h2']+suffix, refGenome['h2'], btindexfile)
							pipeit(call,2,readMapPipe)
							os.system(call)
							pipeit('done h2.',1)
				
				umm = not hashet and mismatches['SE'] or 0
				
				unm = "%s%s.unmapped.fq"%(unmapdir, UID)
				
				if paired:
					# Heterozygosity mapping (Diploid = 1)
					# Align the unmapped reads to each of 2 arbitrary haploid reference indices
					if hashet:
						Nceiling = '0,0' # allow 0 missing or heterozygous characters in valid alignments for first round
						# this guarantees SNP sites will not be covered
						
						call = '%sbowtie2 %s -N %s --n-ceil %s -p %s -x "%s" -1 "%s" -2 "%s" -S "%s" --un-conc "%s" 2>> "%s"'%(btdir,parmstr, umm, Nceiling, nsubprocs, refGenomeBase, fastqFiles[0], fastqFiles[1], btmapdir+UID+'.sam', unm, readMapPipeFile)
						pipeit('- Mapping paired end data for %s...'%(UID))
						pipeit(call,2,readMapPipe)
						if not DRYRUN and (not os.access(btmapdir+UID+'.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						
						# relax missing sequence
						Nceiling = '0,0.15' # default value
						
						unm1 = "%s%s.unmapped.1.fq"%(unmapdir, UID); unm2 = "%s%s.unmapped.2.fq"%(unmapdir, UID)
						
						call = '%sbowtie2 %s --no-head -N %s --n-ceil %s -p %s -x "%s" -1 "%s" -2 "%s" -S "%s" --un-conc "%s" 2>> "%s"'%(btdir,parmstr, mismatches['SE'], Nceiling, nsubprocs, refGenome['h1'], unm1, unm2, btmapdir+UID+'.h1.sam', unmapdir+UID+'.h1.unmapped.fq', readMapPipeFile)
						pipeit(call,2,readMapPipe)
						pipeit('- Mapping remaining paired end data to haplotype 1 for %s...'%(UID))
						printList([call], file=readMapPipeFile, ioMethod='a')
						if not DRYRUN and (not os.access(btmapdir+UID+'.h1.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						
						call = '%sbowtie2 %s --no-head -N %s --n-ceil %s -p %s -x "%s" -1 "%s" -2 "%s" -S "%s" --un-conc "%s" 2>> "%s"'%(btdir,parmstr, mismatches['SE'], Nceiling, nsubprocs, refGenome['h2'], unm1, unm2, btmapdir+UID+'.h2.sam', unmapdir+UID+'.h2.unmapped.fq', readMapPipeFile)
						pipeit(call,2,readMapPipe)
						pipeit('- Mapping remaining paired end data to haplotype 2 for %s...'%(UID))
						printList([call], file=readMapPipeFile, ioMethod='a')
						if not DRYRUN and (not os.access(btmapdir+UID+'.h2.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						
						# remove the partial files
						os.system('rm -f "%s"'%(unm1))
						os.system('rm -f "%s"'%(unm2))
					else: # (Haploid reference genome, paired end reads)
						# relax missing sequence
						Nceiling = '0,0.15' # default value
						# first align reads to reference genome as paired end
						call = '%sbowtie2 %s -N %s --n-ceil %s -p %s -x "%s" -1 "%s" -2 "%s" -S "%s" --un-conc "%s" 2>> "%s"'%(btdir,parmstr, umm, Nceiling, nsubprocs, refGenomeBase, fastqFiles[0], fastqFiles[1], btmapdir+UID+'.sam', unmapdir+UID+'.unmapped.fq', readMapPipeFile)
						pipeit('- Mapping paired end data for %s...'%(UID))
						pipeit(call,2,readMapPipe)
						if not DRYRUN and (not os.access(btmapdir+UID+'.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						
						
						# take the unmapped paired reads and align using tophat to transcriptome
						unm1 = "%s%s.unmapped.1.fq"%(unmapdir, UID); unm2 = "%s%s.unmapped.2.fq"%(unmapdir, UID)
						unfun = unmapdir+UID+'.unmapped.12.fq'
						if ALIGNTOTRANSCRIPTOME:
							call = '%stophat2 -g %s %s --no-convert-bam --GTF "%s" -p %s -o "%s" "%s" "%s" "%s" 2>> "%s"'%(thdir, maxnumalign, TOPHATFLAGS, gfffile, nsubprocs, btmapdir+UID+'.tophat/', refGenomeBase, unm1, unm2, zdir+'read_mapping.log')
							pipeit('- Mapping remaining paired end reads to transcriptome for %s...'%(UID))
							pipeit(call,2,readMapPipe)
							if not DRYRUN and (not os.access(btmapdir+UID+'.tophat/accepted_hits.sam', os.F_OK) or REALIGN): os.system(call)
							# cat the tophat log output to master log
							os.system('cat "%s" >> "%s"'%(btmapdir+UID+'.tophat/align_summary.txt', readMapPipeFile))
							pipeit('done.',1)
							
							# sys.exit('HOLD FOR NOW - DEBUG')
							# tophat returns a single unmapped reads file (unmapped.bam); convert this to fastq
							# unm = unmapdir+UID+'.tophat.unmapped.fq'
							call = 'bamToFastq -i "%s" -fq "%s"'%(btmapdir+UID+'.tophat/unmapped.bam', unfun)
							pipeit('- Converting tophat unmapped file to fastq for %s...'%(UID))
							pipeit(call,2,readMapPipe)
							if not DRYRUN and (not os.access(unfun, os.F_OK) or REALIGN): os.system(call)
							pipeit('done.',1)
						else:
							# no tophat; merge the 2 unmapped bowtie files and push through bowtie as singletons
							# unfun = "%s%s.unmapped.12.fq"%(unmapdir, UID)
							pipeit('- Combining _1 and _2 unmapped reads %s...'%(UID))
							if not DRYRUN and os.access(unm2, os.F_OK) or REALIGN: 
								# os.system('cat "%s" "%s" >> "%s"'%(unm1,unm2,unm))
								os.system('cat "%s" >> "%s"'%(unm2,unm1))
								os.system('mv "%s" "%s"'%(unm1,unfun))
								# remove the partial files
								os.system('rm -f "%s"'%(unm2))
							pipeit('done.', 1)
						
						# FINISH THIS PART HERE - MAKE SURE BOTH RNA and regular are covered
						# run tophat unmapped reads through bowtie once more
						# unm = unmapdir+UID+'.remain.unmapped.singleton.fq'
						call = '%sbowtie2 %s -N %s --n-ceil %s -p %s %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir, sinparmstr, umm, Nceiling, nsubprocs, BOWTIEFLAGS, refGenomeBase, unfun, btmapdir+UID+'.singleton.sam', unm, readMapPipeFile)

						pipeit('- Mapping remaining reads as singletons for %s...'%(UID))
						pipeit(call,2,readMapPipe)
						if not DRYRUN and (not os.access(btmapdir+UID+'.singleton.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						

						# remove the pre-unmapped file
						# os.system('rm -f "%s"'%(unfun))
						
						# iterative read trimming
						for utrim5,utrim3 in ITER_TRIM:
							punm = unm
							unm = unmapdir+UID+'.unmapped.%s-%s.fq'%(utrim5,utrim3)
							ftrim = btmapdir+UID+'.%s-%s.sam'%(utrim5,utrim3)
							sinparmstr2 = sinparmstr+' --trim5 %s --trim3 %s'%(utrim5,utrim3)
							call = '%sbowtie2 %s -N %s --n-ceil %s -p %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir, sinparmstr2, umm, Nceiling, nsubprocs, refGenomeBase, punm, ftrim, unm, readMapPipeFile)
							pipeit('- Mapping with trimming %s-%s for %s...'%(utrim5, utrim3, UID))
							pipeit(call,2,readMapPipe)
							if not DRYRUN and (not os.access(ftrim, os.F_OK) or REALIGN): os.system(call)
							if os.access(punm, os.F_OK): os.system('rm "%s"'%(punm))
							pipeit('done.',1)
						
				else: # (SE reads only)
					if len(fastqFiles) > 2:
						print >> sys.stderr, 'Warning: Found more than 2 matching read files. Using default s_lane_x_sequence.txt format.'
						fastqFiles = [path+'s_%s_1_sequence.txt'%(lane), path+'s_%s_2_sequence.txt'%(lane)]
					elif len(fastqFiles) != 1:
						print >> sys.stderr, 'Check read file names or path. Did not find 1 file for sample %s in dir "%s": %s'%(lane, path, fastqFiles)
					
					# Heterozygosity mapping
					# Align the unmapped reads to each of 2 arbitrary haploid reference indices
					if hashet:
						Nceiling = '0,0' # allow 0 missing or heterozygous characters in valid alignments for first round
						# this guarantees SNP sites will not be covered
						
						call = '%sbowtie2 %s -N %s --n-ceil %s -p %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir, sinparmstr, umm, Nceiling, nsubprocs, refGenomeBase, fastqFiles[0], btmapdir+UID+'.sam', unm, readMapPipeFile)
						pipeit('- Mapping singleton data for %s...'%(UID))
						pipeit(call,2,readMapPipe)
						if not DRYRUN and (not os.access(btmapdir+UID+'.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						
						# relax missing sequence
						Nceiling = '0,0.15' # default value
						
						unm12 = "%s%s.unmapped.fq"%(unmapdir, UID)
						unm = {}
						
						for hh in ['h1','h2']:
							unm[hh] = unmapdir+UID+'.%s.unmapped.fq'%(hh)
							
							call = '%sbowtie2 %s --no-head -N %s --n-ceil %s -p %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir,sinparmstr, mismatches['SE'], Nceiling, nsubprocs, refGenome[hh], unm12, btmapdir+UID+'.%s.sam'%(hh), unm[hh], readMapPipeFile)
							pipeit(call,2,readMapPipe)
							pipeit('- Mapping remaining singleton data to haplotype %s for %s...'%(hh, UID))
							printList([call], file=readMapPipeFile, ioMethod='a')
							if not DRYRUN and (not os.access(btmapdir+UID+'.%s.sam'%(hh), os.F_OK) or REALIGN): os.system(call)
							pipeit('done.',1)
						
						# map to transcriptome
						if ALIGNTOTRANSCRIPTOME:
							for hh in ['h1','h2']:
								# add suffix to reference
								reftmp = refGenome[hh]+'.fasta'
								if not os.access(reftmp, os.F_OK):
									if os.access(refGenome[hh]+'.fa', os.F_OK):
										reftmp = refGenome[hh]+'.fa'
									else:
										sys.exit('NO REFERENCE FOUND: %s'%(refGenome[hh]))
									
								call = '%stophat2 -g %s %s --no-convert-bam --GTF "%s" -p %s -o "%s" "%s" "%s" 2>> "%s"'%(thdir, maxnumalign, TOPHATFLAGS, gfffile, nsubprocs, btmapdir+UID+'.tophat.%s/'%(hh), refGenome[hh], unm[hh], zdir+'read_mapping.log')
								pipeit('- Mapping %s single end reads to transcriptome for %s...'%(hh, UID))
								pipeit(call,2,readMapPipe)
								if not DRYRUN and (not os.access(btmapdir+UID+'.tophat.%s/accepted_hits.sam'%(hh), os.F_OK) or REALIGN): 
									os.system(call)
									# cat the tophat log output to master log
									os.system('cat "%s" >> "%s"'%(btmapdir+UID+'.tophat.%s/align_summary.txt'%(hh), readMapPipeFile))
								pipeit('done.',1)
								
								# convert unmapped.bam to fastq file
								unm[hh] = unmapdir+UID+'.tophat.%s.unmapped.fq'%(hh)
								call = 'bamToFastq -i "%s" -fq "%s"'%(btmapdir+UID+'.tophat.%s/unmapped.bam'%(hh), unm[hh])
								pipeit('- Converting tophat.%s unmapped file to fastq for %s...'%(hh,UID))
								pipeit(call,2,readMapPipe)
								if not DRYRUN and (not os.access(unm[hh], os.F_OK) or REALIGN): os.system(call)
								pipeit('done.',1)
						
						# iterative read trimming
						punm = {}
						for utrim5,utrim3 in ITER_TRIM:
							for hh in ['h1','h2']:
								punm[hh] = unm[hh]
								unm[hh] = unmapdir+UID+'.unmapped.%s.%s-%s.fq'%(hh,utrim5,utrim3)
								ftrim = btmapdir+UID+'.%s.%s-%s.sam'%(hh,utrim5,utrim3)
								sinparmstr2 = sinparmstr+' --trim5 %s --trim3 %s'%(utrim5,utrim3)
								call = '%sbowtie2 %s -N %s --n-ceil %s -p %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir, sinparmstr2, umm, Nceiling, nsubprocs, refGenome[hh], punm[hh], ftrim, unm[hh], readMapPipeFile)
								pipeit('- Mapping %s with trimming %s-%s for %s...'%(hh, utrim5, utrim3, UID))
								pipeit(call,2,readMapPipe)
								if not DRYRUN and (not os.access(ftrim, os.F_OK) or REALIGN): os.system(call)
								if os.access(punm[hh], os.F_OK): os.system('rm "%s"'%(punm[hh]))
								pipeit('done.',1)
						
						if os.access(unm12, os.F_OK): os.system('rm "%s"'%(unm12))
						
					else: # (Haploid reference genome)
						# relax missing sequence
						Nceiling = '0,0.15' # default value
						
						call = '%sbowtie2 %s -N %s --n-ceil %s -p %s %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir, sinparmstr, umm, Nceiling, nsubprocs, BOWTIEFLAGS, refGenomeBase, fastqFiles[0], btmapdir+UID+'.sam', unm, readMapPipeFile)
						pipeit('- Mapping single end reads for %s...'%(UID))
						pipeit(call,2,readMapPipe)
						if not DRYRUN and (not os.access(btmapdir+UID+'.sam', os.F_OK) or REALIGN): os.system(call)
						pipeit('done.',1)
						
						# map to transcriptome
						if ALIGNTOTRANSCRIPTOME:
							call = '%stophat2 -g %s %s --no-convert-bam --GTF "%s" -p %s -o "%s" "%s" "%s" 2>> "%s"'%(thdir, maxnumalign, TOPHATFLAGS, gfffile, nsubprocs, btmapdir+UID+'.tophat/', refGenomeBase, unm, zdir+'read_mapping.log')
							pipeit('- Mapping single end reads to transcriptome for %s...'%(UID))
							pipeit(call,2,readMapPipe)
							if not DRYRUN and (not os.access(btmapdir+UID+'.tophat/accepted_hits.sam', os.F_OK) or REALIGN): os.system(call)
							# cat the tophat log output to master log
							os.system('cat "%s" >> "%s"'%(btmapdir+UID+'.tophat/align_summary.txt', readMapPipeFile))
							pipeit('done.',1)
							
							# convert unmapped.bam to fastq file
							unm = unmapdir+UID+'.tophat.unmapped.fq'
							call = 'bamToFastq -i "%s" -fq "%s"'%(btmapdir+UID+'.tophat/unmapped.bam', unm)
							pipeit('- Converting tophat unmapped file to fastq for %s...'%(UID))
							pipeit(call,2,readMapPipe)
							if not DRYRUN and (not os.access(unm, os.F_OK) or REALIGN): os.system(call)
							pipeit('done.',1)
							
						# iterative read trimming
						# if ITER_TRIM:
						for utrim5,utrim3 in ITER_TRIM:
							punm = unm
							unm = unmapdir+UID+'.unmapped.%s-%s.fq'%(utrim5,utrim3)
							ftrim = btmapdir+UID+'.%s-%s.sam'%(utrim5,utrim3)
							sinparmstr2 = sinparmstr+' --trim5 %s --trim3 %s'%(utrim5,utrim3)
							call = '%sbowtie2 %s -N %s --n-ceil %s -p %s -x "%s" -U "%s" -S "%s" --un "%s" 2>> "%s"'%(btdir, sinparmstr2, umm, Nceiling, nsubprocs, refGenomeBase, punm, ftrim, unm, readMapPipeFile)
							pipeit('- Mapping with trimming %s-%s for %s...'%(utrim5, utrim3, UID))
							pipeit(call,2,readMapPipe)
							if not DRYRUN and (not os.access(ftrim, os.F_OK) or REALIGN): os.system(call)
							if os.access(punm, os.F_OK): os.system('rm "%s"'%(punm))
							pipeit('done.',1)
						# else:
						#	move the tophat unmapped file to unmapped directory
						#
					
			# bowtie/tophat alignment tools
			elif paired:
				# bowtie 1 for alignment
				
				if ALIGNTOOL == 'bowtie':
					call = '%sbowtie %s -p %s "%s" -1 "%s" -2 "%s" "%s%s.sam" --un "%s%s_unfq.fq" 2>> "%s"'%(btdir,parmstr, nsubprocs, refGenome, fastqFiles[0], fastqFiles[1], btmapdir, UID, btmapdir, UID, zdir+'read_mapping.log')
				
				elif ALIGNTOOL == 'tophat':
					createdir(btmapdir+'tophat-%s'%(UID))
					call = '%stophat -N %s %s -g %s -p %s -o "%s" "%s" "%s","%s" 2>> "%s"'%(thdir, mismatches['PE'], TOPHATFLAGS, maxnumalign, nsubprocs, btmapdir+'tophat-%s'%(UID), refGenome, fastqFiles[0],fastqFiles[1], zdir+'read_mapping.log')
				
				elif ALIGNTOOL == 'tophat2':
					# check if user specified fragment size
					fragstr = ''# '-r 50'
					if diamidx: 
						# print 'index', diamidx, seqRuns[sampleID][(run,lane)]
						fraglen = int(seqRuns[sampleID][(run,lane)][-1])
						# need to fix this
						# fragstr = '-r %s'%(fraglen)
						
					createdir(btmapdir+'tophat-%s'%(UID))
					if len(fastqFiles) == 0: sys.exit("no input files found. Exiting.")
					call = '%stophat2 -g %s %s -p %s -o "%s" "%s" "%s" "%s" %s 2>> "%s"'%(thdir, maxnumalign, TOPHATFLAGS, nsubprocs, btmapdir+'tophat-%s'%(UID), refGenomeBase, fastqFiles[0], fastqFiles[1], fragstr, zdir+'read_mapping.log')
				
				
				pipeit('- Mapping paired end data for %s...'%(UID),1)
				printList([call], file=zdir+'read_mapping.log', ioMethod='a')
				
				if not DRYRUN and (not os.access(btmapdir+UID+'.sam', os.F_OK) or REALIGN): os.system(call)
				
				# set this up to process the remaining singletons as well, and then merge the output files...
				# outfile is _unfq.fq
				
				if ALIGNTOOL == 'bowtie': # map the remaining singleton reads
					pipeit('- Mapping left-over singleton data for %s...'%(UID),1)
					for lane in [1,2]:
						call2 = '%sbowtie %s -p %s "%s" "%s" "%s" --un "%s" 2>> "%s"'%(btdir, parmstr, nsubprocs, refGenome, btmapdir+UID+'_unfq_%s.fq'%(lane), btmapdir+UID+'_singles_%s.sam'%(lane), unmapdir+UID+'_unmapped_%s.fq'%(lane), zdir+'read_mapping.log')
						printList([call2], file=zdir+'read_mapping.log', ioMethod='a')
						if not DRYRUN and ((not MAPSARECOMPLETE and not os.access(btmapdir+UID+'_singles_%s.sam'%(lane), os.F_OK)) or REALIGN): os.system(call2)
						
						# remove the input fq file, contents of which are in singles_map and unmapped
						os.system('rm -f "%s"'%(btmapdir+UID+'_unfq_%s.fq'%(lane)))
						
				elif ALIGNTOOL == 'bowtie2':
					# rename the --un files from the PE run as unmapped and move them into the unmapped folder
					for lane in [1,2]:
						call = 'mv %s %s'%(btmapdir+UID+'_unfq_%s.fq'%(lane), unmapdir+UID+'.unmapped_%s.fq'%(lane))
						os.system(call)
			else:
				if not COLORSPACE:
					if len(fastqFiles) > 2:
						print >> sys.stderr, 'Warning: Found more than 2 matching read files. Using default s_lane_x_sequence.txt format.'
						fastqFiles = [path+'s_%s_1_sequence.txt'%(lane), path+'s_%s_2_sequence.txt'%(lane)]
					elif len(fastqFiles) != 1:
						print >> sys.stderr, 'Check read file names or path. Did not find 1 file for sample %s in dir "%s": %s'%(lane, path, fastqFiles)
					
				if ALIGNTOOL == 'bowtie':
					call = '%sbowtie %s -p %s "%s" "%s" "%s%s.sam" --un "%s%s_unmapped.fq" 2>> "%s"'%(btdir, sinparmstr, nsubprocs, refGenome, fastqFiles[0], btmapdir, UID, unmapdir, UID, zdir+'read_mapping.log')
				elif ALIGNTOOL == 'tophat':
					createdir(btmapdir+'tophat-%s'%(UID))
					# datstr = ' '.join(['"%s"'%(fqf) for fqf in fastqFiles])
						
					call = '%stophat -N %s %s -g %s -p %s -o "%s" "%s" "%s" 2>> "%s"'%(thdir, mismatches['SE'], TOPHATFLAGS, maxnumalign, nsubprocs, btmapdir+'tophat-%s'%(UID), refGenome, fastqFiles[0], zdir+'read_mapping.log')
				elif ALIGNTOOL == 'tophat2':
					createdir(btmapdir+'tophat-%s'%(UID))
					# datstr = ' '.join(['"%s"'%(fqf) for fqf in fastqFiles])
					call = '%stophat2 -g %s %s -p %s -o "%s" "%s" "%s" 2>> "%s"'%(thdir, maxnumalign, TOPHATFLAGS, nsubprocs, btmapdir+'tophat-%s'%(UID), refGenomeBase, fastqFiles[0], zdir+'read_mapping.log')
					
				if COLORSPACE:
					
					if ALIGNTOOL == 'tophat':
						sys.exit('tophat mapping in colorspace not yet implemented.')
					
					# trim any .txt files
					fastqFiles = filter(lambda x: x.count('.csfasta') or x.count('.qual'), fastqFiles)
						
					if len(fastqFiles) != 2:
						print >> sys.stderr, 'Check read file names or path. Did not find 2 files (.csfasta and .qual) for sample %s in file "%s": %s'%(lane,path, fastqFiles)
							
					# call = '%sbowtie %s -C -p %s "%s" -f "%ss_%s_sequence.csfasta" -Q "%ss_%s_sequence.qual" "%s%s.map" --un "%s%s_unmapped.fq" --snpfrac %s'%(btdir, sinparmstr, nprocs, refGenome, path, lane, path,lane, btmapdir, UID, unmapdir, UID, colorspacetheta)
					call = '%sbowtie %s -C -p %s "%s" -f "%s" -Q "%s" "%s%s.sam" --un "%s%s_unmapped.fq" --snpfrac %s 2>> "%s"'%(btdir, sinparmstr, nsubprocs, refGenome, fastqFiles[0], fastqFiles[1], btmapdir, UID, unmapdir, UID, colorspacetheta, zdir+'read_mapping.log')
						
						
				if not DRYRUN and ((not MAPSARECOMPLETE and (not os.access(btmapdir+UID+'.sam', os.F_OK))) or REALIGN):
					if 'tophat' in ALIGNTOOL and os.access(btmapdir+'tophat-%s/accepted_hits.bam'%(UID), os.F_OK): continue
					pipeit('- Assembling original singleton data for %s...'%(UID),1)
					printList([call], file=zdir+'read_mapping.log', ioMethod='a')
					os.system(call)
					
				
		
		
		
			
		print '- DONE', sampleID
	
	print '=> Done read mapping for all specified samples.'
	print
	if not ORGANIZEREADS and not LIKCALLING and not PEAKCALLING and not JOINPIECES and not ALIGNTOTRANSCRIPTOME and not RUNCUFFLINKS:
		sys.exit('=> Done for now. Use the -O and -index flags to compute organize and index read maps for downstream analysis.')
	

# =========================================================
# Part IB: Assay raw read coverage / Insert size stats
# =========================================================
if COUNTREADS:
	print '- Loading reference genome...'
	G = readFasta(bfaref)
	# chromosome total lengths (bp)
	chromlengths = dict( map(lambda k: (k, len(G[k]['seq'])), G) )
	master_chroms = chromlengths.keys()
	
	readbphead = ['SampleID', 'Filename', 'Lane', 'PE', 'Reads', 'Bases', 'Coverage']
	readbptab = []
	readlen = 40
	print 'Computing number of reads for samples...'
	for sampleID in sortseqruns:
		for run,lane in seqRuns[sampleID].keys():
			refGenome,paired,path = seqRuns[sampleID][(run,lane)][:3]
			# count up reads and bases
			if paired:
				for pe in [1]: # don't need the second because it has same number of reads
					call = "wc -l %s > %s"%(path+'s_'+lane+'_'+str(pe)+'_sequence.txt', outdir+'tmpread.txt')
					print 'call', call
					os.system(call)
					line = readList(outdir+'tmpread.txt')
					nreads, crap = line[0].split(' ')
					nreads = float(nreads)/2.0 # each read takes up 2 rows in file
					row = [sampleID, run, lane, pe,'%.5e'%(float(nreads)), '%.5e'%(float(nreads)*readlen), '%.2f'%(float(nreads)*readlen/float(2*sum(chromlengths.values())))]
					readbptab += [row]
				os.system('rm "%s"'%(outdir+'tmpread.txt'))
	printFormattedTable(readbptab, readbphead, file=zdir+'Raw reads summary-%s.txt'%(zkey))
	
	sys.exit('Done counting raw reads. Exiting...')
	
	# remaining is deprecated - organize reads returns mapped totals
	
	mapreadshead = ['SampleID', 'Filename', 'Lane', 'PE Reads', 'Rem Singletons', 'PE Bases', 'Rem Singleton Bases', 'PE Coverage', 'Rem Sing. Cov.']
	mapreadstab = []
	
	# count up number of mapped reads and bases
	# -----------------------------------------
	for sampleID in sortseqruns:
		for run,lane in seqRuns[sampleID].keys():
			refGenome,paired,path = seqRuns[sampleID][(run,lane)][:3]
			UID = '%s-%s-%s'%(sampleID,run,lane)
			
			call = "wc -l %s > %s"%(btmapdir+UID+'.sam', outdir+'tmpread.txt')
			print 'call', call
			os.system(call)
			line = readList(outdir+'tmpread.txt')
			nreads, crap = line[0].split(' ')
			
			# this is still a work in progress
			snreads = 0
			if paired and USESINGLETONS:
				for pe in ['1','2']:
					call = "wc -l %s > %s"%(btmapdir+UID+'_singles_'+pe+'.sam', outdir+'tmpread.txt')
					print 'call', call
					os.system(call)
					line = readList(outdir+'tmpread.txt')
					tmpnreads, crap = line[0].split(' ')
					snreads += int(tmpnreads)
			pecov = '%.2f'%(float(nreads)*readlen/float(2*sum(chromlengths.values())))
			sincov = '%.2f'%(float(snreads)*readlen/float(2*sum(chromlengths.values())))
			row = [sampleID, run, lane, '%.2e'%(float(nreads)), '%.2e'%(float(snreads)),  '%.2e'%(float(nreads)*readlen),'%.2e'%(float(snreads)*readlen), pecov, sincov]
			mapreadstab += [row]
			
	os.system('rm "%s"'%(outdir+'tmpread.txt'))
	
	printFormattedTable(mapreadstab, mapreadshead, file=zdir+'Mapped reads summary-%s.txt'%(zkey))
	sys.exit('Done counting raw and mapped reads. Exiting...')
#
if GETINSERTSIZE:
	# save an insert size histogram for each sampleID and a summary table providing .95, .99, and min/max stats
	# this only works in bowtie saved output as sam maps
	
	print 'Computing PE insert size distributions...'
	peinsize = []
	
	for sampleID in sortseqruns:
		insvalues = []
		
		for run,lane in seqRuns[sampleID].keys():
			refGenome,paired,path = seqRuns[sampleID][(run,lane)][:3]
			UID = '%s-%s-%s'%(sampleID,run,lane)
			
			# only used for paired end reads
			if paired:
				# get data from the sam/bam file
				fh = open(btmapdir+UID+'.sam')
				for line in fh.readlines():
					if line[0]=='@': continue
					
					row = line.strip('\n').split('\t')
					# print 'row has', len(row)
					# print row
					insvalue = abs(int(row[8]))
					if insvalue != 0:
						insvalues += [insvalue]
				fh.close()
		
		insvalues = filter(lambda x: x!=nastr and True, insvalues)
		
		print '- %s has %s values'%(sampleID, len(insvalues))
		print '  saving histogram...'
		hist = histogram(insvalues, nbins=25, yfreq=True)
		printList(hist, file=zinsdir+sampleID+' histogram.txt')
		
		print '  calculating statistics...'
		ansd = sd(insvalues)
		meanins = avg(insvalues)
		sdL = meanins - ansd; sdR = meanins + ansd
		
		ninefiveL = percentile(insvalues, .025)
		ninefiveR = percentile(insvalues, .975)
		
		ninenineL = percentile(insvalues, .005)
		ninenineR = percentile(insvalues, .995)
		
		mn = min(insvalues); mx = max(insvalues)
		
		perow = [sampleID, meanins, sdL, sdR, ninefiveL, ninefiveR, ninenineL, ninenineR, mn, mx]
		
		print perow
		
		peinsize += [perow]
	
	header = ['ID', 'Mean', 'minSD', 'maxSD', 'min.95', 'max.95', 'min.99', 'max.99', 'min', 'max']
	printTable(peinsize, header, file=peinsfile)
	
	print '=> DONE.\nRerun glimmrHet.py using the -a or -ra flags to utilize this information.\n'
	sys.exit()


# ==========================================================
# Part II: generate final read files grouped by multiplicity
# ==========================================================
crdfile = zdir+'Chromosome_read_distribution.txt'
if ORGANIZEREADS:
	pipeit('READ ORGANIZATION',1)
	
	if not EXORGANIZE: 
		d,r,c = readTable(zdir+'readFilesToSampleIDs.txt',header=0)
		rawfiles = dict([(r[i],map(lambda x: x.strip(), d[i])) for i in range(len(d))])
	
	# touch new summary tables
	pmafile = zdir+'QC_alignment_mismatches%s.txt'%(ubamap)
	if not os.access(pmafile, os.F_OK):
		header=['Sample', 'Files', 'Diploid', 'Total reads', 'Mapped reads', 'Proportion reads mapped', 'Mapped alignments', 'Perfect alignments', 'Mismatched alignments', 'Total alignments', 'Frac perfect', 'Frac mismatch', 'Avg mismatch/aln', 'Alns w/mismatch', 'Aln multiplicity', 'Deg alns/read']
		record('\t'.join(header), pmafile, 'w')
	
	almfile = zdir+'QC_alignment_multiplicity%s.txt'%(ubamap)
	if not os.access(almfile, os.F_OK):
		record('Sample: SE and PE multiplicity distributions (alignments, reads)', almfile, 'w')
	
	if not os.access(crdfile, os.F_OK):
		record('SampleID\tDiploid\tChrom\tUnique SE\tUnique PE\tNon-unique SE\tNon-unique PE', crdfile, 'w')
	
	currentsample = 0
	for sampleID in sorted(seqRuns.keys()):
		hashet = None
		for run,lane in seqRuns[sampleID].keys():
			refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
		hapdip = 'hap'
		if hashet: hapdip = 'dip'
		
		sys.stdout.write('=> Sample %s\n'%(sampleID)); sys.stdout.flush()
		
		if not REORGANIZEREADS:
			if not PESECOMBO and os.access(btmapdir+sampleID+'.PE.unique.map', os.F_OK):
				sys.stdout.write('  Already processed.\n'); sys.stdout.flush(); continue
			elif PESECOMBO and os.access(btmapdir+sampleID+'.unique.map', os.F_OK):
				sys.stdout.write('  Already processed.\n'); sys.stdout.flush(); continue
			
		
		# record read multiplicity - number of times each read appears
		# ------------------------------------------------------------
		# pipeit('- Assessing read multiplicity...',1)
		
		# per sampleID variables of interest
		npmms = {hapdip:0}; nrows = {hapdip:0} # count up total rows and total rows with mismatches
		nqmms = {hapdip:0} # average number of mismatches per alignment
		hqmms = {hapdip:{}} # histogram of number of alignments x number of mismatches
		# assco = {hapdip:0}; nrows = {hapdip:0} # average bowtie2 alignment score
		totreads = {hapdip:{'SE':0, 'PE':0, 'Total':0}}
		totaln = {hapdip:{'SE':0, 'PE':0, 'Total':0}}
		sidHist = {hapdip:{'SE':{}, 'PE':{}, 'Total':{}}} # aggregate multiplicity histogram
		
		ulinesat = {hapdip:{'PE':{}, 'SE':{}}}
		dlinesat = {hapdip:{'PE':{}, 'SE':{}}}
		isizes = {hapdip:[]}
		
		# for each sampleID save new SAM-format read files
		funi = {hapdip:{}, hapdip:{}}
		fdeg = {hapdip:{}, hapdip:{}}
		
		
		if ALLELESPECIFIC:
			funi = {'h1':{}, 'h2':{}}
			fdeg = {'h1':{}, 'h2':{}}
			ulinesat = {'h1':{'PE':{}, 'SE':{}}, 'h2':{'PE':{}, 'SE':{}}}
			dlinesat = {'h1':{'PE':{}, 'SE':{}}, 'h2':{'PE':{}, 'SE':{}}}
			npmms = {'h1':0, 'h2':0}; nrows = {'h1':0, 'h2':0} # count up total rows and total rows with mismatches
			nqmms = {'h1':0, 'h2':0} # average number of mismatches per alignment
			hqmms = {'h1':{}, 'h2':{}} # histogram of number of alignments x number of mismatches
			totreads = {'h1':{'SE':0, 'PE':0, 'Total':0}, 'h2':{'SE':0, 'PE':0, 'Total':0}}
			totaln = {'h1':{'SE':0, 'PE':0, 'Total':0}, 'h2':{'SE':0, 'PE':0, 'Total':0}}
			sidHist = {'h1':{'SE':{}, 'PE':{}, 'Total':{}}, 'h2':{'SE':{}, 'PE':{}, 'Total':{}}} # aggregate multiplicity histogram
			isizes = {'h1':[], 'h2':[]}
			
			stmp = open(btmapdir+sampleID+'.h1.unique.map', 'w')
			funi['h1'] = {'SE':stmp, 'PE':stmp}
			stmp = open(btmapdir+sampleID+'.h1.degenerate.map', 'w')
			fdeg['h1'] = {'SE':stmp, 'PE':stmp}
			
			stmp = open(btmapdir+sampleID+'.h2.unique.map', 'w')
			funi['h2'] = {'SE':stmp, 'PE':stmp}
			stmp = open(btmapdir+sampleID+'.h2.degenerate.map', 'w')
			fdeg['h2'] = {'SE':stmp, 'PE':stmp}
		else:
			if UDPESECOMBO:
				stmp = open(btmapdir+sampleID+'.map', 'w')
				funi[hapdip] = {'SE':stmp, 'PE':stmp}
				fdeg[hapdip] = {'SE':stmp, 'PE':stmp}
			elif PESECOMBO:
				stmp = open(btmapdir+sampleID+'.unique.map', 'w')
				funi[hapdip] = {'SE':stmp, 'PE':stmp}
				stmp = open(btmapdir+sampleID+'.degenerate.map', 'w')
				fdeg[hapdip] = {'SE':stmp, 'PE':stmp}
			else:
				# Partition reads into 2 multiplicity groups
				funi[hapdip] = {'SE':open(btmapdir+sampleID+'.SE.unique.map', 'w'), 'PE':open(btmapdir+sampleID+'.PE.unique.map', 'w')}
				fdeg[hapdip] = {'SE':open(btmapdir+sampleID+'.SE.degenerate.map', 'w'), 'PE':open(btmapdir+sampleID+'.PE.degenerate.map', 'w')}
		
		HEADER = {'SE':None,'PE':None}
		HEADERUNI = {'SE':"",'PE':""}
		HEADERDEG = {'SE':"",'PE':""}
		
 		# the list of (run,lane) pairs of samples corresponding to this sampleID
		thesamples = seqRuns[sampleID].keys(); thesamples.sort()
		
		# only need one sample just as a dummy here
		if EXORGANIZE: thesamples = [thesamples[0]]
		
		for run,lane in thesamples:
			# multiplicity counting by read
			refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
			UID = '%s-%s-%s'%(sampleID,run,lane)
			# if user wants to organize a preexisting map
			if EXORGANIZE: UID = sampleID
			
			finalbads = open(unmapdir+UID+'.badmapping.txt', 'w')
			
			rlp = '%s,%s-'%(run,lane) # run,lane prefix for all reads
			currentsample += 1
			
			# what map files are relevant for this sampleID?
			fpairs = []
			
			if not EXORGANIZE:
				pipeit('- [DIPLOID=%s, ALLELESPECIFIC=%s]'%(bool(hashet),ALLELESPECIFIC),1)
				if hashet and ALLELESPECIFIC:
					if ALIGNTOOL != 'bowtie2': sys.exit('allele specific mapping only compatible with bowtie2.')
					
					fpairs += [('h1','SE', btmapdir+UID+'.h1.sorted.sam')]
					fpairs += [('h2','SE', btmapdir+UID+'.h2.sorted.sam')]
					
					if not os.access(btmapdir+UID+'.h1.sorted.sam', os.F_OK) or not os.access(btmapdir+UID+'.h2.sorted.sam', os.F_OK):
						# bowtie has issue printing unaligned reads
						# merge and sort the three files...
						pipeit('- Merging haploid and diploid readmaps...')
						
						for hh in ['h1','h2']:
							hhtarget = btmapdir+UID+'.%s.merge.sam'%(hh)
							os.system('cat "%s" "%s" > "%s"'%(btmapdir+UID+'.%s.sam'%(hh), btmapdir+UID+'.sam', hhtarget))
							if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(btmapdir+UID+'.%s.sam'%(hh)))
							
							if ALIGNTOTRANSCRIPTOME:
								f = btmapdir+UID+'.tophat.%s/accepted_hits.sam'%(hh)
								call = 'cat "%s" >> "%s"'%(f, hhtarget)
								os.system(call)
								if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(f))
								
							for utrim5,utrim3 in ITER_TRIM:
								f = btmapdir+UID+'.%s.%s-%s.sam'%(hh,utrim5,utrim3)
								os.system('cat "%s" >> "%s"'%(f, hhtarget))
								if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(f))
								
						if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(btmapdir+UID+'.sam'))
						pipeit('done.',1)
						
						pipeit('- Sorting merged readmap...')
						os.system('sort -k 1,1 -S 50% --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, btmapdir+UID+'.h1.merge.sam', btmapdir+UID+'.h1.sorted.sam'))
						os.system('sort -k 1,1 -S 50% --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, btmapdir+UID+'.h2.merge.sam', btmapdir+UID+'.h2.sorted.sam'))
						
						if not REORGANIZEREADS and not KEEPBOWTIE: 
							os.system('rm -f "%s"'%(btmapdir+UID+'.h1.merge.sam'))
							os.system('rm -f "%s"'%(btmapdir+UID+'.h2.merge.sam'))
						pipeit('done.',1)
					else:
						pipeit('- Found existing files: %s, %s'%(btmapdir+UID+'.h1.sorted.sam', btmapdir+UID+'.h2.sorted.sam'),1)
					
				else:
					if ALIGNTOOL == 'bowtie2':
						fpairs += [(hapdip, 'SE', btmapdir+UID+'.sorted.sam')]
						
						if not os.access(btmapdir+UID+'.sorted.sam', os.F_OK):
							if hashet:
								# bowtie has issue printing unaligned reads
								# merge and sort the three files...
								if len(ITER_TRIM): pipeit('- Merging haploid and diploid and trimmed readmaps...')
								else: pipeit('- Merging haploid and diploid readmaps...')
								
								filestomerge = []
								filestomerge += [btmapdir+UID+'.h1.sam']
								filestomerge += [btmapdir+UID+'.h2.sam']
								
								# ALLELE SPECIFIC HERE
								if ALIGNTOTRANSCRIPTOME:
									# add transcriptome map
									filestomerge += [btmapdir+UID+'.tophat.%s/accepted_hits.sam'%('h1')]
									filestomerge += [btmapdir+UID+'.tophat.%s/accepted_hits.sam'%('h2')]
									
								for utrim5,utrim3 in ITER_TRIM:
									filestomerge += [btmapdir+UID+'.%s.%s-%s.sam'%('h1',utrim5,utrim3)]
									filestomerge += [btmapdir+UID+'.%s.%s-%s.sam'%('h2',utrim5,utrim3)]
									
								os.system('cp "%s" "%s"'%(btmapdir+UID+'.sam', btmapdir+UID+'.merge.sam'))
								
								for f in filestomerge:
									call = 'cat "%s" >> "%s"'%(f, btmapdir+UID+'.merge.sam')
									os.system(call)
									if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(f))
									
								# if not REORGANIZEREADS and not KEEPBOWTIE: 
								# 	os.system('rm -f "%s"'%(btmapdir+UID+'.sam'))
								# 	os.system('rm -f "%s"'%(btmapdir+UID+'.h1.sam'))
								# 	os.system('rm -f "%s"'%(btmapdir+UID+'.h2.sam'))
								pipeit('done.',1)
								
								pipeit('- Sorting merged readmap...')
								os.system('sort -k 1,1 -S 50%'+' --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, btmapdir+UID+'.merge.sam', btmapdir+UID+'.sorted.sam'))
								if not REORGANIZEREADS and not KEEPBOWTIE: 
									# os.system('rm -f "%s"'%(btmapdir+UID+'.sam'))
									os.system('rm -f "%s"'%(btmapdir+UID+'.merge.sam'))
								pipeit('done.',1)
							else:
								if len(ITER_TRIM): pipeit('- Merging haploid and diploid and trimmed readmaps...')
								# else: pipeit('- Merging haploid and diploid readmaps...')
								
								if not len(ITER_TRIM) and not ALIGNTOTRANSCRIPTOME:
									# haploid file is the one samfile we have; just rename
									os.system('mv "%s" "%s"'%(btmapdir+UID+'.sam',btmapdir+UID+'.sorted.sam'))
								else:
									os.system('cp "%s" "%s"'%(btmapdir+UID+'.sam',btmapdir+UID+'.merge.sam'))
									
									if paired:
										os.system('cat "%s" >> "%s"'%(btmapdir+UID+'.singleton.sam',btmapdir+UID+'.merge.sam'))
										if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(btmapdir+UID+'.singleton.sam'))
										
									if ALIGNTOTRANSCRIPTOME:
										# merge genome and transcriptome maps
										fff = btmapdir+UID+'.tophat/accepted_hits.sam'
										call = 'cat "%s" >> "%s"'%(fff, btmapdir+UID+'.merge.sam')
										os.system(call)
										if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(fff))
										
									for utrim5,utrim3 in ITER_TRIM:
										newf = btmapdir+UID+'.%s-%s.sam'%(utrim5,utrim3)
										os.system('cat "%s" >> "%s"'%(btmapdir+UID+'.%s-%s.sam'%(utrim5,utrim3),btmapdir+UID+'.merge.sam'))
										if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(newf))
										
									pipeit('- Sorting merged readmap...')
									os.system('sort -k 1,1 -S 50%'+' --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, btmapdir+UID+'.merge.sam', btmapdir+UID+'.sorted.sam'))
									if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(btmapdir+UID+'.merge.sam'))
									pipeit('done.',1)
									
						else:
							pipeit('- Found existing file: %s'%(btmapdir+UID+'.sorted.sam'),1)
					else:
						if paired == 1:
							# note the PE goes last so that the SAM header is correctly saved to the map files
							# filenames for singleton leftovers from PE
							fpairs += [(hapdip,'SE', btmapdir+UID+'_singles_1.sam')]
							fpairs += [(hapdip,'SE', btmapdir+UID+'_singles_2.sam')]
							# filename for standard PE run
							fpairs += [(hapdip,'PE', btmapdir+UID+'.sam')]
						else:
							# filename for singleton run
							fpairs += [(hapdip,'SE', btmapdir+UID+'.sam')]
			
			else:
				qfile = btmapdir+UID+'.sam'
				qsfile = btmapdir+UID+'.sorted.sam'
				if EXPATH:
					qfile = EXPATH+UID+'.sam'
					qsfile = EXPATH+UID+'.sorted.sam'
				
				fpairs = [(hapdip, 'SE', qsfile)]
				
				if not EXOSORT: qfile = qsfile
				
				pipeit('- Looking for file: %s...'%(qfile))
				if not os.access(qfile, os.F_OK): sys.exit('Cannot find file!')
				else: pipeit('FOUND.',1)
				
				if EXOSORT:
					pipeit('- Sorting external readmap...')
					os.system('sort -k 1,1 -S 50%'+' --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, qfile, qsfile))
					if not REORGANIZEREADS and not KEEPBOWTIE: os.system('rm -f "%s"'%(qfile))
					pipeit('done.',1)
				else:
					pipeit('- User provided pre-sorted map file.', 1)
			
			# special nounique map to store discarded best-noguess alignments
			if STRATEGY=='bestno': fbestno = open(unmapdir+UID+'.bestno_discard.sam', 'w')
			
			# for each map file, partition into organized map files
			for hapdip,pese,fname in fpairs:
				# Stream-organization method for read organization
				# because the reads are sorted by id, we can easily count number of alignments
				# based on number of rows
				
				keyid = None # id of read under consideration
				keypese = None # if this read maps as SE or PE
				keypeseint = 0 # 0 = SE, 1 = PE
				ALN = {} # store alignment rows for a single hitid
				# MMS = [] # store mismatches per alignment row
				# UNIQUEPARTALN = [] # store alignment rows - but only readid, pos, ...quality
				
				# histogram summarizing multiplicity distribution
				RIDm = dict(map(lambda x: ((0,x),0), xrange(1,maxnumalign+1)))
				RIDm.update(dict(map(lambda x: ((1,x),0), xrange(1,maxnumalign+1))))
				
				if os.access(fname, os.F_OK):
					parts = fname.split('/')
					pipeit('- Sample %s: Counting alignments for %s...'%(currentsample, parts[len(parts)-1]))
					fh = open(fname)
					count = 0
					
					# FETCHING THE HEADER MAY HAVE A BUG HERE
					for row in fh:
						if count % 5000000 == 0: pipeit('%.1e,'%(count))
						count += 1
						if count >= 1000000: break # is this some kind of failsafe?
						# reprint sam-format headers, but only once
						if (HEADER[pese]==None or HEADER[pese]==fname) and (row[0:3]=='@PG' or row[0:3]=='@HD' or row[0:3]=='@SQ'):
							# don't print header twice
							if UDPESECOMBO and pese == 'PE':
								# print>>funi[pese], row[:-1]
								HEADER[pese] = fname
								HEADERUNI[pese] += row#[:-1]
							elif PESECOMBO:
								if pese == 'PE':
									HEADER[pese] = fname
								else:
									# print>>funi[pese], row[:-1]
									# print>>fdeg[pese], row[:-1]
									HEADER[pese] = fname
									HEADERUNI[pese] += row#[:-1]
									HEADERDEG[pese] += row#[:-1]
							else:
								# print>>funi[pese], row[:-1]
								# print>>fdeg[pese], row[:-1]
								HEADER[pese] = fname
								HEADERUNI[pese] += row#[:-1]
								HEADERDEG[pese] += row#[:-1]
							continue
							
						elif row[0]=='@': continue
						# done fetching header
						
						# get first read
						try:
							lst = row[:-1].split('\t')
							hitid = lst[0]; chrom = lst[2]
							FLAG = int(lst[1])
							rem = ' '.join(lst[11:len(lst)])
							hasmate = lst[6]; insertsize=int(lst[8])
						except Exception, e:
							pipeit('Exception %s: %s'%(e, row[:-1]))
							continue
						
						break
					fh.close()
					
					# now as we find reads, group by hitid to filter duplicates
					fh = open(fname)
					count = 0
					
					# iterate through the merged+sorted map file
					try:
						# group all entries with identical hitid and mapping location
						hdbuffer = []
						currtrip = None
						focaltrip = None
						focallst = None
						focalrow = None
						
						while True:
							if count % 5000000 == 0: pipeit('%.1e,'%(count))
							count += 1
							
							# print 'newline | focal', focaltrip, focallst
							
							while focaltrip == currtrip:
								row = fh.next()
								if row[0]=='@': continue
								
								try:
									lst = row[:-1].split('\t')
									hitid = lst[0]; chrom = lst[2]
									FLAG = int(lst[1])
									rem = ' '.join(lst[11:len(lst)])
									hasmate = lst[6]; insertsize=int(lst[8])
								except Exception, e:
									pipeit('Exception %s: %s'%(e, row[:-1]))
									continue
								
								# print 'looking:', lst
								
								# toss unmapped reads
								if chrom == '*': continue
								
								# if no valid alignment
								if FLAG&(0x4): pipeit(row,0,finalbads); continue
								
								# if marked as mapped but CIGAR string unavailable, then toss
								if lst[5] == '*': pipeit(row,0,finalbads); continue
								
								# if marked as paired but the mate is mapped to same location, then toss
								# if FLAG&(0x1) and lst[3] == lst[7]: pipeit(row,0,finalbads); continue
								# print 'pass'
								
								# we have a valid read alignment
								focaltrip = (hitid,chrom,lst[3])
								focallst = lst
								focalrow = row
								if not currtrip: currtrip = focaltrip
								
								if focaltrip == currtrip: hdbuffer += [(focaltrip,focallst,focalrow)]
								# else: 
								# 	print 'mismatch'
								# 	print 'focal', focaltrip
								# 	print 'curr', currtrip
									
							# if we broke out of the above loop, we have filled the buffer
							# remove duplicates and process this batch of reads
							# ----------------------------------
							# print 'purging buffer', len(hdbuffer)
							
							# work with last row added
							ftrip,flst,frow = hdbuffer.pop()
							FLAG = int(flst[1])
							
							# resolving diploid dual-mapping issue
							# since we map to 2 reference genomes, we can have 2 records to the same location
							# order of interest
							# 1. paired
							# 2. unpaired with fewest mismatches
							# since files are sorted, keep a buffer of recent alignments
							# print 'saving part', frow
							
							keeprows = []
							# check if current read is paired and has + template length
							if FLAG&(0x1):
								if int(flst[8]) != 0:
									# keep this read and toss the one in the buffer
									# print 'have good; toss existing', lst
									for crap,lstb,rowb in hdbuffer: pipeit(rowb,0,finalbads)
									keeprows += [frow]
								# 0 template length - but check if there are any competing paired alignments
								elif not len(hdbuffer):
									# save this read
									keeprows += [(frow)]
								else:
									# toss this read
									pipeit(frow,0,finalbads)
									# keep the others
									for crap,lstb,rowb in hdbuffer: keeprows += [rowb]
							else: # unpaired
								# unpaired read - choose among equivalent alignments
								# print 'keeping', frow
								# for tri in hdbuffer: print 'keeping', tri[2]
								# for now, keep all
								keeprows += [frow]
								for crap,lstb,rowb in hdbuffer: keeprows += [rowb]
								
							# set the new current triple
							currtrip = focaltrip
							# set new buffer
							hdbuffer = [(focaltrip,focallst,focalrow)]
							# ----------------------------------
							
							
							# finally, organize and store this batch of curated reads
							# -------------------------------------------------------
							for row in keeprows:
								lst = row[:-1].split('\t')
								hitid = lst[0]; chrom = lst[2]
								FLAG = int(lst[1])
								rem = ' '.join(lst[11:len(lst)])
								hasmate = lst[6]; insertsize=int(lst[8])
								
								preinfo = ''.join([rlp,'\t'.join(lst[:11])])
								prelinen = None
								if hasmate != '*' and (hitid[-2:] == '/1' or hitid[-2:] == '/2'):
									hitid = hitid[:-2]
									lst[0] = hitid
									prelinen = ''.join([rlp,'\t'.join(lst)])
								else:
									prelinen = ''.join([rlp,row[:-1]]) # prepend run,lane information
									
									
								# don't double count, so only record positive values
								if count < 3: count = 3
								# only keep a fraction of them to save on memory
								if hasmate != '*' and insertsize >= 0 and random.random() < 2./ln(count): isizes[hapdip] += [insertsize]
								
								nmismatches = re.split('.*NM:i:(\d+).*', rem)
								nmismatches = (len(nmismatches)>1 and int(nmismatches[1])) or 0
								
								# alignment score?
								alnscore = re.split('.*AS:i:(\d+).*', rem)
								alnscore = (len(alnscore)>1 and int(alnscore[1])) or 0
								
								# sam format now eliminates /1 and /2 for PE reads that map as PE
								# thus we need to use flag info for this
								# However, older bowtie versions keep the /1 and /2, so trim
								
								if keyid == None:
									keyid = hitid
									keypese = hasmate == '*' and 'SE' or 'PE'
									keypeseint = keypese == 'PE' and 1 or 0
									ALN[preinfo] = (prelinen,alnscore)
									# MMS = [nmismatches]
									# UNIQUEPARTALN = [preinfo]
								# distinguish between single and paired
								elif keyid == hitid:
									# ALN += [prelinen]
									# MMS += [nmismatches]
									# UNIQUEPARTALN += [preinfo]
									# print 'testing second hit', preinfo, nmismatches
									if preinfo in ALN:
										# print 'testing second hit', preinfo, nmismatches
										if alnscore > ALN[preinfo][1]: ALN[preinfo] = (prelinen,alnscore)
									else:
										ALN[preinfo] = (prelinen,alnscore)
								
								elif keyid != hitid:
									# we encountered the next read
									KEEPREAD = True
							
									# bowtie2 may have a bug where unaligned reads are printed to every file
									# make sure only unique lines are printed
									ALN,MMS = unzip(ALN.values())
							
									# save previous read/alignments to the appropriate file
									# -----------------------------------------------------
									# if PE then write to PE.unique if there are only 2 entries
									if keypese == 'PE' and len(ALN) == 2:
										print >>funi[hapdip][pese], ALN[0]
										print >>funi[hapdip][pese], ALN[1]
										# get our own line count for proper sam file format
										try: ulinesat[hapdip][keypese][chrom] += 2
										except KeyError: ulinesat[hapdip][keypese][chrom] = 2
									elif keypese == 'SE' and len(ALN) == 1:
										print >>funi[hapdip][pese], ALN[0]
										# get our own line count for proper sam file format
										try: ulinesat[hapdip][keypese][chrom] += 1
										except KeyError: ulinesat[hapdip][keypese][chrom] = 1
									elif keypese == 'PE' and pese == 'SE' and len(ALN) == 1:
										print >>funi[hapdip][pese], ALN[0]
										# get our own line count for proper sam file format
										try: ulinesat[hapdip][keypese][chrom] += 1
										except KeyError: ulinesat[hapdip][keypese][chrom] = 1
									else:
								
										# determine if one alignment has fewest mismatches
										# if not, then toss all of them
										# this works for SE reads
										if STRATEGY=='bestno':
											# what about PE reads? 
											# need to group mismatch pairs by reads
											# this assumes mate pairs are interleaved
											omms = MMS
											if keypese == 'PE':
												MMS = map(lambda x,y: (x,y), MMS[:len(MMS)/2],MMS[len(MMS)/2:])
										
											if MMS.count(min(MMS)) > 1: KEEPREAD = False
									
											if KEEPREAD:
												# store the single best read (by mismatches)
												print >>funi[hapdip][pese], ALN[argmin(MMS)]
												# get our own line count for proper sam file format
												try: ulinesat[hapdip][keypese][chrom] += 1
												except KeyError: ulinesat[hapdip][keypese][chrom] = 1
											else:
												for li in xrange(len(ALN)):
													print >>fbestno, '\t'.join([ALN[li],str(omms[li])])
										else:
											# always keep read - store all alignments here
											for linen in ALN:
												print >>fdeg[hapdip][pese], linen
												# get our own line count for proper sam file format
												try: dlinesat[hapdip][keypese][chrom] += 1
												except KeyError: dlinesat[hapdip][keypese][chrom] = 1
							
									if KEEPREAD:
										# record histogram information
										if keypese == 'PE' and pese == 'SE' and len(ALN) == 1:
											try: RIDm[(keypeseint,len(ALN))] += 1
											except KeyError: RIDm[(keypeseint,len(ALN))] = 1
										elif keypese == 'PE': 
											try: RIDm[(keypeseint,len(ALN)/2)] += 2
											except KeyError: RIDm[(keypeseint,len(ALN)/2)] = 2
										else:
											try: RIDm[(keypeseint,len(ALN))] += 1
											except KeyError: RIDm[(keypeseint,len(ALN))] = 1
								
										# record mismatches
										for nmm in MMS:
											# npmms += (len(mmtmp)>1 and int(mmtmp[1])>0 and 1) or 0
											npmms[hapdip] += nmm > 0 and 1 or 0
											nqmms[hapdip] += nmm
											try: hqmms[hapdip][nmm] += 1
											except KeyError: hqmms[hapdip][nmm] = 1
											nrows[hapdip] += 1
								
									# -----------------------------------------------------
							
									# create new target
									keyid = hitid
									keypese = hasmate == '*' and 'SE' or 'PE'
									keypeseint = keypese == 'PE' and 1 or 0
									# ALN = [prelinen]
									# MMS = [nmismatches]
									ALN = {preinfo:(prelinen,nmismatches)}
					except StopIteration: pass
					
					# finally, spit out the last read(s) in the file (copy of above code)
					# -------------------------------------------------------------------
					
					# work with last row added
					try:
						ftrip,flst,frow = hdbuffer.pop()
					except IndexError:
						pipeit('Error: empty read map file?',1)
						continue
					
					FLAG = int(flst[1])
					
					keeprows = []
					# check if current read is paired and has + template length
					if FLAG&(0x1):
						if int(flst[8]) != 0:
							# keep this read and toss the one in the buffer
							# print 'have good; toss existing', lst
							for crap,lstb,rowb in hdbuffer: pipeit(rowb,0,finalbads)
							keeprows += [frow]
						# 0 template length - but check if there are any competing paired alignments
						elif not len(hdbuffer):
							# save this read
							keeprows += [(frow)]
						else:
							# toss this read
							pipeit(frow,0,finalbads)
							# keep the others
							for crap,lstb,rowb in hdbuffer: keeprows += [rowb]
					else: # unpaired
						# unpaired read - choose among equivalent alignments
						# for now, keep all
						keeprows += [frow]
						for crap,lstb,rowb in hdbuffer: keeprows += [rowb]
					
					# finally, organize and store this final batch of curated reads
					# -------------------------------------------------------
					for row in keeprows:
						lst = row[:-1].split('\t')
						hitid = lst[0]; chrom = lst[2]
						FLAG = int(lst[1])
						rem = ' '.join(lst[11:len(lst)])
						hasmate = lst[6]; insertsize=int(lst[8])
								
						preinfo = ''.join([rlp,'\t'.join(lst[:11])])
						prelinen = None
						if hasmate != '*' and (hitid[-2:] == '/1' or hitid[-2:] == '/2'):
							hitid = hitid[:-2]
							lst[0] = hitid
							prelinen = ''.join([rlp,'\t'.join(lst)])
						else:
							prelinen = ''.join([rlp,row[:-1]]) # prepend run,lane information
									
									
						# don't double count, so only record positive values
						if count < 3: count = 3
						# only keep a fraction of them to save on memory
						if hasmate != '*' and insertsize >= 0 and random.random() < 2./ln(count): isizes[hapdip] += [insertsize]
								
						nmismatches = re.split('.*NM:i:(\d+).*', rem)
						nmismatches = (len(nmismatches)>1 and int(nmismatches[1])) or 0
								
						# alignment score?
						alnscore = re.split('.*AS:i:(\d+).*', rem)
						alnscore = (len(alnscore)>1 and int(alnscore[1])) or 0
								
						# sam format now eliminates /1 and /2 for PE reads that map as PE
						# thus we need to use flag info for this
						# However, older bowtie versions keep the /1 and /2, so trim
								
						if keyid == None:
							keyid = hitid
							keypese = hasmate == '*' and 'SE' or 'PE'
							keypeseint = keypese == 'PE' and 1 or 0
							ALN[preinfo] = (prelinen,alnscore)
							# MMS = [nmismatches]
							# UNIQUEPARTALN = [preinfo]
						# distinguish between single and paired
						elif keyid == hitid:
							# ALN += [prelinen]
							# MMS += [nmismatches]
							# UNIQUEPARTALN += [preinfo]
							# print 'testing second hit', preinfo, nmismatches
							if preinfo in ALN:
								# print 'testing second hit', preinfo, nmismatches
								if alnscore > ALN[preinfo][1]: ALN[preinfo] = (prelinen,alnscore)
							else:
								ALN[preinfo] = (prelinen,alnscore)
								
						elif keyid != hitid:
							# we encountered the next read
							KEEPREAD = True
							
							# bowtie2 may have a bug where unaligned reads are printed to every file
							# make sure only unique lines are printed
							ALN,MMS = unzip(ALN.values())
							
							# save previous read/alignments to the appropriate file
							# -----------------------------------------------------
							# if PE then write to PE.unique if there are only 2 entries
							if keypese == 'PE' and len(ALN) == 2:
								print >>funi[hapdip][pese], ALN[0]
								print >>funi[hapdip][pese], ALN[1]
								# get our own line count for proper sam file format
								try: ulinesat[hapdip][keypese][chrom] += 2
								except KeyError: ulinesat[hapdip][keypese][chrom] = 2
							elif keypese == 'SE' and len(ALN) == 1:
								print >>funi[hapdip][pese], ALN[0]
								# get our own line count for proper sam file format
								try: ulinesat[hapdip][keypese][chrom] += 1
								except KeyError: ulinesat[hapdip][keypese][chrom] = 1
							elif keypese == 'PE' and pese == 'SE' and len(ALN) == 1:
								print >>funi[hapdip][pese], ALN[0]
								# get our own line count for proper sam file format
								try: ulinesat[hapdip][keypese][chrom] += 1
								except KeyError: ulinesat[hapdip][keypese][chrom] = 1
							else:
								
								# determine if one alignment has fewest mismatches
								# if not, then toss all of them
								# this works for SE reads
								if STRATEGY=='bestno':
									# what about PE reads? 
									# need to group mismatch pairs by reads
									# this assumes mate pairs are interleaved
									omms = MMS
									if keypese == 'PE':
										MMS = map(lambda x,y: (x,y), MMS[:len(MMS)/2],MMS[len(MMS)/2:])
										
									if MMS.count(min(MMS)) > 1: KEEPREAD = False
									
									if KEEPREAD:
										# store the single best read (by mismatches)
										print >>funi[hapdip][pese], ALN[argmin(MMS)]
										# get our own line count for proper sam file format
										try: ulinesat[hapdip][keypese][chrom] += 1
										except KeyError: ulinesat[hapdip][keypese][chrom] = 1
									else:
										for li in xrange(len(ALN)):
											print >>fbestno, '\t'.join([ALN[li],str(omms[li])])
								else:
									# always keep read - store all alignments here
									for linen in ALN:
										print >>fdeg[hapdip][pese], linen
										# get our own line count for proper sam file format
										try: dlinesat[hapdip][keypese][chrom] += 1
										except KeyError: dlinesat[hapdip][keypese][chrom] = 1
							
							if KEEPREAD:
								# record histogram information
								if keypese == 'PE' and pese == 'SE' and len(ALN) == 1:
									try: RIDm[(keypeseint,len(ALN))] += 1
									except KeyError: RIDm[(keypeseint,len(ALN))] = 1
								elif keypese == 'PE': 
									try: RIDm[(keypeseint,len(ALN)/2)] += 2
									except KeyError: RIDm[(keypeseint,len(ALN)/2)] = 2
								else:
									try: RIDm[(keypeseint,len(ALN))] += 1
									except KeyError: RIDm[(keypeseint,len(ALN))] = 1
								
								# record mismatches
								for nmm in MMS:
									# npmms += (len(mmtmp)>1 and int(mmtmp[1])>0 and 1) or 0
									npmms[hapdip] += nmm > 0 and 1 or 0
									nqmms[hapdip] += nmm
									try: hqmms[hapdip][nmm] += 1
									except KeyError: hqmms[hapdip][nmm] = 1
									nrows[hapdip] += 1
								
							# -----------------------------------------------------
							
							# create new target
							keyid = hitid
							keypese = hasmate == '*' and 'SE' or 'PE'
							keypeseint = keypese == 'PE' and 1 or 0
							# ALN = [prelinen]
							# MMS = [nmismatches]
							ALN = {preinfo:(prelinen,nmismatches)}
					# ================================================================
					
					# catch exception of 0 mapped reads
					try: ALN,MMS = unzip(ALN.values())
					except Exception,e: continue
					
					# we encountered the next read
					KEEPREAD = True
					
					# save previous read/alignments to the appropriate file
					# -----------------------------------------------------
					# if PE then write to PE.unique if there are only 2 entries
					if keypese == 'PE' and len(ALN) == 2:
						print >>funi[hapdip][pese], ALN[0]
						print >>funi[hapdip][pese], ALN[1]
						# get our own line count for proper sam file format
						try: ulinesat[hapdip][keypese][chrom] += 2
						except KeyError: ulinesat[hapdip][keypese][chrom] = 2
					elif keypese == 'SE' and len(ALN) == 1:
						print >>funi[hapdip][pese], ALN[0]
						# get our own line count for proper sam file format
						try: ulinesat[hapdip][keypese][chrom] += 1
						except KeyError: ulinesat[hapdip][keypese][chrom] = 1
					elif keypese == 'PE' and pese == 'SE' and len(ALN) == 1:
						print >>funi[hapdip][pese], ALN[0]
						# get our own line count for proper sam file format
						try: ulinesat[hapdip][keypese][chrom] += 1
						except KeyError: ulinesat[hapdip][keypese][chrom] = 1
					else:
						# determine if one alignment has fewest mismatches
						# if not, then toss all of them
						# this works for SE reads
						if STRATEGY=='bestno':
							# what about PE reads? 
							# need to group mismatch pairs by reads
							# this assumes mate pairs are interleaved
							omms = MMS
							if keypese == 'PE':
								MMS = map(lambda x,y: (x,y), MMS[:len(MMS)/2],MMS[len(MMS)/2:])
								
								if MMS.count(min(MMS)) > 1: KEEPREAD = False
							
								if KEEPREAD:
									# store the single best read (by mismatches)
									print >>funi[hapdip][pese], ALN[argmin(MMS)]
									# get our own line count for proper sam file format
									try: ulinesat[hapdip][keypese][chrom] += 1
									except KeyError: ulinesat[hapdip][keypese][chrom] = 1
								else:
									for li in xrange(len(ALN)):
										print >>fbestno, '\t'.join([ALN[li],str(omms[li])])
								
						else:
							# always keep read - store all alignments here
							for linen in ALN:
								print >>fdeg[hapdip][pese], linen
								# get our own line count for proper sam file format
								try: dlinesat[hapdip][keypese][chrom] += 1
								except KeyError: dlinesat[hapdip][keypese][chrom] = 1
					
					if KEEPREAD:
						# record histogram information
						if keypese == 'PE' and pese == 'SE' and len(ALN) == 1:
							try: RIDm[(keypeseint,len(ALN))] += 1
							except KeyError: RIDm[(keypeseint,len(ALN))] = 1
						elif keypese == 'PE': 
							try: RIDm[(keypeseint,len(ALN)/2)] += 2
							except KeyError: RIDm[(keypeseint,len(ALN)/2)] = 2
						else:
							try: RIDm[(keypeseint,len(ALN))] += 1
							except KeyError: RIDm[(keypeseint,len(ALN))] = 1
						
						# record mismatches
						for nmm in MMS:
							# npmms += (len(mmtmp)>1 and int(mmtmp[1])>0 and 1) or 0
							npmms[hapdip] += nmm > 0 and 1 or 0
							nqmms[hapdip] += nmm
							try: hqmms[hapdip][nmm] += 1
							except KeyError: hqmms[hapdip][nmm] = 1
							nrows[hapdip] += 1
					# -------------------------------------------------------------------
					
					fh.close()
					sys.stdout.write('done (n=%s).\n'%(count)); sys.stdout.flush()
				else: pipeit('Cannot access %s'%(fname),1)
				# end loop reading read info
				
				# Aggregate the multiplicity histogram information
				for (keypeseint,val),vct in RIDm.iteritems():
					keypese = keypeseint == 0 and 'SE' or 'PE'
					try: sidHist[hapdip][keypese][val] += vct
					except KeyError: sidHist[hapdip][keypese][val] = vct
					
					try: sidHist[hapdip]['Total'][val] += vct
					except KeyError: sidHist[hapdip]['Total'][val] = vct
					
					totreads[hapdip][keypese] += vct
					totaln[hapdip][keypese] += val*vct
				del RIDm
				
			# finished processing component read map file for this map file
			if STRATEGY=='bestno': fbestno.close()
			finalbads.close()
		
		
		hds = [hapdip]
		if ALLELESPECIFIC: hds = ['h1','h2']
		
		for hapdip in hds:
			for pese in ['PE','SE']:
				funi[hapdip][pese].close(); fdeg[hapdip][pese].close()
		
		# Sort the unique maps
		if ALLELESPECIFIC:
			# if user wants AS, then this must be bowtie2 so there is no SE/PE division option
			for hd in ['h1','h2']:
				stmp = btmapdir+sampleID+'.%s.unique.map'%(hd)
				funi = {'SE':stmp, 'PE':stmp}
				stmp = btmapdir+sampleID+'.%s.degenerate.map'%(hd)
				fdeg = {'SE':stmp, 'PE':stmp}
				
				pipeit('- Sorting %s...'%(getFilename(funi['SE'])))
				os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s "%s" > "%s"'%(nsubprocs, funi['SE'], funi['SE']+'.sorted'))
				os.system('mv "%s" "%s"'%(funi['SE']+'.sorted', funi['SE']))
				
				# DONE SORTING
				# Add the header into files
				pipeit('Adding SAM header...')
				
				fho = open(funi['SE'][:-4]+'.header.txt', 'w')
				print >>fho, HEADERUNI['SE'],
				fho.close()
				# concatenate data onto header
				os.system('cat "%s" >> "%s"'%(funi['SE'], funi['SE'][:-4]+'.header.txt'))
				# rename
				os.system('mv "%s" "%s"'%(funi['SE'][:-4]+'.header.txt', funi['SE']))
				
				fho = open(fdeg['SE'][:-4]+'.header.txt', 'w')
				print >>fho, HEADERDEG['SE'],
				fho.close()
				# concatenate data onto header
				os.system('cat "%s" >> "%s"'%(fdeg['SE'], fdeg['SE'][:-4]+'.header.txt'))
				# rename
				os.system('mv "%s" "%s"'%(fdeg['SE'][:-4]+'.header.txt', fdeg['SE']))
				pipeit('done.',1)
		else:
			if UDPESECOMBO:
				stmp = btmapdir+sampleID+'.map'
				funi = {'SE':stmp, 'PE':stmp}
				fdeg = {'SE':stmp, 'PE':stmp}
				
				pipeit('- Sorting %s...'%(getFilename(funi['SE'])))
				os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, funi['SE'], funi['SE']+'.sorted'))
				os.system('mv "%s" "%s"'%(funi['SE']+'.sorted', funi['SE']))
				pipeit('done.',1)
			
			elif PESECOMBO:
				stmp = btmapdir+sampleID+'.unique.map'
				funi = {'SE':stmp, 'PE':stmp}
				stmp = btmapdir+sampleID+'.degenerate.map'
				fdeg = {'SE':stmp, 'PE':stmp}
			
				pipeit('- Sorting %s...'%(getFilename(funi['SE'])))
				os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, funi['SE'], funi['SE']+'.sorted'))
				os.system('mv "%s" "%s"'%(funi['SE']+'.sorted', funi['SE']))
				pipeit('done.',1)
			
			else:
				# Partition reads into 2 multiplicity groups
				funi = {'SE':btmapdir+sampleID+'.SE.unique.map', 'PE':btmapdir+sampleID+'.PE.unique.map'}
				fdeg = {'SE':btmapdir+sampleID+'.SE.degenerate.map', 'PE':btmapdir+sampleID+'.PE.degenerate.map'}
			
				for pese in ['SE','PE']:
					pipeit('- Sorting %s...'%(getFilename(funi[pese])))
					os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s -T "%s ""%s" > "%s"'%(nsubprocs, outdir, funi[pese], funi[pese]+'.sorted'))
					os.system('mv "%s" "%s"'%(funi[pese]+'.sorted', funi[pese]))
					pipeit('done.',1)
			# DONE SORTING
			
			# Add the header into files
			pipeit('- Adding SAM header...')
			if UDPESECOMBO:
				fho = open(funi['SE'][:-4]+'.header.txt', 'w')
				print >>fho, HEADERUNI['SE'],
				fho.close()
				# concatenate data onto header
				os.system('cat "%s" >> "%s"'%(funi['SE'], funi['SE'][:-4]+'.header.txt'))
				# rename
				os.system('mv "%s" "%s"'%(funi['SE'][:-4]+'.header.txt', funi['SE']))
			
			elif PESECOMBO:
				fho = open(funi['SE'][:-4]+'.header.txt', 'w')
				print >>fho, HEADERUNI['SE'],
				fho.close()
				# concatenate data onto header
				os.system('cat "%s" >> "%s"'%(funi['SE'], funi['SE'][:-4]+'.header.txt'))
				# rename
				os.system('mv "%s" "%s"'%(funi['SE'][:-4]+'.header.txt', funi['SE']))
			
				fho = open(fdeg['SE'][:-4]+'.header.txt', 'w')
				print >>fho, HEADERDEG['SE'],
				fho.close()
				# concatenate data onto header
				os.system('cat "%s" >> "%s"'%(fdeg['SE'], fdeg['SE'][:-4]+'.header.txt'))
				# rename
				os.system('mv "%s" "%s"'%(fdeg['SE'][:-4]+'.header.txt', fdeg['SE']))
			
			else:
				for pese in ['SE', 'PE']:
					fho = open(funi[pese][:-4]+'.header.txt', 'w')
					print >>fho, HEADERUNI[pese],
					fho.close()
					# concatenate data onto header
					os.system('cat "%s" >> "%s"'%(funi[pese], funi[pese][:-4]+'.header.txt'))
					# rename
					os.system('mv "%s" "%s"'%(funi[pese][:-4]+'.header.txt', funi[pese]))
				
					fho = open(fdeg[pese][:-4]+'.header.txt', 'w')
					print >>fho, HEADERDEG[pese],
					fho.close()
					# concatenate data onto header
					os.system('cat "%s" >> "%s"'%(fdeg[pese], fdeg[pese][:-4]+'.header.txt'))
					# rename
					os.system('mv "%s" "%s"'%(fdeg[pese][:-4]+'.header.txt', fdeg[pese]))
			pipeit('done.',1)
		
		# ----------------------------------------
		# Finished all map files for this sampleID
		
		# Save results for this sampleID / haploid/diploid
		recfile = almfile
		record(sampleID, recfile, 'a')
		# report total reads for this sample
		lst = []
		totalcount = 0
		if sampleID not in rawfiles:
			pipeit('Warning: will not obtain total raw reads for sample %s.'%(sampleID),1)
		# elif len(filter(lambda x: getSuffix(x) in ['bz2', 'gz','zip'], rawfiles[sampleID])):
		else:
			for f in rawfiles[sampleID]:
				if getSuffix(f) in ['bz2', 'gz','zip']:
					os.system('zcat "%s" | wc -l > "%s"'%(f,zdir+'tmpfile.txt'))
				else:
					os.system('wc -l "%s" > "%s"'%(f,zdir+'tmpfile.txt'))
				
				count = 0
				try:
					info = readList(zdir+'tmpfile.txt')[0].strip().split(' ')
					count = int(int(info.pop(0))/4.) # 4 line fasta file
					lst += ['%s\t%s'%(count,' '.join(info))]
				except IndexError:
					# if we can no longer access the fastq files
					lst += ['%s\t%s'%(count,' '.join(['NA']))]
				totalcount += count
		record('\n'.join(lst),recfile,'a')
		record('%s\tTotal'%(totalcount),recfile,'a')
		os.system('rm -f "%s"'%(zdir+'tmpfile.txt'))
		
		hdkeepm = []
		hdkeepd = []

		row = []
		for hd in hds:
			# move on to next group if there are no mapped reads
			# if nrows[hd] == 0: continue
			
			# save total read counts for this sampleID
			totreads[hd]['Total'] += totreads[hd]['SE']+totreads[hd]['PE']
			totaln[hd]['Total'] += totaln[hd]['SE']+totaln[hd]['PE']
			
			# record number of perfect/mismatch alignments
			perfrac = 0
			nper = nrows[hd] - npmms[hd]
			meanmm = divna(nqmms[hd], nrows[hd])
			
			perfrac = divna(nper,float(nrows[hd]))
			mmfrac = divna(npmms[hd],float(nrows[hd]))
			row += ['%s: There are %s/%s (%s) perfect (%s/%s) alignments'%(hd, nper,nrows[hd],perfrac,maxreadlen,maxreadlen)]
			row += ['%s: There are %s/%s (%s) mismatched (<%s/%s) alignments'%(hd, npmms[hd],nrows[hd],mmfrac,maxreadlen,maxreadlen)]
			row += ['Avg mismatches/alignment = %s'%(meanmm)]
			
			# print hqmms[hd]
			hqstr = ','.join(['%s:%s'%(k,roundna(4)(divna(hqmms[hd][k],nrows[hd]))) for k in sorted(hqmms[hd].keys())])
			if not hqstr: hqstr = nastr
			
			# save read and multiplicity statistics
			sestr = ', '.join('('+','.join(map(str,p))+')' for p in sorted(sidHist[hd]['SE'].items()))
			pestr = ', '.join('('+','.join(map(str,p))+')' for p in sorted(sidHist[hd]['PE'].items()))
			totstr = ', '.join('('+','.join(map(str,p))+')' for p in sorted(sidHist[hd]['Total'].items()))
			if sestr == '': sestr = nastr
			if pestr == '': pestr = nastr
			
			if hd != hapdip: record(hd,recfile,'a')
			record('SE: reads=%s, alignments=%s, alignments/read=%s'%(totreads[hd]['SE'], totaln[hd]['SE'], roundna(3)(divna(totaln[hd]['SE'],totreads[hd]['SE']))), recfile, 'a')
			record('Histogram: '+sestr, recfile, 'a')
			record('PE: reads=%s, alignments=%s, alignments/read=%s'%(totreads[hd]['PE'], totaln[hd]['PE'], roundna(3)(divna(totaln[hd]['PE'], totreads[hd]['PE']))), recfile, 'a')
			record('Histogram: '+pestr, recfile, 'a')
			record('Total: reads=%s, alignments=%s, alignments/read=%s'%(totreads[hd]['Total'], totaln[hd]['Total'], roundna(3)(divna(totaln[hd]['Total'], totreads[hd]['Total']))), recfile, 'a')
			frac = nastr
			if totalcount > 0: frac = round(totreads[hd]['Total']/float(totalcount), 3)
			record('Fraction mapped: %s'%(frac), recfile, 'a')
			record('Histogram: '+totstr, recfile, 'a')
			
			# reHist to summarize multiplicity distribution
			# ---------------------------------------------------
			newbins = [1,2,3,5,10,50]
			# calculate avg alignments per read for non-unique reads only
			# do this for 2,3,4,5
			hist = map(lambda x: map(intna,x.strip('(').strip(')').split(',')), totstr.split(', '))
			ttotreads = sum([foo[1] for foo in hist[1:maxnumalign-1]])
			ttotaln = sum([foo[0]*foo[1] for foo in hist[1:maxnumalign-1]])
			record('Multiplicity: '+'\t'.join(map(str,newbins)), recfile, 'a')
			# raise exception if no reads in map
			try: new = reHist(hist,newbins, freq=1)
			except ValueError: new = zip(newbins,[0 for z in newbins])
			newstr = '\t'.join(map(lambda x: str(round(x[1],5)), new))
			record('     Summary: '+newstr, recfile, 'a')
			sapr = divna(ttotaln,ttotreads)
			record('Deg aln/read: %s'%(sapr), recfile, 'a')
			
			# add summary to summary QC file
			tmpstr = ';'.join(map(lambda x: '%s:%s'%(x[0],x[1]), zip(newbins , map(lambda x: str(round(x[1],5)), new))))
			hdkeepm += [tmpstr]
			hdkeepd += [str(sapr)]

		# save perfect/mismatch reads summary
		info = '\t'.join(map(str, [sampleID, ';'.join(map(getFile, rawfiles[sampleID])), hd, totalcount, totreads[hd]['Total'], roundna(4)(divna(totreads[hd]['Total'],totalcount)), totaln[hd]['Total'], nper, npmms[hd], nrows[hd], roundna(4)(perfrac), roundna(4)(mmfrac), roundna(4)(meanmm), hqstr, ' '.join(hdkeepm), ' '.join(hdkeepd)]))
		record(info, pmafile, 'a')
			

			# ---------------------------------------------------
		record('', recfile, 'a')
		del sidHist # clean up
		printList(row, file=pmadir+'Perfect vs mismatch alignments-%s.txt'%(sampleID))
		
		# Done saving statistics for this sampleID
		pipeit('- Done partitioning unique and degenerate reads.',1)
		
		
		# ==============================================================
		# save the number of SE/PE reads on each chromosome
		fh = open(crdfile, 'a')
		# this prints the number of reads on each chromosome
		# for k in sorted(GRlengths.keys()):
		
		allchr = {}
		for hd in hds:
			for k in sorted(unique(ulinesat[hd]['SE'].keys()+ulinesat[hd]['PE'].keys())): allchr[k] = None
			
		for hd in hds:
			for k in allchr:
				if k == '*': continue # unaligned
				try:
					use = upe = dse = dpe = 0
					try: use = ulinesat[hd]['SE'][k]
					except KeyError: pass
					try: upe = ulinesat[hd]['PE'][k]
					except KeyError: pass
					try: dse = dlinesat[hd]['SE'][k]
					except KeyError: pass
					try: dpe = dlinesat[hd]['PE'][k]
					except KeyError: pass
					print >>fh, '\t'.join(map(str, [sampleID, hd, k, use, upe, dse, dpe]))
				except KeyError:
					print >>fh, '\t'.join(map(str, [sampleID, hd, k, 0, 0, 0, 0]))
				
		fh.close()
		# ==============================================================
		
		
		# =========================================
		# sort read maps (not currently implemented)
		# for pese in ['PE', 'SE']:
		# 	os.system('sort -k 2,3n "%s" > "%s"'%(btmapdir+sampleID+'.'+pese+'.unique.map', btmapdir+sampleID+'.'+pese+'.unique.map'))
		# 	os.system('sort -k 2,3n "%s" > "%s"'%(btmapdir+sampleID+'.'+pese+'.degenerate.map', btmapdir+sampleID+'.'+pese+'.degenerate.map'))
		# alternatively could use samtools for sorting...
		# =========================================
		
		
		pipeit('- Saving empirical insert size distribution...')
		for hd in hds:
			hist = histogram(isizes[hd][0:250000], bins=makeRange(minins,maxins+2,2))
			printTable(hist, header=['Size', 'Count'], file=zinsdir+'%s_%s_insert_size_distribution.txt'%(sampleID,hd))
			pipeit(hd+'...')
		pipeit('done.',1)
		
		# clean up
		del ulinesat; del dlinesat; del isizes
		
		# finally delete the precursor bowtie files
		# -----------------------------------------
		if not KEEPBOWTIE:
			for run,lane in seqRuns[sampleID].keys():
				UID = '%s-%s-%s'%(sampleID,run,lane)
				refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
				thefiles = [btmapdir+UID+'.sam', btmapdir+UID+'.h1.sorted.sam', btmapdir+UID+'.h2.sorted.sam', btmapdir+UID+'.sorted.sam', btmapdir+UID+'_singles_1.sam', btmapdir+UID+'_singles_2.sam']
				# thefiles += [btmapdir+UID+'.sam.sort', btmapdir+UID+'_singles_1.sam.sort', btmapdir+UID+'_singles_2.sam.sort']
				for fn in thefiles: os.system('rm -f "%s"'%(fn))
	
	if CLEARBOWTIE:
		pipeit('- Clearing precusor read map files...',1)
		for sampleID in sorted(seqRuns.keys()):
			for run,lane in seqRuns[sampleID].keys():
				UID = '%s-%s-%s'%(sampleID,run,lane)
				refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
				thefiles = [btmapdir+UID+'.sam', btmapdir+UID+'_singles_1.sam', btmapdir+UID+'_singles_2.sam']
				for fn in thefiles: os.system('rm -f "%s"'%(fn))

	pipeit('DONE READ ORGANIZATION.\n',1)


if FILTERPCRDUPLICATES:
	print '=> FILTERING PCR REPLICATE READS'
	filterdir = badReadDir#zdir+'libraryErrors/'
	createdir(filterdir)
	
	# idea is to filter based on sequence criteria
	# Partition into unique/deg to conserve memory
	# while making an assumption that sequencing error is negligible here
	# default sequence length is set to 36
	
	for sampleID in sortseqruns:
		if not replotFD:
			sys.stdout.write('- Filtering PCR replicate reads: %s...\n'%(sampleID)); sys.stdout.flush()
			if ALLELESPECIFIC:
				filterReplicateReads(btmapdir, filterdir, sampleID, ud='h1.unique', mm=pcrfiltermm, verbose=False, redo=redoFD, seed=pcrfilterlength, T=pcrfilterT)
				filterReplicateReads(btmapdir, filterdir, sampleID, ud='h1.degenerate', mm=pcrfiltermm, verbose=False, redo=redoFD, seed=pcrfilterlength, T=pcrfilterT)
				filterReplicateReads(btmapdir, filterdir, sampleID, ud='h2.unique', mm=pcrfiltermm, verbose=False, redo=redoFD, seed=pcrfilterlength, T=pcrfilterT)
				filterReplicateReads(btmapdir, filterdir, sampleID, ud='h2.degenerate', mm=pcrfiltermm, verbose=False, redo=redoFD, seed=pcrfilterlength, T=pcrfilterT)
			else:
				filterReplicateReads(btmapdir, filterdir, sampleID, ud='unique', mm=pcrfiltermm, verbose=False, redo=redoFD, seed=pcrfilterlength, T=pcrfilterT)
				filterReplicateReads(btmapdir, filterdir, sampleID, ud='degenerate', mm=pcrfiltermm, verbose=False, redo=redoFD, seed=pcrfilterlength, T=pcrfilterT)
		
		if replotFD or andplotFD:
			if ALLELESPECIFIC: uds = ['h1.unique','h1.degenerate','h2.unique','h2.degenerate']
			sys.stdout.write('- Plotting PCR replicate reads: %s...\n'%(sampleID)); sys.stdout.flush()
			# filter to identify clonally unique reads
			filterReplicateReads(btmapdir, filterdir, sampleID, ud=['unique','degenerate'], mm=pcrfiltermm, verbose=False, redo=redoFD,  seed=pcrfilterlength, T=pcrfilterT, replot=True)
		
	print '=> DONE FILTERING\n'

	
if ALIGNTOTRANSCRIPTOME:
	# merge the unique and degenerage files for cufflinks analysis
	# 1. duplicate the unique file
	pipeit('- Merging unique and degenerate files for cufflinks analysis...',1)
	for sampleID in sorted(seqRuns.keys()):
		# for run,lane in seqRuns[sampleID].keys():
		# 	UID = '%s-%s-%s'%(sampleID,run,lane)
		# 	refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
		pipeit('  - Sample: %s...'%(sampleID))
		if ALLELESPECIFIC:
			for hd in ['h1','h2']:
				pipeit('%s,'%(hd))
				utmp = btmapdir+sampleID+'.%s.unique.map'%(hd)
				dtmp = btmapdir+sampleID+'.%s.degenerate.map'%(hd)
				ntmp = btmapdir+sampleID+'.%s.UD.unsort.sam'%(hd)
				ntmps = btmapdir+sampleID+'.%s.UD.sort.sam'%(hd)
				ntmph = btmapdir+sampleID+'.%s.UD.sam'%(hd)
				
				if not DRYRUN and (not os.access(ntmph, os.F_OK) or REORGANIZEREADS):
					# fetch the sam header file from unique
					if not USERPROVIDEDSAM: 
						samheaderfile = btmapdir+sampleID+'.headerfile.txt'
						os.system('samtools view -S -H "%s" > "%s"'%(utmp, samheaderfile))
						# make sure the header file has unique entries only
						lst = readList(samheaderfile)
						samdct = [lst[0]] # index is line order, value is unique entry
						for entry in lst[1:]:
							if entry != samdct[-1]: samdct += [entry]
						printList(samdct,file=samheaderfile)
					else:
						pipeit('(SAM header:%s)...'%(samheaderfile))
						
					if not os.access(ntmps, os.F_OK) or REORGANIZEREADS:
						if not os.access(ntmp, os.F_OK) or REORGANIZEREADS: 
							mergeSamFiles([utmp,dtmp],outfile=ntmp)
						# now must sort for cufflinks
						pipeit('sorting,')
						os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s "%s" > "%s"'%(nsubprocs, ntmp,ntmps))
						os.system('rm -f "%s"'%(ntmp))
					# cat the header
					os.system('cat "%s" "%s" > "%s"'%(samheaderfile, ntmps, ntmph))
					os.system('rm -f "%s"'%(ntmps))
					if not USERPROVIDEDSAM: 
						pipeit('removing sam file %s'%(samheaderfile),1)
						os.system('rm -f "%s"'%(samheaderfile))
						samheaderfile = ''
		else:
			utmp = btmapdir+sampleID+'.unique.map'
			dtmp = btmapdir+sampleID+'.degenerate.map'
			ntmp = btmapdir+sampleID+'.UD.unsort.sam'
			ntmps = btmapdir+sampleID+'.UD.sort.sam'
			ntmph = btmapdir+sampleID+'.UD.sam'
			
			# this is on hold until i get it to be cufflinks compatible
			if not DRYRUN and (not os.access(ntmph, os.F_OK) or REORGANIZEREADS):
				# fetch the sam header file from unique
				if not USERPROVIDEDSAM: 
					samheaderfile = btmapdir+sampleID+'.headerfile.txt'
					os.system('samtools view -S -H "%s" > "%s"'%(utmp, samheaderfile))
					# make sure the header file has unique entries only
					lst = readList(samheaderfile)
					samdct = [lst[0]] # index is line order, value is unique entry
					for entry in lst[1:]:
						if entry != samdct[-1]: samdct += [entry]
					printList(samdct,file=samheaderfile)
				else:
					pipeit('(SAM header:%s)...'%(samheaderfile))
					
				if not os.access(ntmps, os.F_OK) or REORGANIZEREADS:
					# pipeit('%s not found,'%(getLabel(ntmps)))
					if not os.access(ntmp, os.F_OK) or REORGANIZEREADS: 
						# pipeit('%s not found,'%(getLabel(ntmp)))
						mergeSamFiles([utmp,dtmp],outfile=ntmp)
					# now must sort by position for cufflinks
					pipeit('sorting,')
					os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s "%s" > "%s"'%(nsubprocs, ntmp,ntmps))
					os.system('rm -f "%s"'%(ntmp))
				# cat the header
				os.system('cat "%s" "%s" > "%s"'%(samheaderfile, ntmps, ntmph))
				os.system('rm -f "%s"'%(ntmps))
				if not USERPROVIDEDSAM: 
					pipeit('removing sam file %s'%(samheaderfile),1)
					os.system('rm -f "%s"'%(samheaderfile))
					samheaderfile = ''
				
		pipeit('done.',1)
	pipeit('- Done merging sam files.',1)
	
if RUNCUFFLINKS:
	pipeit('- Running cufflinks...',1)
	cuffdir = zdir+'cufflinks/'
	if CUFFLABEL: cuffdir = zdir+'cufflinks.%s/'%(CUFFLABEL)
	createdir(cuffdir)
	
	cuffPipeFile = zdir+'run_cufflinks.log'
	cuffPipe = open(cuffPipeFile, 'a')
	
	for sampleID in sorted(seqRuns.keys()):
		fraglen = fragDiameter4Sample[sampleID]
		
		RG = []
		LT = []

		for run,lane in seqRuns[sampleID].keys():
			refGenome,hashet,paired,path = seqRuns[sampleID][(run,lane)][:4]
			RG += [refGenome]
			# print 'items', seqRuns[sampleID][(run,lane)]
			LT += [seqRuns[sampleID][(run,lane)][-1]]
		LT = unique(LT)

		RG = unique(RG)
		if len(RG) > 1: sys.exit('Multiple reference genomes provided. Use -r to specify one file or update mapfile.\n%s'%(RG))
		refGenome = RG[0]

		if len(LT) > 1:
			sys.exit('mixing libraries of multiple chemistries: %s'%(LT))
		libtype = LT[0]
		pipeit('  - Sample: %s [frag length=%s, lib=%s]...'%(sampleID,fraglen,libtype))
		if ALLELESPECIFIC:
			for hd in ['h1','h2']:
				sys.exit("not yet implemented")
				# pipeit('%s,'%(hd))
				# ntmp = btmapdir+sampleID+'.%s.UD.sam'%(hd)
				# if not DRYRUN and (not os.access(cuffdir+'cufflinks.%s.%s'%(hd,sampleID), os.F_OK) or REDOCUFFLINKS):
				# 	# RUN HERE
				# 	sys.exit("not yet implemented")
		else:
			cuffdir2 = cuffdir+'cufflinks.%s/'%(sampleID)
			createdir(cuffdir2)
			
			if not DRYRUN and (not os.access(cuffdir2+'isoforms.fpkm_tracking', os.F_OK) or REDOCUFFLINKS):
				# find refGenome suffix
				refGenomeFull = refGenome
				if not os.access(refGenome,os.F_OK):
					candidates = getFiles(getPath(refGenome),include=getLabel(refGenome))
					final = filter(lambda x: x in [refGenome+'.fasta',refGenome+'.fa'], candidates)
					# print 'final', final
					if len(final) > 1:
						pipeit('Warning: multiple reference genomes found: %s. Using the first entry.'%(final),1)
					refGenomeFull = final[0]
				if refGenomeFull == None:
					sys.exit('ERROR --cufflinks no reference genome provided. Use -r flag.\n')
					
				options = '-o "%s" -p %s -m %s --library-type %s -b "%s" %s'%(cuffdir2,nsubprocs,fraglen,libtype,refGenomeFull,CUFFLINKSFLAGS)
				if gfffile: options += ' --GTF "%s"'%(gfffile)
				
				if maxnumalign > 1: 
					pipeit('Cufflinks multi-read correction ENabled (-d %s)'%(maxnumalign),1)
					options += ' -u'
				else:
					pipeit('Cufflinks multi-read correction DISabled (-d %s)'%(maxnumalign),1)
				call = '%scufflinks %s "%s"'%(cufflinksdir, options, btmapdir+sampleID+'.UD.sam')
				print '\nCall>', call
				pipeit(call,2,cuffPipe)
				os.system(call)
			elif not DRYRUN:
				pipeit('already processed.',1)
				
		pipeit('done.',1)
	
	pipeit('- Merging cufflinks files...')
	# merge FPKM columns into single matrix
	labels = []; pairs = []
	for sampleID in sorted(seqRuns.keys()):
		cuffdir2 = cuffdir+'cufflinks.%s/'%(sampleID)
		outfile = cuffdir2+'isoforms.fpkm_tracking'
		try:
			d,r,c = readTable(outfile,rownames=0)
			idx = c.index('FPKM')
			fpkm = [(row[0],row[idx]) for row in d if len(row) >= idx]
			pairs += [fpkm]
			labels += [sampleID]
		except IOError: pass
	
	# outer join - fill with NA
	names, mat = joinSet(pairs,inner=0)
	printTable(mat,['Gene']+labels, names, file=cuffdir+'%s.fpkm.txt'%(zkey))
	pipeit('- Done cufflinks analysis.',1)
	
	cuffPipe.close()


# link filenames if user supplies their own unique map file
if LINKUNIQUE:
	pipeit('- Linking existing unique maps...',1)
	for sampleID in sortseqruns:
		sourcefile = None
		tmp = '%s%s%s'%(btmapdir,sampleID,mapsuffix)
		if os.access(tmp, os.F_OK): sourcefile = tmp
		else: 
			pipeit('  Error: cannot find file %s'%(tmp),1)
			pipeit('  Check --readmapdir <PATH/TO/MAP.FILE>',1)
			continue
		targetfile = '%s%s.SE.unique.map'%(btmapdir,sampleID)
		if os.access(targetfile, os.F_OK): 
			pipeit('  %s already exists.'%(targetfile),1)
			break
		pipeit('- Linking map file %s to %s'%(sourcefile.split('/')[-1], targetfile.split('/')[-1]),1)
		os.system('ln -s "%s" "%s"'%(sourcefile,targetfile))




# THIS IS PROBLEMATIC FOR THE HEADER!
if INDEXMAPFILES or INDEXUNI:
	# sort unique read maps
	pipeit('- Sorting and indexing unique maps...', 1)
	
	peses = ['PE.', 'SE.']
	if PESECOMBO or UDPESECOMBO: peses = ['']
	
	ANYSAMPLES = False
	for sampleID in sorted(seqRuns.keys()):
		for pese in peses:
			bases = [btmapdir+sampleID+'.'+pese+'unique.map']
			if ALLELESPECIFIC:
				bases = [btmapdir+sampleID+'.%s.%sunique.map'%('h1',pese), btmapdir+sampleID+'.%s.%sunique.map'%('h2',pese)]
				
			for base in bases:
				# print 'base', base
				if not os.access(base, os.F_OK): continue
				if not os.access(base+'.index', os.F_OK) or RESORT:
					ANYSAMPLES = True
					pipeit('  %s...'%(getLabel(base)))
					if not INDEXONLY:
						# os.system('mv "%s" "%s"'%(base, base+'.unsort'))
						# os.system('sort -k 3,3 -k 4,4n "%s" > "%s"'%(base+'.unsort', base))
						os.system('sort -k 3,3 -k 4,4n -S 50%'+' --parallel=%s -T "%s" "%s" > "%s"'%(nsubprocs, outdir, base, base+'.sorted'))
						os.system('mv "%s" "%s"'%(base+'.sorted', base))
						# os.system('rm -f "%s"'%(base+'.unsort')) # let user delete these
					
						# use samtools if available - avoids sorting the header
					indexFile(base, base+'.index')
					pipeit('done',1)
			
	if not ANYSAMPLES: pipeit('  No unique maps need indexing.',1)
	else: pipeit('- Done sorting and indexing unique maps. Manually remove *.unsort files if desired.', 1)

if INDEXMAPFILES or INDEXDEG:
	pipeit('- Indexing degenerate maps...', 1)
	
	peses = ['PE.', 'SE.']
	if PESECOMBO or UDPESECOMBO: peses = ['']
	
	ANYSAMPLES = False
	for sampleID in sorted(seqRuns.keys()):
		for pese in peses:
			bases = [btmapdir+sampleID+'.'+pese+'degenerate.map']
			afiles = [badReadDir+'%s_degenerate_redundantIDs.txt'%(sampleID)]
			tags = [sampleID+'.'+pese+'degenerate.map']
			if ALLELESPECIFIC:
				bases = [btmapdir+sampleID+'.%s.%sdegenerate.map'%('h1',pese), btmapdir+sampleID+'.%s.%sdegenerate.map'%('h2',pese)]
				afiles = [badReadDir+'%s_h1.degenerate_redundantIDs.txt'%(sampleID),badReadDir+'%s_h2.degenerate_redundantIDs.txt'%(sampleID)]
				tags = [sampleID+'.%s.%sdegenerate.map'%('h1',pese),sampleID+'.%s.%sdegenerate.map'%('h2',pese)]
				
			for z in range(len(bases)):
				base = bases[z]
				afile = afiles[z]
				tag = tags[z]
				if not os.access(base, os.F_OK): continue
				if not os.access(base+'.index', os.F_OK) or RESORT:
					pipeit('  %s...'%(getLabel(base)))
					ANYSAMPLES = True
					badIDs = {}
					if REDUNDANCYFILTER:
						pipeit('load PCR duplicates...')
						if not os.access(afile, os.F_OK): pipeit('Warning: redundancy filter ON but file not found: %s'%(afile),1)
						badIDs = loadPCRduplicates(afile)
					indexDegFile(base, tag, badIDs)
					pipeit('done.',1)
	if not ANYSAMPLES: pipeit('  No degenerate maps need indexing.',1)
	else: pipeit('- Done indexing degenerate maps.', 1)


if ESTIMATEGLOBALERROR:
	print ' - Estimating global sequencing error by phred scores...'
	
	# dictionary keyed by ascii character returning phred q-score
	# ascii 33 to 126 - full range
	fastqchars = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '\{', '\|', '\}', '~']
	A2Q = dict(zip(fastqchars, range(len(fastqchars))))
	
	
	log = lambda x: math.log(x, 10)
	# convert between phred score and p-value
	Q2P = lambda e: (e == 0 and 1-1e-16) or 10**(-e/10.) # avoid 1
	P2Q = lambda p: -10*log(p)
	
	# Compute the probability of no error
	A2P = dict( [ascii, Q2P(A2Q[ascii])] for ascii in A2Q.keys() )
	
	
	print 'Estimating global per-nucleotide sequencing error rate...'
	errorsum = 0
	ecount = 0
	# errors = []
	for mapfile in getFiles(btmapdir, include='unique'):
		print '- Loading', mapfile
		fh = open(mapfile)
		count = 1
		for line in fh:
			if count % 1000000 == 0: print '     reached', count
			count+=1
			
			try:
				lst = line[:-1].split('\t')
				chrom = lst[2]
				if chrom == '*': continue
				
				for x in map(lambda x: A2Q[x], lst[10]):
					errorsum += x
					ecount += 1
					# errors += [x]
				
			# exception if a line beginning with @
			except IndexError: pass
			except ValueError: pass
		fh.close()
	
	print 'Mean error over %s positions is'%(ecount), errorsum/float(ecount), 'or', Q2P(errorsum/float(ecount))
	# print 'Median error is', median(errors)
	print
	sys.exit()
# DONE


# # process SAM files via samtools
# # ------------------------------
# # print 'PART IIB: Convert SAM read maps--->BAM maps'
# for sampleID in sorted(seqRuns.keys()):
# 	for pese in ['PE', 'SE']:
# 		for ud in ['unique', 'degenerate']:
# 			# first check if we already processed the file
# 			if not os.access(btmapdir+sampleID+'.'+pese+'.'+ud+'.bam.bai', os.F_OK):
# 				pipeit(' - Samtools processing: %s\n'%(sampleID+'.'+pese+'.'+ud+'.map'))
# 				# convert SAM map file to BAM file
# 				# print 'call:', 'samtools view -bS "%s" > "%s"'%(btmapdir+sampleID+'.'+pese+'.'+ud+'.map', btmapdir+sampleID+'.'+pese+'.'+ud+'.bam')
# 				os.system('samtools view -bS "%s" > "%s"'%(btmapdir+sampleID+'.'+pese+'.'+ud+'.map', btmapdir+sampleID+'.'+pese+'.'+ud+'.bam'))
# 				# sort the BAM file
# 				pys.sort(btmapdir+sampleID+'.'+pese+'.'+ud+'.bam', btmapdir+sampleID+'.'+pese+'.'+ud)
# 				# index the BAM files for random retrieval
# 				pys.index(btmapdir+sampleID+'.'+pese+'.'+ud+'.bam')
# 				
# 				# delete the sam files
# 				if not KEEPSAM:
# 					os.system('rm -f "%s"'%(btmapdir+sampleID+'.'+pese+'.'+ud+'.map'))
# # print 'DONE PART IIB: SAM--->BAM conversion'



GR = None # we need the reference genome here...
flg = None
NEEDTOLOADREFERENCE = True


if not LIKCALLING and not PEAKCALLING and not JOINPIECES: sys.exit('DONE for now. Use the --lik and --post flags to compute likelihood and posterior probabilites for PTM binding.')

sampleIDs = sortseqruns

# Fill in prior probabilities
Pr = {}
for sampleID in sampleIDs:
	prior = []
	for run,lane in seqRuns[sampleID].keys():
		prior += [seqRuns[sampleID][(run,lane)][4]]
	if len(unique(prior)) > 1:
		sys.exit('Different prior probabilities are specified for the same sampleID. Check map file.')
	Pr[sampleID] = prior[0]

# associate sampleID with internal alias
sample2alias = {}
alias2sample = {}
inpsample = ''; nucsample = ''; ptmsample = ''

for sampleID in sampleIDs:
	aka = []
	for run,lane in seqRuns[sampleID].keys():
		# [refGenome, paired, path, prior, alias]
		aka += [seqRuns[sampleID][(run,lane)][5]]
	if len(unique(aka)) > 1:
		sys.exit('Different sample aliases are specified for the same sampleID. Check map file.')
	
	# equivalence
	if aka[0] == 'tf': aka[0] = 'nuc'

	sample2alias[sampleID] = aka[0]
	alias2sample[aka[0]] = sampleID
	

	if aka[0] == 'ptm': ptmsample = sampleID
	elif aka[0] == 'inp': inpsample = sampleID
	elif aka[0] == 'nuc': nucsample = sampleID
	elif aka[0] == 'atac' or aka[0] == 'dnase': nucsample = sampleID
	elif aka[0] == 'atacinp' or aka[0] == 'dnaseinp': inpsample = sampleID

if inpsample == 'input': inpsample == 'inp'

aliases = [sample2alias[x] for x in sampleIDs]

nchannels = len(Pr.keys())

# only deal with priors if using posterior
if PEAKCALLING:
	# if prior is not specified, use uninformative prior
	if nastr in Pr.values():
		pipeit('Prior probabilities not specified...using uninformative prior.')
		if nchannels == 2: Pr = dict(map(lambda k: (k,0.5), Pr.keys()))
		elif nchannels == 3: Pr = dict(map(lambda k: (k,1/3.), Pr.keys()))

	# fix input prior if only 2 channels are specified
	if nchannels == 1:
		if nucsample: Pr[inpsample] = 1-Pr[nucsample]
		else: sys.exit('Unknown model spec in adjusting prior probabilities: %s'%(Pr.keys()))
		
	elif nchannels == 2:
		if ptmsample and nucsample: Pr[nucsample] = 1-Pr[ptmsample]
		elif ptmsample and inpsample: Pr[inpsample] = 1-Pr[ptmsample]
		elif nucsample and inpsample: Pr[inpsample] = 1-Pr[nucsample]
		else: sys.exit('Unknown model spec in adjusting prior probabilities: %s'%(Pr.keys()))
		
	elif nchannels == 3:
		# need to make 3 changes
	
		# 1. Correct the nuc prior to be nuc - ptm
		# this makes the priors mutually exclusive (ie naked vs marked)
		try:
			Pr[nucsample] -= Pr[ptmsample]
		except KeyError:
			print Pr, 'hey', nucsample, ptmsample
			sys.exit('glimmrHet.py Error correcting prior probabilities: inp, nuc, ptm not properly specified in map file.')
	
		# 2. Normalize priors for relative antibody efficiency bias
		# idea is to increase the prior binding for the ptm and decrease for nuc
		# given the relative efficiency of both antibodies
		# while keeping the total prior probability Pnuc+Pptm constant
		K = (Pr[ptmsample]/Pr[nucsample])/Pxab
		x = (K*Pr[nucsample] - Pr[ptmsample])/(K+1)
		Pr[ptmsample] += x; Pr[nucsample] -= x
	
		# 3. Fix input if totals are screwy
		if Pr[inpsample] + Pr[nucsample] + Pr[ptmsample] > 1:
			Pr[inpsample] = 1 - (Pr[nucsample] + Pr[ptmsample])
		

lPerr = log(Perr) # log of the antibody efficiency
lPxab = log(Pxab) # log of the cross antibody efficiency
# lQ = math.log(Qsample,10)
lPr = dict( map(lambda k: (k,log(Pr[k])), Pr.keys()) )


# --------------------------
# fix sampleIDs to include allele specific analysis if specified
# haploidSampleIDs = sampleIDs
if ALLELESPECIFIC:
	haploidSampleIDs = []
	for sampleID in sampleIDs:
		haploidSampleIDs += [sampleID+'.h1']
		haploidSampleIDs += [sampleID+'.h2']
	sampleIDs = haploidSampleIDs
	
# first check that unique index files exist
if not PEAKCALLING:
	badsamples = []
	for sampleID in sampleIDs:
		if PESECOMBO:
			if not os.access('%s%s.unique.map'%(btmapdir,sampleID), os.F_OK):
				badsamples += [sampleID]
		else:
			if USESINGLETONS:
				if not os.access('%s%s.SE.unique.map'%(btmapdir,sampleID), os.F_OK):
					badsamples += [sampleID]
			else:
				if not os.access('%s%s.PE.unique.map'%(btmapdir,sampleID), os.F_OK):
					badsamples += [sampleID]
	
	if len(badsamples):
		pipeit('Error: the following maps must be indexed using --index:\n%s'%(', '.join(badsamples)),1)
	for s in badsamples: sampleIDs.pop(sampleIDs.index(s))
	
createdir(resdir); createdir(tmpdir)

theptmdir = resdir
if DRYRUN: theptmdir = tmpdir

completedir = zdir+'completed/'
createdir(completedir)


if LIKCALLING:
	
	# if not bfaref: sys.exit('No reference fasta file provided. Try again using --ref <PATH/TO/REFERENCE.FASTA>.')
	
	# SUBPROCESS CALLING for likelihood calculation
	# ---------------------------------------------
	jobs = []
	remstr = ''
	if REDOLIKCALLING: remstr += ' --relik'
	if STRATEGY == 'unique': remstr += ' --uniq'
	if STRATEGY == 'best': remstr += ' --best'
	if STRATEGY == 'bestno': remstr += ' --bestno'
	
	if LIPieces: remstr += ' --lip %s'%(LIPieces)
	cct = 1
	
	# if len(sampleIDs)==1 and nsubprocs > 1: xflag = 'subchrom'
	
	# EX: -rest chrY -x 10
	if len(seqbounds.keys()) == 1 and nsubprocs > 1: xflag = 'subchrom'
	# print 'TEST', xflag, 'threads', nsubprocs; sys.exit()

	# sampleIDs accounts for --ids if user supplied them
	NUM_IDS_TO_PROCESS = 0
	for sampleID in sampleIDs:
		dsampleID = ALLELESPECIFIC and '.'.join(sampleID.split('.')[:-1]) or sampleID
		if not bfaref: bfaref = seqRuns[dsampleID].values()[0][0]
		
		# may need to add suffix to reference file
		if not os.access(bfaref, os.F_OK):
			if os.access(bfaref+'.fa', os.F_OK):
				bfaref = bfaref+'.fa'
			elif os.access(bfaref+'.fasta', os.F_OK):
				bfaref = bfaref+'.fasta'
			else: sys.exit('Cannot access reference file listed in map file: %s'%(bfaref))
		
		
		print '- Computing likelihoods for', sampleID, 'against', getLabel(bfaref)
		
		# check if output file exists before rescheduling
		# print 'LOOKING', resdir+sampleID+'.mother.txt', os.access(resdir+sampleID+'.mother.txt', os.F_OK)
		# sys.exit()
		if REDOLIKCALLING or (not os.access(resdir+sampleID+'.mother'+boundslabel+'.txt', os.F_OK) and not os.access(resdir+sampleID+'.mother.txt', os.F_OK)):
			# set up jobs
			
			thislogdir = logdir+'%s.logs/'%(sampleID); createdir(thislogdir)
			
			# get the chroms from header in map file
			chrmlen = None
			if headerformat == 'sam':
				if PESECOMBO: chrmlen = getChromsFromHeader(btmapdir+sampleID+'.unique.map', fileformat=headerformat)
				elif USESINGLETONS: chrmlen = getChromsFromHeader(btmapdir+sampleID+'.SE.unique.map', fileformat=headerformat)
				else: chrmlen = getChromsFromHeader(btmapdir+sampleID+'.PE.unique.map', fileformat=headerformat)
			elif headerformat == 'fasta':
				chrmlen = getChromsFromHeader(bfaref, fileformat=headerformat)
			
			copydir = chrmlen.keys()
			# if user provides bounds, distribute them
			if len(seqbounds):
				pipeit(' - Distributing processes by user-supplied boundaries file, per sampleID...',1)
				copydir = seqbounds.keys()
			
			# now proceed to genotype each reference genomic position
			npieces = nsubprocs
			nchroms = len(copydir)
			
			if xflag == 'bychrom': npieces = nchroms
			elif xflag != 'subchrom' and npieces > nchroms: npieces = nchroms
			copydir.sort() # must be sorted!
			
			if nchroms == 0:
				sys.exit('Cannot find header. Check map directory: %s'%(btmapdir+sampleID))
				
			# ideally we want to split up evenly by summing the genome length
			# and then breaking up individual chromosomes
			# now that read maps are not loaded into memory we can do this easily - only factor is reference genome
			
			# add flags: --bychrom --bysample --auto (default)
			# --bychrom is --auto at the moment
			# no need for --by sample really
			
			
			# parallelize within chromosome for each chromosome for each sampleID
			if xflag == 'subchrom':
				
				ssstr = ''
				if sampleSubsample[dsampleID]: ssstr = ' --subsample %s'%(sampleSubsample[dsampleID])
				idstr = ' --sampleid %s '%(sampleID)
				cPID = sampleID
				asstr = ALLELESPECIFIC and '-as' or ''

				call = '%s %s --pid "%s" --m "%s" -i "%s" -o "%s" --key "%s" --d %s%s --rl %s --expphred %s --lambda %s --diameter %s %s %s %s'%(python, glpath, cPID, mapfile, btmapdir, outdir, zkey, maxnumalign, remstr, maxreadlen, expPhred, Lambda, fragDiameter4Sample[dsampleID], ssstr, idstr, asstr)
				
				if not HUGLY: call += ' --ref "%s"'%(bfaref)
				if USESWALLOW: call = '%s %s'%(swallowpath, call)
				if DEGREADSTREAM: call += ' --stream --buffer %s'%(bufferlimit)
				if not REDUNDANCYFILTER: call += ' --nopcr'
				else: call += ' --pcrdir "%s"'%(badReadDir)
				if PESECOMBO: call += ' --combo'
				if USESGE: call += ' --pythonpath "%s"'%(SGEPYTHONPATH)
				
				# continue to break up each chromosome into pieces
				# this is straightforward if we are not using SGE, so start with local multiprocessing
				# take each scaffold/chromosome at a time and break up into x pieces
				chromgroups = []
				for achr in copydir:
					# determine piece length
					# maxlen = chrmlen[achr]
					# maxlen = seqbounds[achr][1] != nastr and seqbounds[achr][1] or chrmlen[achr]
					
					# if user broke it down, use it
					if seqbounds[achr][0][1] != nastr:
						for x in seqbounds[achr]:
							chromgroups += [[achr]+[x[0],x[1]]]
					else:
						maxlen = chrmlen[achr]
						intervals = makeRange(0,maxlen,n=npieces+1)
						pieces = [[int(intervals[i-1]),int(intervals[i])-1] for i in xrange(1,len(intervals))]
						# just make sure first and last positions are correct
						pieces[0][0] = 0; pieces[-1][1] = maxlen
						# print 'P', achr, pieces, chrmlen[achr]
						for p in pieces: chromgroups += [[achr]+list(p)]
					
				pipeit('- Automatically distributing into --x %s jobs x %s chromosomes = %s jobs...'%(npieces,len(copydir), npieces*len(copydir)),1)
				
				# generate a job per chromosome interval
				for achr,start,stop in chromgroups:
					boundslabeltmp = '%s,%s,%s'%(achr,start,stop)
					boundstr = ' --rest %s'%(boundslabeltmp)
					logstr = ' >> "%s"'%(thislogdir+'%s_%s_%s.log'%(sampleID,boundslabeltmp,tmpstr))
					
					refstr = ''
					if HUGLY:
						if achr in bounds2Ref: 
							refstr = ' --ref "%s"'%(bounds2Ref[achr])
						else:
							refstr = ' --ref "%s"'%(bfaref)
					
					if VERBOSE: pipeit('   Call %s: %s'%(cct, call+refstr+boundstr+logstr),1)
					jobs += [call+refstr+boundstr+logstr]
					cct += 1
			
			else:
				pipeit('- Automatically distributing into -x %s jobs per sampleID (%s sampleIDs)...'%(npieces,len(sampleIDs)),1)
				
				restdir = zdir+'boundsfiles/'
				createdir(restdir)
				
				# split chromosomes into LIPieces loading groups
				chromgroups = []
				itemspergroup = math.ceil(nchroms / float(npieces))
				# if xflag == 'crosschrom': itemspergroup = nchroms # why?
				
				while len(copydir):
					currgrp = []
					while len(currgrp) < itemspergroup and len(copydir):
						currgrp += [copydir.pop(0)]
					chromgroups += [currgrp]
				
				# WHY ARE BOUNDS NOT TAKEN INTO ACCOUNT HERE
				
				# generate a job per chromosome group
				for chromgrp in chromgroups:
					# create a bounds file with these rows
					# (can't just concatenate and call because there is a character limit for unix shell commands)
					
					ssstr = ''
					if sampleSubsample[dsampleID]: ssstr = ' --subsample %s'%(sampleSubsample[dsampleID])
					idstr = ' --sampleid %s '%(sampleID)
					
					boundslabeltmp = ''
					cPID = sampleID
					boundstr = ''
					
					if nsubprocs > 1 and xflag != 'bychrom':
						boundslabeltmp = '%s--%s'%(chromgrp[0],chromgrp[-1])
						if len(chromgrp) == 1: boundslabeltmp = '%s'%(chromgrp[0])
						boundslabeltmp = boundslabeltmp.strip()
						boundsfiletmp = restdir+'%s.txt'%(boundslabeltmp)
						printList(chromgrp, file=boundsfiletmp)
						cPID += '.'+boundslabeltmp
						
						boundstr = '--rest "%s"'%(boundsfiletmp)
						if not len(boundsfiletmp): boundstr = ''
						
					elif nsubprocs > 1 and xflag == 'bychrom':
						# len(chromgrp) == 1 for this case
						boundstr = '--rest "%s"'%(chromgrp[0].strip())
						cPID += '.'+chromgrp[0].strip()
						
					elif boundsfile:
						# just reprint the string
						boundstr = '--rest "%s"'%(boundsfile)
						cPID += boundslabel
					
					# print xflag, "SDFL:KJSDLKJSFD', boundstr", boundstr, xflag
					# sys.exit()
					# before adding this as a job, just check to see if it is already completed
					if not REDOLIKCALLING and os.access(resdir+sampleID+'.mother.'+chromgrp[0].strip()+'.txt', os.F_OK): 
						pipeit('Done %s, '%(cPID))
						continue
					
					asstr = ALLELESPECIFIC and '-as' or ''
					
					call = '%s %s --pid "%s" --ref "%s" --m "%s" -i "%s" -o "%s" --key "%s" %s --d %s%s --rl %s --expphred %s --lambda %s --diameter %s %s %s %s'%(python, glpath,cPID, bfaref, mapfile, btmapdir, outdir, zkey, boundstr, maxnumalign, remstr, maxreadlen, expPhred, Lambda, fragDiameter4Sample[dsampleID], ssstr, idstr, asstr)
					
					if USESWALLOW: call = '%s %s'%(swallowpath, call)
					if DEGREADSTREAM: call += ' --stream --buffer %s'%(bufferlimit)
					if not REDUNDANCYFILTER: call += ' --nopcr'
					else: call += ' --pcrdir "%s"'%(badReadDir)
					if PESECOMBO: call += ' --combo'
					if USESGE: call += ' --pythonpath "%s"'%(SGEPYTHONPATH)
					
					logstr = ' > "%s"'%(thislogdir+'%s.log'%(cPID))
					if VERBOSE: pipeit('\n   Call %s: %s'%(cct, call+logstr),1)
					jobs += [call+logstr]
					cct += 1
	
			NUM_IDS_TO_PROCESS += 1
		else:
			pipeit('Already finished '+sampleID+'.mother'+boundslabel+'.txt',1)
	# if not len(jobs) and : sys.exit('- No jobs to schedule. Exiting...')
	
	if len(jobs):
		pipeit('- Now launching a total of %s jobs using %s parallel threads.'%(len(jobs),nsubprocs),1)
		pipeit('- Results will be saved to "%s".'%(loglikfile),1)
		printList(jobs,delim='\n\n',file=loglikfile) # save to log file
		
		# set up 
		# if os.access(completedir, os.F_OK) and len(getFiles(completedir)) > 0:
			# os.system('rm -f "%slikelihood"*'%(completedir))
		
		# this should be totally modular
		njobs = len(jobs) # initial job count
		nlaunched = 0
		running_jobs = []
		ncomplete = 0
		nkilled = 0
		if USESGE:
			while len(jobs):
				pipeit('- Scheduling job %s/%s (%s completed, %s killed | x%s processors)'%(nlaunched+1, njobs, ncomplete, nkilled, nsubprocs),1)
				# and we can launch another job
				# fetch the next job
				call = jobs.pop(0)
				# and execute it
				nlaunched += 1
				os.system('qsub -b y %s -cwd %s' % (SGEMEM, call))
				time.sleep(stagger_launch)
		else:
			while len(jobs):
				# as long as there are more jobs to be processed
				if nlaunched - ncomplete - nkilled < nsubprocs:
					pipeit('- Scheduling job %s/%s (%s completed, %s killed | x%s processors)'%(nlaunched+1, njobs, ncomplete, nkilled, nsubprocs),1)
					# and we can launch another job
					# fetch the next job
					call = jobs.pop(0)
					# and execute it
					nlaunched += 1
					job = subprocess.Popen(call, shell=True)
					running_jobs += [job]
					time.sleep(stagger_launch)
				else:
					# check if existing job has been killed
					# nkilled = 0
					# # check to see if an existing job has finished
					# ncomplete = len(getFiles(completedir, include='likelihood'))
					
					nkilled = 0
					marked = []
					for ji in xrange(len(running_jobs)): 
						if running_jobs[ji].poll() != None: 
							# print 'poll killed', rj[ji].poll(), len(running_jobs)
							ncomplete += 1
							marked += [running_jobs[ji]]
					for m in marked:
						running_jobs.pop(running_jobs.index(m))
					
					time.sleep(1)
	
			if not KEEPALIVE and njobs > 1 and nsubprocs > 1:
				sys.exit('\n- Likelihood processes have been launched. Check log_files/ directory to monitor progress.\n- Or run glimmrChecklist.py -i <output_dir>/log_files/<sampleID>.logs/')
				# When jobs have finished, execute glimmrHet.py with --merge[+] <sample1,...,samplek>\nto merge partial files into complete likelihood files.\nThen use the -p flag to compute posterior probabilities.
			else:
				for job in running_jobs: job.wait() # halt code until all jobs are finished.
	



# rejoin subchrom pieces
if JOINPIECES:# and xflag == 'subchrom':
	pipeit('- Joining subchrom likelihood pieces...',1)
	JOINFLAG = False
	for sampleID in sampleIDs:
		pipeit('  %s...'%(sampleID))
		# figure out the files
		tmpmap = btmapdir+sampleID+'.unique.map'
		if PESECOMBO: tmpmap = btmapdir+sampleID+'.unique.map'
		elif USESINGLETONS: tmpmap = btmapdir+sampleID+'.SE.unique.map'
		else: tmpmap = btmapdir+sampleID+'.PE.unique.map'

		copydir = getChromsFromHeader(tmpmap).keys()
		if len(seqbounds): copydir = seqbounds.keys()
		
		chromfiles = ['%s.mother.%s.txt'%(sampleID,achr) for achr in sorted(copydir)]
		finalfile = resdir+'%s.mother.txt'%(sampleID)
		if not os.access(finalfile, os.F_OK):
			for achr in sorted(copydir):
				pipeit('%s,'%(achr))
				motherfile = '%s.mother.%s.txt'%(sampleID,achr)
				chrpieces = getFiles(resdir, include='%s.mother.%s,'%(sampleID,achr), exclude=chromfiles)
				
				if os.access(resdir+motherfile, os.F_OK) and len(chrpieces)==0: continue
				elif len(chrpieces) == 0: 
					pipeit('No files with prefix %s.mother.%s,'%(sampleID,achr),1)
					JOINFLAG = True
					continue
				
				# put these into sorted order - by chrom then by position
				# parts = map(lambda x: [int(x.split('/')[-1].split('.')[2].split(',')[1]),x], chrpieces)
				# print 'PARTS', printTable(sorted(parts))
				parts = map(lambda x: [int(x.split('/')[-1].split('.')[-2].split(',')[1]),x], chrpieces)
				# print 'NEW PARTS', printTable(sorted(newparts))
				# print 'TESTING', sorted(parts) == sorted(newparts)
				# sys.exit()
				chrmother = resdir+'%s.mother.%s.txt'%(sampleID,achr)
				fh = open(chrmother,'w'); fh.close()
				
				tmpfiles = []
				for start,piece in sorted(parts):
					# print 'concatenating', piece.split('/')[-1]
					os.system('cat "%s" >> "%s"'%(piece,chrmother))
					tmpfiles += [piece]
				
				# if CLEANJOINPIECES:
				# do this always
				for tmpf in tmpfiles: os.system('rm -f "%s"'%(tmpf))
			
			if FORCEJOIN: JOINFLAG = False
			if len(copydir) > 1 and not JOINFLAG:
				pipeit('- Final merge...')
				# merge chromosomes into single file
				
				fh = open(finalfile,'w'); fh.close()
				tmpfiles = []
				for achr in sorted(copydir):
					chromfile = '%s%s.mother.%s.txt'%(resdir,sampleID,achr)
					os.system('cat "%s" >> "%s"'%(chromfile, finalfile))
					tmpfiles += [chromfile]
				if CLEANJOINPIECES:
					for tmpf in tmpfiles: os.system('rm -f "%s"'%(tmpf))
				pipeit('done',1)
		else:
			pipeit('already done',1)
	pipeit('- Finished merging the "--subchrom" mother likelihood files',1)
	pipeit('- Results saved to "%s"'%(resdir),1)
	if JOINFLAG:
		pipeit('- Did not merge chrom files into final mother file due to potential missing files. Use -forcejoin to override.')
	if not CLEANJOINPIECES:
		pipeit('- NB, you can now manually delete the partial likelihood files (samplex.mother.chrx.start,stop,+.txt).\nFor example: rm sample.mother.chrY,*.txt',1)
	
	
# DONE LIKELIHOOD COMPUTATIONS

if not PEAKCALLING: sys.exit('')

# pipeit('fix allele specific issues for posterior! what were these again?',1)

# NOW COMPUTE POSTERIOR PROBABILITIES AS SPECIFIED IN MAP FILE
# in order to multiprocess posteriors, need to nest the following in a function
# then we make calls to glimmrHet.py for posterior which calls the function
# here just launch jobs and then wait
# ------------------------------------------------------
jobs = []

# files that contain the resulting posterior probabilities and statistics
# outfile = resdir+'%s_model_calls%s.txt'%(keyid,boundslabel)
# keyfile = resdir+'%s.posterior%s.txt'%(keyid,boundslabel)
# keylogfile = postlogdir+'%s.posterior%s.log'%(keyid,boundslabel)

# create the jobs - loop through mapfile combos
# if user provided specific ids via --ids, then use them
# otherwise process the sample combinations in the map file
if not len(singlesampleid): # use these as the post combos
	# postcombos = [singlesampleid]
	cct = 1
	for combo in postcombos:
		# recall glimmr process with same arguments just add the ids
		# and the antibody cross efficiency parameter
		localPxab = 1.
		if len(combo) == 4: localPxab = combo.pop()
		# print 'combo', combo, 'pxab', localPxab
		idstr = ','.join(map(lambda x: '"'+x.strip()+'"', combo))
		call = ' '.join(sys.argv)+ ' --ids '+idstr+' --px %s'%(localPxab)
		call += ' --pid %s'%('-'.join(combo))
		# if boundsfile: call += ' --rest "%s"'%(boundsfile)
		# print 'HEY', idstr, 'combo', combo
		call += ' > "%s"'%(postlogdir+'%s.posterior.log'%(','.join(sorted(combo))))
		
		if USESWALLOW: call = '%s %s'%(swallowpath, call)
		
		print 'Call %s: %s'%(cct, call)
		jobs += [call]
		cct += 1
	pipeit('',1)
	
	# set up 
	# if os.access(completedir, os.F_OK) and len(getFiles(completedir)) > 0:
		# os.system('rm -f "%sposterior"*'%(completedir))
	
	njobs = len(jobs) # initial job count
	nlaunched = 0
	running_jobs = []
	ncomplete = 0
	nkilled = 0
	while len(jobs):
		# as long as there are more jobs to be processed
		if nlaunched - ncomplete - nkilled < nsubprocs:
			pipeit('- Scheduling job %s/%s (%s completed, %s killed | x%s processors)'%(nlaunched+1, njobs, ncomplete, nkilled, nsubprocs),1)
			
			# and we can launch another job
			# fetch the next job
			call = jobs.pop(0)
			# and execute it
			job = subprocess.Popen(call, shell=True)
			running_jobs += [job]
			nlaunched += 1
			time.sleep(1)
		else:
			# check if existing job has been killed
			nkilled = 0
			for job in running_jobs: 
				if job.poll() != None: nkilled += 1
			# if nkilled > 0: print '%s jobs killed'%(nkilled)
			# check to see if an existing job has finished
			ncomplete = len(getFiles(completedir, include='posterior'))
			time.sleep(10)
		
	# for job in running_jobs: job.wait() # halt code until all jobs are finished.
	# ----------------------------------
	
	if len(postcombos) and not KEEPALIVE:
		sys.exit('glimmrHet.py: Finished scheduling posterior probability jobs.')
	elif len(singlesampleid) > 3: 
		sys.exit('Specify only 3 channels using --ids or --sampleids flag to compute posterior probs.')

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# define the key sampleid
keyid = ''
if ptmsample: keyid = ptmsample
elif nucsample: keyid = nucsample
postfile = theptmdir+'%s.posterior%s.txt'%(keyid,boundslabel)

numHypotheses = 0
totalLoci = 0

if not os.access(postfile, os.F_OK) or REDOPEAKCALLING:
	covbiasdir = zdir+'coverage_bias/'
	createdir(covbiasdir)
	
	# shorthand
	SID = sampleIDs
	S2A = sample2alias
	A2S = alias2sample
	
	motherfiles = {}
	for sampleID in sampleIDs:
		motherfile = theptmdir+sampleID+'.mother'
		if boundsfile: motherfile += boundslabel+'.txt'
		else: motherfile += '.txt'
		motherfiles[sampleID] = motherfile
	
	# modeltype selection
	mdec = {1:'inp,nuc', 2:'nuc,ptm', 3:'inp,ptm', 4:'inp,nuc,ptm', 6:'atac,atacinp', 7:'dnase,dnaseinp'}
	modeltype = 0
	
	if nchannels == 2 and 'inp' in A2S and A2S['inp'] in SID and 'nuc' in A2S and A2S['nuc'] in SID: modeltype = 1
	elif nchannels == 2 and 'nuc' in A2S and A2S['nuc'] in SID and 'ptm' in A2S and A2S['ptm'] in SID: modeltype = 2
	elif nchannels == 2 and 'inp' in A2S and A2S['inp'] and 'ptm' in A2S and A2S['ptm'] in SID: modeltype = 3
	elif nchannels == 3: modeltype = 4
	
	elif nchannels == 2 and 'atac' in A2S and A2S['atac'] in SID and 'atacinp' in A2S and A2S['atacinp'] in SID: modeltype = 6
	elif nchannels == 2 and 'dnase' in A2S and A2S['dnase'] in SID and 'dnaseinp' in A2S and A2S['dnaseinp'] in SID: modeltype = 7
	elif nchannels == 1 and 'atac' in alias2sample and alias2sample['atac'] in sampleIDs:
		modeltype = 6
		nchannels = 2
		inpsample = alias2sample['atac']+'.shuffle'
		sampleIDs += [inpsample]
		motherfile = theptmdir+inpsample+'.mother'
		if boundsfile: motherfile += boundslabel+'.txt'
		else: motherfile += '.txt'
		motherfiles[inpsample] = motherfile
		
		sampleReadLengths[inpsample] = sampleReadLengths[nucsample]
		fragDiameter4Sample[inpsample] = fragDiameter4Sample[nucsample]
		
		Pr[inpsample] = 1-Pr[nucsample]
		sample2alias[inpsample] = 'atacinp'
		alias2sample['atacinp'] = inpsample
	elif nchannels == 1 and 'dnase' in alias2sample and alias2sample['dnase'] in sampleIDs:
		modeltype = 7
		nchannels = 2
		inpsample = alias2sample['dnase']+'.shuffle'
		sampleIDs += [inpsample]
		motherfile = theptmdir+inpsample+'.mother'
		if boundsfile: motherfile += boundslabel+'.txt'
		else: motherfile += '.txt'
		motherfiles[inpsample] = motherfile
		
		sampleReadLengths[inpsample] = sampleReadLengths[nucsample]
		fragDiameter4Sample[inpsample] = fragDiameter4Sample[nucsample]
		
		Pr[inpsample] = 1-Pr[nucsample]
		sample2alias[inpsample] = 'dnaseinp'
		alias2sample['dnaseinp'] = inpsample
		
	
	# shorthand
	SID = sampleIDs
	S2A = sample2alias
	
	A2S = alias2sample
	
	if modeltype == 0: 
		sys.exit('\nError: Model specification in map file is incorrect. Supply 2 or 3 ids with --ids.\nProvided: %s'%(sampleIDs))
	
	
	
	nsids = len(SID)
	rlc = range(nsids)
	
	inpidx = None;
	if inpsample: inpidx = sampleIDs.index(inpsample)
	nucidx = None
	if nucsample: nucidx = sampleIDs.index(nucsample)
	ptmidx = None
	if ptmsample: ptmidx = sampleIDs.index(ptmsample)
	keyidx = sampleIDs.index(keyid)
	
	# files that contain the resulting posterior probabilities and statistics
	# outfile = resdir+'%s_model_calls%s.txt'%(keyid,boundslabel)
	keyfile = resdir+'%s.posterior%s.txt'%(keyid,boundslabel)
	keylogfile = postlogdir+'%s.posterior%s.log'%(keyid,boundslabel)
	
	outlogfh = open(keylogfile,'a')
	
	pipeit('- Model choice: %s (%s)'%(modeltype,mdec[modeltype]),1)
	pipeit('- Model choice: %s (%s)'%(modeltype,mdec[modeltype]),1,outlogfh)
	
	# Confirm that the likelihood files exist
	HAVEALLLIKS = True
	for s in SID:
		if not os.access(motherfiles[s], os.F_OK):
			HAVEALLLIKS = False
			if modeltype == 6:
				sys.stderr.write('  - WARNING: Cannot access shuffled likelihood file for ATAC-seq sample: %s.\n    Use glimmrShuffleLikelihood.py to generate it.\n'%(motherfiles[s])); sys.stderr.flush()
			elif modeltype == 7:
				sys.stderr.write('  - WARNING: Cannot access shuffled likelihood file for DNase-seq sample: %s.\n    Use glimmrShuffleLikelihood.py to generate it.\n'%(motherfiles[s])); sys.stderr.flush()
			else:
				sys.stderr.write('  - WARNING: Cannot access likelihood file %s\n'%(motherfiles[s])); sys.stderr.flush()
			
		else: pipeit('    Found likelihood file %s'%(motherfiles[s].split('/')[-1]),1)
	if not HAVEALLLIKS: sys.exit('Error: missing likelihood files. Call glimmrHet.py using the -s flag.')
	
	# Check if reference genome is already in memory
	NEEDTOLOADREFERENCE = True
	if GR != None: NEEDTOLOADREFERENCE = False
	
	flg = None
	if NEEDTOLOADREFERENCE:
		if not bfaref:
			# choose as reference the entry in the map file corresponding to the ptm sample
			bfaref = seqRuns[keyid].values()[0][0]
		
		# may need to add suffix
		if not os.access(bfaref, os.F_OK):
			if os.access(bfaref+'.fa', os.F_OK):
				bfaref = bfaref+'.fa'
			elif os.access(bfaref+'.fasta', os.F_OK):
				bfaref = bfaref+'.fasta'
			else: sys.exit('Cannot access reference file listed in map file: %s'%(bfaref))
		
		sys.stdout.write('- Loading reference genome: %s '%(bfaref))
		sys.stdout.flush()
		GR = readFasta(bfaref, TOUPPER=True)
		flg = sum(map(lambda x: len(x['seq']), GR.values()))
		sys.stdout.write('N=%s loci, M=%s chromosomes\n'%(flg,len(GR)))
		
		pipeit('- Reference genome: %s N=%s loci, M=%s chromosomes'%(bfaref,flg,len(GR)),1,outlogfh)
		
	# get chromosome lengths
	chromlengths = dict( map(lambda k: (k, len(GR[k]['seq'])), GR) )
	
	# Set up genotyping bounds
	GRNA = {}
	bounds = seqbounds
	if len(bounds.keys()):
		flg = 0
		print '- We have genotyping bounds...'
		# print bounds
		
		null = [nastr, 'NA']
		for tmpk in GR.keys():
			master = GR[tmpk]['seq'][:]
			
			SEEN = False
			for k in [tmpk, tmpk.lower(), tmpk.upper()]:
				if k in bounds and not SEEN:
					SEEN = True
					lastpos = 0
					# may be many positions of interest
					for low,high,flag in sorted(bounds[k]):
						if flag == '+':
							# remove positions not specified in a range
						
							if low not in null and high not in null:
								# we have a range
								tmps = master
								tmp = ''.join(map(lambda x: '*', tmps[0:low]))
								tmp += tmps[low:high+1]
								tmp += ''.join(map(lambda x: '*', tmps[high+1:len(tmps)]))
								master = tmp
							elif low not in null:
								# we have a single position of interest
								# since we've sorted these positions, 
						
								# erase between lastpos and low
								tmp = master[0:lastpos+1]
								tmp += ''.join(map(lambda x: '*', master[lastpos+1:low]))
								tmp += master[low]
								tmp += master[low+1:len(master)]
								master = tmp
								lastpos = low
							# else just leave the entire chromosome intact
						elif flag == '-':
							# set ranges to NA
							if low not in null and high not in null:
								tmps = master
								tmp = tmps[0:low]
								tmp += ''.join(map(lambda x: '*', tmps[low:high+1]))
								tmp += tmps[high+1:len(tmps)]
								master = tmp
						
					# now I can finish erasing the end
					if lastpos:
						master = master[0:lastpos] + ''.join(map(lambda x: '*', master[lastpos:len(master)]))
				
					# write back
					GRNA[tmpk] = {'id':tmpk, 'seq':master}
					flg += len(master)
	else: GRNA = GR
	del GR # no longer need the reference
	
	
	pipeit('- Computing posterior probabilities...',1)
	pipeit('- Computing posterior probabilities...',1,outlogfh)
	
	# Initialize variables for posterior computation
	# ==============================================
	# define the PTM bias coefficient
	# relative weight given to PTM due to the fact that
	# nucleosome reads may contain PTM
	Zptm = None; lZptm = None
	if nucsample and ptmsample:
		# what is the proportion of ptm among any histone?
		Zptm = max(Pr[ptmsample]/(Pr[nucsample]+Pr[ptmsample]), Perr)
		lZptm = log(Zptm)
		# print 'test', Perr, Pr[ptmsample], Pr[nucsample], Pr[ptmsample]/Pr[nucsample], Zptm
	
	
	# display the calibrated parameters
	pipeit('- Pxef = %s | Perr = %s | PTM bias = %s'%(Pxab, Perr, Zptm), 1)
	pipeit('- Pxef = %s | Perr = %s | PTM bias = %s'%(Pxab, Perr, Zptm), 1,outlogfh)
	pipeit('- Priors',1)
	pipeit('- Priors',1,outlogfh)
	
	for k in sorted(Pr.keys()): 
		pipeit('  %s: %.3f'%(k,Pr[k]),1)
		pipeit('  %s: %.3f'%(k,Pr[k]),1,outlogfh)
	pipeit('',1)
	
	
	
	bS2A = {}
	for anid in allSeqRuns.keys():
		for lst in allSeqRuns[anid].values(): bS2A[anid] = lst[-1]
	
	# between-scale normalization
	if DOBNORM and not len(betweenSID):
		# betweenBase is a list of types specified by -bnorm at the CLI (so major,minor,male)
		if swapBase:
			for anid in SID: betweenSID += map(lambda x: anid.replace(swapBase,x), betweenBase)
		
		# figure out which betweenBase is found in the keyid sample
		else:
			# print 'keyid', keyid
			# print 'base', betweenBase
			# whichsample = filter(lambda x: keyid.count(x), betweenBase)
			# whichsample = filter(lambda x: x in keyid, betweenBase)
			
			# find the longest length true match
			tmp = filter(lambda x: x in keyid, betweenBase)
			whichsample = []
			for x in tmp:
				if not len(whichsample):
					whichsample = [x]
				elif len(x) > len(whichsample[0]):
					whichsample = [x]
				else:
					whichsample += [x]
			
			if len(whichsample) != 1: sys.exit('Error: auto between-sample norm specified but keyid not matched: %s vs %s, result = %s'%(keyid, betweenBase,whichsample))
			
			# sys.exit('have %s found %s'%(keyid,whichsample))
			
			swapBase = whichsample[0]
			betweenBase.pop(betweenBase.index(swapBase))
			# print 'ok', betweenBase
			for anid in SID: betweenSID += map(lambda x: anid.replace(swapBase,x), betweenBase)
	
	# print 'test!', betweenSID
	
	if DOBNORM: 
		if len(betweenSID):
			pipeit('- Samples for between-sample coverage normalization: %s'%(betweenSID), 1)
			pipeit('- Samples for between-sample coverage normalization: %s'%(betweenSID), 1,outlogfh)
		else:
			sys.exit('Error: between-sample norm specified but no samples specified.')
	
	
	# READ LENGTH BIAS NORMALIZATION
	# ==============================
	# likelihoods are computed for each sample, which may have different read lengths
	# so need to correct the log likelihood scores at each position by the ratio of read lengths
	# between the shortest read length and sample read length
	# get smallest read length over samples
	mrl = min([sampleReadLengths[k] for k in SID])
	# compute ratio of smallest read length (largest likelihood) to sample length (larger likelihood)
	Krl = dict(map(lambda k: (k,mrl/float(sampleReadLengths[k])), SID))
	pipeit('- Read length corrections: %s'%(Krl),1)
	pipeit('- Read length corrections: %s'%(Krl),1,outlogfh)
	
	# DIAMETER BIAS NORMALIZATION
	# ===========================
	mdiam = min([fragDiameter4Sample[k] for k in SID]) # smallest diameter
	Kdiam = dict(map(lambda k: (k,mdiam/float(fragDiameter4Sample[k])), SID))
	pipeit('- Fragment diameter corrections: %s'%(Kdiam),1)
	pipeit('- Fragment diameter corrections: %s'%(Kdiam),1,outlogfh)
	
	# Composite length + bias correction (reduces computations)
	Kll = dict(map(lambda k: (k,Krl[k]*Kdiam[k]), SID))
	pipeit('- Composite correction: %s'%(Kll),1)
	pipeit('- Composite correction: %s'%(Kll),1,outlogfh)
	# convert to list for quick access
	Kll = [Kll[anid] for anid in SID]
	
	
	# COVERAGE BIAS NORMALIZATION
	# ===========================
	# Count total reads per channel to normalize for coverage bias
	pipeit('- Counting total mapped reads...',1)
	pipeit('- Counting total mapped reads...',1,outlogfh)
	
	# This should be done per chromosome (or better, per megabase or so) to avoid scale bias
	# These values are generated during mapping
	totalReads = {}; betTotalReads = {}
	ssubtmp = {}
	if os.access(crdfile, os.F_OK):
		pipeit('- Loading %s'%(crdfile.split('/')[-1]),1)
		pipeit('- Loading %s'%(crdfile.split('/')[-1]),1,outlogfh)
		d,r,c = readTable(crdfile, header=1)
		
		hashapdipcol = False
		countstartidx = 1
		if len(c) == 6: hashapdipcol = True; countstartidx = 2
		
		# print 'test!', hashapdipcol,countstartidx
		
		for i in range(len(d)):
			samplename = r[i]
			# depends on whether the file has hap/dip column
			# if not, assumes haploid
			if hashapdipcol: currscaf = d[i][1]
			else: currscaf = d[i][0]
			
			if samplename in SID:
				# convert this to inp, nuc, or ptm
				if S2A[samplename] not in totalReads: totalReads[S2A[samplename]] = {currscaf:0}
				
				# get the total number of reads for this row
				trchr = 0 # total read count for this scaffold/chromosome
				
				# breakdown of read counts
				countUSE, countUPE, countDSE, countDPE = map(floatna, d[i][countstartidx:])
				
				if not UNIQ: # use U and D
					if USESINGLETONS:
						# use SE and PE
						trchr = sumna([countUSE, countUPE, countDSE, countDPE])
					else:
						# only PE
						trchr = sumna([countUPE, countDPE])
				else: # use U only
					if USESINGLETONS:
						# use SE and PE
						trchr = sumna([countUSE, countUPE])
					else:
						# only PE
						trchr = countUPE != nastr and countUPE or 0
						
				# store and correct for sample-specific read subsampling (chrom file doesn't account for subsampling)
				totalReads[S2A[samplename]][currscaf] = trchr * sampleSubsample[samplename]
				ssubtmp[S2A[samplename]] = sampleSubsample[samplename]
			
			# store between sample IDs as well
			elif samplename in betweenSID:
				# convert this to inp, nuc, or ptm
				if bS2A[samplename] not in betTotalReads: betTotalReads[bS2A[samplename]] = {currscaf:0}
				
				# get the total number of reads for this row
				trchr = 0 # total read count for this scaffold/chromosome
				
				# breakdown of read counts
				try:
					countUSE, countUPE, countDSE, countDPE = map(floatna, d[i][2:])
				except ValueError:
					# if this happens it is because the file is from older version of glimmr without the hap/dip column
					countUSE, countUPE, countDSE, countDPE = map(floatna, d[i][1:])
				
				if not UNIQ: # use U and D
					if USESINGLETONS:
						# use SE and PE
						trchr = sumna([countUSE, countUPE, countDSE, countDPE])
					else:
						# only PE
						trchr = sumna([countUPE, countDPE])
				else: # use U only
					if USESINGLETONS:
						# use SE and PE
						trchr = sumna([countUSE, countUPE])
					else:
						# only PE
						trchr = countUPE != nastr and countUPE or 0
						
				# store and correct for sample-specific read subsampling (chrom file doesn't account for subsampling)
				# note we are ACCUMULATING READS FOR BACKGROUND, UNLIKE THE SECTION ABOVE
				try: betTotalReads[bS2A[samplename]][currscaf] += trchr * sampleSubsample[samplename]
				except KeyError: betTotalReads[bS2A[samplename]][currscaf] = trchr * sampleSubsample[samplename]
	else: print >> sys.stderr, 'Warning: could not access file: %s'%(crdfile)
	
	# link atac shuffle
	for key in ['atac','dnase']:
		if key in totalReads and '%sinp'%(key) not in totalReads: totalReads['%sinp'%(key)] = totalReads[key]
		if key in ssubtmp and '%sinp'%(key) not in ssubtmp: ssubtmp['%sinp'%(key)] = ssubtmp[key]
		if key in betTotalReads: betTotalReads['%sinp'%(key)] = betTotalReads[key]

	pipeit('- Subsample rate: %s'%(ssubtmp),1)
	pipeit('- Subsample rate: %s'%(ssubtmp),1,outlogfh)
	
	# Compute total reads for each sample
	# ===================================
	GTReads = {} # = {'inp':0, 'nuc':0, 'ptm':0}
	for alias in totalReads.keys(): GTReads[alias] = sum(totalReads[alias].values())
	if len(GTReads) == 0: 
		pipeit('Warning: no coverage entries found. Ensure sampleIDs match. Proceeding with no coverage correction.',1)
		pipeit('Warning: no coverage entries found. Ensure sampleIDs match. Proceeding with no coverage correction.',1,outlogfh)
	else: 
		pipeit('- GTreads: %s'%(GTReads),1)
		pipeit('- GTreads: %s'%(GTReads),1,outlogfh)
	
	# Compute total reads for each background sample
	# ==============================================
	bGTReads = {} # = {'inp':0, 'nuc':0, 'ptm':0}
	for alias in betTotalReads.keys(): bGTReads[alias] = sum(betTotalReads[alias].values())
	if DOBNORM:
		if len(bGTReads) == 0: 
			pipeit('Warning: no coverage entries found. Ensure sampleIDs match. Proceeding with no coverage correction.',1)
		else: pipeit('- bGTreads: %s'%(bGTReads),1)
	
	
	# Estimate global coverage bias: probabilistic weights
	
	# why 2 times the total?
	# ttr = 2*float(sum(GTReads.values()))
	ttr = float(sum(GTReads.values()))
	Pgbias = {}
	pri = prn = prp = None
	
	if modeltype == 1:
		if ttr > 0:
			pri = ttr/GTReads['inp']
			prn = ttr/GTReads['nuc']
			Pgbias = {'inp':pri/(pri+prn), 'nuc':prn/(pri+prn)}
		else:
			Pgbias = {'inp':1/2., 'nuc':1/2.}
	elif modeltype == 2:
		if ttr > 0:
			prn = ttr/GTReads['nuc']
			prp = ttr/GTReads['ptm']
			Pgbias = {'nuc':prn/(prn+prp), 'ptm':prp/(prn+prp)}
		else:
			Pgbias = {'nuc':1/2., 'ptm':1/2.}
	elif modeltype == 3:
		if ttr > 0:
			pri = ttr/GTReads['inp']
			prp = ttr/GTReads['ptm']
			Pgbias = {'inp':pri/(pri+prp), 'ptm':prp/(pri+prp)}
		else:
			Pgbias = {'inp':1/2., 'ptm':1/2.}
	elif modeltype == 4:
		if ttr > 0:
			pri = ttr/GTReads['inp']
			prn = ttr/GTReads['nuc']
			prp = ttr/GTReads['ptm']
			Pgbias = {'inp':pri/(pri+prn+prp), 'nuc':prn/(pri+prn+prp), 'ptm':prp/(pri+prn+prp)}
		else: 
			Pgbias = {'inp':1/3., 'nuc':1/3., 'ptm':1/3.}
	elif modeltype == 5:
		if ttr > 0:
			prn = ttr/GTReads['atac']
			Pgbias = {'atac':1.}
		else:
			Pgbias = {'atac':1.}
	elif modeltype == 6:
		if ttr > 0:
			pri = ttr/GTReads['atacinp']
			prn = ttr/GTReads['atac']
			Pgbias = {'atacinp':pri/(pri+prn), 'atac':prn/(pri+prn)}
		else:
			Pgbias = {'atacinp':1/2., 'atac':1/2.}
	elif modeltype == 7:
		if ttr > 0:
			pri = ttr/GTReads['dnaseinp']
			prn = ttr/GTReads['dnase']
			Pgbias = {'dnaseinp':pri/(pri+prn), 'dnase':prn/(pri+prn)}
		else:
			Pgbias = {'dnaseinp':1/2., 'dnase':1/2.}
	
	# Scale the global within-sample read bias correction with between-sample correction
	# ==================================================================================
	# for each alias inp, nuc, ptm, we correct the alias for the sample of interest
	# by creating a weight that is the inverse of the proportion of this sample's reads
	# compared to all reads
	
	pipeit('- Global coverage correction',1)
	pipeit('- Global coverage correction',1,outlogfh)
	pipeit('  Within-sample coverage weights:  %s'%(' '.join(['%s=%.3f'%(x,y) for x,y in Pgbias.items()])), 1)
	pipeit('  Within-sample coverage weights:  %s'%(' '.join(['%s=%.3f'%(x,y) for x,y in Pgbias.items()])), 1,outlogfh)
	
	# Compute global between-sample coverage weights
	# ==============================================
	if DOBNORM:
		# if no other samples are given, defaults to 1 which is no change
		betWeight = {}
		try:
			bias = {'inp':1,'nuc':1,'ptm':1} # bias towards background?
			totaltmp = dict([ (aka,GTReads[aka]+bias[aka]*bGTReads[aka]) for aka in GTReads.keys()])
			for aka in totaltmp.keys():
				# weight is the fraction of total coverage from other samples
				# betWeight[aka] = (totaltmp[aka] - GTReads[aka])/float(totaltmp[aka])
				betWeight[aka] = bGTReads[aka]/float(totaltmp[aka])
		except KeyError:
			pipeit('- Not all samples have data for this protein/mark. Defaulting to no bewteen-sample normaliztion.',1)
			pipeit('- Not all samples have data for this protein/mark. Defaulting to no bewteen-sample normaliztion.',1,outlogfh)
			# this will occur if not all conditions have data for this mark/sample
			for aka in GTReads.keys():
				betWeight[aka] = 1/float(len(GTReads.keys()))
		pipeit('  Between-sample weights: %s'%(' '.join(['%s=%.3f'%(x,y) for x,y in betWeight.items()])), 1)
		pipeit('  Between-sample weights: %s'%(' '.join(['%s=%.3f'%(x,y) for x,y in betWeight.items()])), 1,outlogfh)
	else:
		pipeit('- No between-sample normalization', 1)
		pipeit('- No between-sample normalization', 1,outlogfh)
	
	
	totalLoci = 0     # count number of reference positions passed
	numHypotheses = 0 # count up number of hypotheses we make
	numKeyHypotheses = 0 # count up number of hypotheses where MAP indicates keyid
	missingLik = dict( map(lambda x: (x,0), sampleIDs) )
	
	# now proceed to genotype each reference genomic position
	chromosomesOfInterest = sorted(GRNA.keys())
	nchroms = len(chromosomesOfInterest)
	
	# open handles to each of the sample likelihood files
	fh = dict( map(lambda x: (x,open(motherfiles[x])), sampleIDs) )
	
	
	# If multiple ptm samples are called against the same input and nuc data we don't need to resave the input and nuc data...
	# outfh = open(outfile, 'w')
	outfh = open(keyfile, 'w')
	
	
	# pseudolog = log(.99) # impute with a perfect read (no mismatches)
	pseudolog = log(1.0 - 10**(-expPhred/10.)) # expected likelihood assuming no mismsatches
	pseudocount = 1 # for missing data
	# want to set pseudocount to the genomewide average coverage (for input only)
	
	# keep track of current position for each likelihood file, per chromosome
	pos = dict( map(lambda x: (x,-1), sampleIDs) )
	chrm = dict( map(lambda x: (x,None), sampleIDs) )
	# keep track of current line for each likelihood file
	line = dict( map(lambda x: (x,None), sampleIDs) )
	
	matched = dict( map(lambda x: (x,False), sampleIDs) )
	
	# save the relative bias probabilities to this directory
	covbiasfile = covbiasdir+'coverage_bias_%s.txt'%(','.join(sorted(A2S.values())))
	cbfh = open(covbiasfile,'w'); cbfh.close()
	record('chrom\tinp\tnuc\tptm', covbiasfile, 'a')
	
	# # need to save global priors as separate variable here
	# gPr = Pr
	
	# Loop over chromosome groups
	# tmpidx = chromosomesOfInterest.index('scaffold244')
	# for chrom in chromosomesOfInterest[tmpidx:]:
	for chrom in chromosomesOfInterest:
		chromlen = chromlengths[chrom]
		
		# compute scaffold specific read bias correction
		# this is an unweighted average between global bias and chromosome-specific bias
		Pbias = {}
		ttr = 0
		try: ttr = float(sum([totalReads[k][chrom] for k in totalReads.keys()]))
		except KeyError: pass
		
		if ttr > 0:
			try:
				if modeltype == 1:
					pri = ttr/totalReads['inp'][chrom]; prn = ttr/totalReads['nuc'][chrom]
					Pbias = {'inp':pri/(pri+prn), 'nuc':prn/(pri+prn)}
					# Combine this scaffold-bias with global-bias probability
					Pbias = {'inp':(Pbias['inp']+Pgbias['inp']), 'nuc':(Pbias['nuc']+Pgbias['nuc'])}
					# Rescale to probability
					Pt = float(sum(Pbias.values()))
					Pbias = {'inp':Pbias['inp']/Pt, 'nuc':Pbias['nuc']/Pt}
					# print 'global bias    : inp: %.3f, nuc: %.3f'%(Pgbias['inp'], Pgbias['nuc'])
					# print 'local bias     : inp: %.3f, nuc: %.3f'%(totalReads['inp'][chrom], totalReads['nuc'][chrom])
					pipeit('- %s cov. bias: inp: %.3f, nuc: %.3f'%(chrom, Pbias['inp'], Pbias['nuc']),1)
					pipeit('- %s cov. bias: inp: %.3f, nuc: %.3f'%(chrom, Pbias['inp'], Pbias['nuc']),1,outlogfh)
					
				elif modeltype == 2:
					prn = ttr/totalReads['nuc'][chrom]; prp = ttr/totalReads['ptm'][chrom]
					Pbias = {'nuc':prn/(prn+prp), 'ptm':prp/(prn+prp)}
					Pbias = {'nuc':Pbias['nuc']+Pgbias['nuc'], 'ptm':Pbias['ptm']+Pgbias['ptm']}
					Pt = float(sum(Pbias.values()))
					Pbias = {'nuc':Pbias['nuc']/Pt, 'ptm':Pbias['ptm']/Pt}
					pipeit('- %s cov. bias: nuc: %.3f, ptm: %.3f'%(chrom, Pbias['nuc'], Pbias['ptm']),1)
					pipeit('- %s cov. bias: nuc: %.3f, ptm: %.3f'%(chrom, Pbias['nuc'], Pbias['ptm']),1,outlogfh)
					
				elif modeltype == 3:
					pri = ttr/totalReads['inp'][chrom]; prp = ttr/totalReads['ptm'][chrom]
					Pbias = {'inp':pri/(pri+prp), 'ptm':prp/(pri+prp)}
					Pbias = {'inp':Pbias['inp']+Pgbias['inp'], 'ptm':Pbias['ptm']+Pgbias['ptm']}
					Pt = float(sum(Pbias.values()))
					Pbias = {'inp':Pbias['inp']/Pt, 'ptm':Pbias['ptm']/Pt}
					pipeit('- %s cov. bias: inp: %.3f, ptm: %.3f'%(chrom, Pbias['inp'], Pbias['ptm']),1)
					pipeit('- %s cov. bias: inp: %.3f, ptm: %.3f'%(chrom, Pbias['inp'], Pbias['ptm']),1,outlogfh)
					
				elif modeltype == 4:
					pri = ttr/totalReads['inp'][chrom]; prn = ttr/totalReads['nuc'][chrom]; prp = ttr/totalReads['ptm'][chrom]
					Pbias = {'inp':pri/(pri+prn+prp), 'nuc':prn/(pri+prn+prp), 'ptm':prp/(pri+prn+prp)}
					Pbias = {'inp':Pbias['inp']+Pgbias['inp'], 'nuc':Pbias['nuc']+Pgbias['nuc'], 'ptm':Pbias['ptm']+Pgbias['ptm']}
					Pt = float(sum(Pbias.values()))
					Pbias = {'inp':Pbias['inp']/Pt, 'nuc':Pbias['nuc']/Pt, 'ptm':Pbias['ptm']/Pt}
					pipeit('- %s cov. bias: inp: %.3f, nuc: %.3f, ptm: %.3f'%(chrom, Pbias['inp'], Pbias['nuc'], Pbias['ptm']),1)
					pipeit('- %s cov. bias: inp: %.3f, nuc: %.3f, ptm: %.3f'%(chrom, Pbias['inp'], Pbias['nuc'], Pbias['ptm']),1,outlogfh)
			
				elif modeltype == 5:
					# only 1 sample so no correction is needed
					Pbias = {'atac':1}
					# Combine this scaffold-bias with global-bias probability
					Pbias = {'atac':(Pbias['atac']+Pgbias['atac'])}
					# Rescale to probability
					Pt = float(sum(Pbias.values()))
					Pbias = {'atac':Pbias['atac']/Pt}
					pipeit('- %s cov. bias: atac: %.3f'%(chrom, Pbias['atac']),1)
					pipeit('- %s cov. bias: atac: %.3f'%(chrom, Pbias['atac']),1,outlogfh)
				
				elif modeltype == 6:
					pri = ttr/totalReads['atacinp'][chrom]; prn = ttr/totalReads['atac'][chrom]
					Pbias = {'atacinp':pri/(pri+prn), 'atac':prn/(pri+prn)}
					# Combine this scaffold-bias with global-bias probability
					Pbias = {'atacinp':(Pbias['atacinp']+Pgbias['atacinp']), 'atac':(Pbias['atac']+Pgbias['atac'])}
					# Rescale to probability
					Pt = float(sum(Pbias.values()))
					Pbias = {'atacinp':Pbias['atacinp']/Pt, 'atac':Pbias['atac']/Pt}
					# print 'global bias    : inp: %.3f, nuc: %.3f'%(Pgbias['inp'], Pgbias['nuc'])
					# print 'local bias     : inp: %.3f, nuc: %.3f'%(totalReads['inp'][chrom], totalReads['nuc'][chrom])
					pipeit('- %s cov. bias: atacinp: %.3f, atac: %.3f'%(chrom, Pbias['atacinp'], Pbias['atac']),1)
					pipeit('- %s cov. bias: atacinp: %.3f, atac: %.3f'%(chrom, Pbias['atacinp'], Pbias['atac']),1,outlogfh)

				elif modeltype == 7:
					pri = ttr/totalReads['dnaseinp'][chrom]; prn = ttr/totalReads['dnase'][chrom]
					Pbias = {'dnaseinp':pri/(pri+prn), 'dnase':prn/(pri+prn)}
					# Combine this scaffold-bias with global-bias probability
					Pbias = {'dnaseinp':(Pbias['dnaseinp']+Pgbias['dnaseinp']), 'dnase':(Pbias['dnase']+Pgbias['dnase'])}
					# Rescale to probability
					Pt = float(sum(Pbias.values()))
					Pbias = {'dnaseinp':Pbias['dnaseinp']/Pt, 'dnase':Pbias['dnase']/Pt}
					# print 'global bias    : inp: %.3f, nuc: %.3f'%(Pgbias['inp'], Pgbias['nuc'])
					# print 'local bias     : inp: %.3f, nuc: %.3f'%(totalReads['inp'][chrom], totalReads['nuc'][chrom])
					pipeit('- %s cov. bias: dnaseinp: %.3f, nuc: %.3f'%(chrom, Pbias['dnaseinp'], Pbias['dnase']),1)
					pipeit('- %s cov. bias: dnaseinp: %.3f, nuc: %.3f'%(chrom, Pbias['dnaseinp'], Pbias['dnase']),1,outlogfh)
					
			except KeyError: 
				print >> sys.stderr, 'Warning: could not locate coverage information for %s from chromosome %s. Check file: %s. Defaulting to no bias.'%(sampleID, chrom, crdfile)
				Pbias = Pgbias
			
			except ZeroDivisionError:
				print >> sys.stderr, 'Warning: could not locate coverage information for %s from chromosome %s. Check file: %s. Defaulting to no bias.'%(sampleID, chrom, crdfile)
				Pbias = Pgbias
				
		else:
			# no reads for this scaffold; default to global bias
			
			# ATAC
			if nchannels == 1 and 'atac' in A2S and A2S['atac'] in SID:
				Pbias = {'atac':Pgbias['atac']}
				print '- %s default cov. bias: atac: %.3f'%(chrom, Pbias['atac'])
			elif nchannels == 2 and 'atac' in A2S and 'atacinp' in A2S and A2S['atacinp'] in SID and A2S['atac'] in SID:
				Pbias = {'atacinp':Pgbias['atacinp'], 'atac':Pgbias['atac']}
				print '- %s default cov. bias: atacinp: %.3f, atac: %.3f'%(chrom, Pbias['atacinp'], Pbias['atac'])
			
			# DNase
			elif nchannels == 1 and 'dnase' in A2S and A2S['dnase'] in SID:
				Pbias = {'dnase':Pgbias['dnase']}
				print '- %s default cov. bias: dnase: %.3f'%(chrom, Pbias['dnase'])
			elif nchannels == 2 and 'dnase' in A2S and 'dnaseinp' in A2S and A2S['dnaseinp'] in SID and A2S['dnase'] in SID:
				Pbias = {'dnaseinp':Pgbias['dnaseinp'], 'dnase':Pgbias['dnase']}
				print '- %s default cov. bias: dnaseinp: %.3f, dnase: %.3f'%(chrom, Pbias['dnaseinp'], Pbias['dnase'])
			
			# ChIP
			elif nchannels == 2 and 'inp' in A2S and 'nuc' in A2S and A2S['inp'] in SID and A2S['nuc'] in SID:
				Pbias = {'inp':Pgbias['inp'], 'nuc':Pgbias['nuc']}
				print '- %s default cov. bias: inp: %.3f, nuc: %.3f'%(chrom, Pbias['inp'], Pbias['nuc'])
			elif nchannels == 2 and 'ptm' in A2S and 'nuc' in A2S and  A2S['nuc'] in SID and A2S['ptm'] in SID:
				Pbias = {'ptm':Pgbias['ptm'], 'nuc':Pgbias['nuc']}
				print '- %s default cov. bias: nuc: %.3f, ptm: %.3f'%(chrom, Pbias['nuc'], Pbias['ptm'])
			elif nchannels == 2 and 'inp' in A2S and 'ptm' in A2S and  A2S['inp'] and A2S['ptm'] in SID:
				Pbias = {'inp':Pgbias['inp'], 'ptm':Pgbias['ptm']}
				print '- %s default cov. bias: inp: %.3f, ptm: %.3f'%(chrom, Pbias['inp'], Pbias['ptm'])
			elif nchannels == 3:
				Pbias = {'inp':Pgbias['inp'], 'nuc':Pgbias['nuc'], 'ptm':Pgbias['ptm']}
				print '- %s default cov. bias: inp: %.3f, nuc: %.3f, ptm: %.3f'%(chrom, Pbias['inp'], Pbias['nuc'], Pbias['ptm'])
			
		
		if DOBNORM:
			# Apply between sample normalization, and renormalize this to 1
			# =============================================================
			for aka in betWeight.keys(): Pbias[aka] += betWeight[aka]
			ttmp = float(sum(Pbias.values()))
			for aka in betWeight.keys(): Pbias[aka] = Pbias[aka]/ttmp
			pipeit('- Final local coverage correction: %s'%(' '.join(['%s=%.3f'%(x,y) for x,y in Pbias.items()])), 1)
			pipeit('- Final local coverage correction: %s'%(' '.join(['%s=%.3f'%(x,y) for x,y in Pbias.items()])), 1,outlogfh)
			
		# and print out a summary
		pbk = Pbias.keys(); pbk.sort()
		info = '\t'.join([chrom]+[str(Pbias[k]) for k in pbk])
		record(info, covbiasfile, 'a')
		
		# # Make a local copy in case we make scaffold-specific priors
		# Pr = dict(map(lambda x: (x,gPr[x]), gPr.keys())) # make a copy of the global prior values
		
		# scale prior by coverage bias... NO NO NO this doesn't do enough to control
		# Pr = dict(map(lambda x: (x, Pr[x] * Pbias[S2A[x]]), Pr.keys()))
		# Pt = float(sum(Pr.values())) # Rescale priors to sum to 1
		# Pr = dict( map(lambda k: (k, Pr[k]/Pt), Pr.keys()) )
		
		# coverage correction/bias proportions
		Kcov = [Pbias[S2A[anid]] for anid in SID]
		
		MDRunLength = 0; MDRunStart = 0
		for Gi in xrange(chromlen):
			totalLoci += 1     # count number of reference positions passed
			ref = None # reference genotype
			if GRNA[chrom]['seq'][Gi] == '*': 
				pipeit('ERROR| MISSING REFERENCE NT: %s'(GRNA[chrom]['seq'][Gi-5:Gi+5]), 1, outlogfh)
				continue # can't resequence without a reference
			
			frac = 100*(totalLoci/float(flg))
			if Gi % 100000 == 0: print '   %s|| %s %s/%s | Hypos=%s/%s | %.2f'%(keyid, chrom, Gi,chromlen, numKeyHypotheses, numHypotheses, frac)+'% complete'
			
			# get the log likelihood values for different binding outcomes
			lik = [[] for c in sampleIDs]
			info = ['-,-,-,-' for c in sampleIDs]
			depth = dict([(c,0) for c in sampleIDs])
			totalmultiplicity = [nastr for c in sampleIDs] # total number of alignments to other loci over all reads and channels
			dmean = [nastr for c in sampleIDs] # harmonic mean of total multiplicity per channel, ie avg other alignments
			MISSINGDATAFLAG = [0 for c in sampleIDs]
			
			# pipeit('NEXT POSITION: chrom=%s pos=%s (numHypos=%s)'%(chrom,Gi,numHypotheses), 1, outlogfh)
			
			# iterate over likelihood files
			for ci in rlc:
				c = sampleIDs[ci]
				Clik = None
				Cdep = 0
				tmul = 0
				Cd = 0
				Cref = 'N'
				
				# print 'chrom=%s, pos=%s, ci=%s, sample=%s'%(chrom, Gi, ci, c)
				# print '     chrm[c]=%s, pos[c]=%s'%(chrm[c], pos[c])
				
				# enclose the whole thing in try block to catch mal-formed lines
				try:
					# determine if this likelihood file contains info for this position
					if chrm[c] == chrom and pos[c] == Gi:
						# pipeit('  %s: match %s %s'%(c, chrom, Gi), 1, outlogfh)
						row = line[c][:-1].split('\t')
						Cchrom, Cpos, Cref, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:8]
						# print '- local', c, Cchrom, Cpos, '=', Gi
						Clik = float(ClikY)
						chrm[c] = Cchrom
						pos[c] = int(Cpos)
						info[ci] = '%s,%s,%s,%s'%(Cdep,Cmul,Ctmul,Cd)
						tmul = int(Ctmul)
						Cd = float(Cd)
						if Cchrom == chrom: matched[c] = True
						else: matched[c] = False
						# if Cchrom != chrom: break
					elif chrm[c] == chrom:
						# if 'scaffold' in chrom: pipeit('  %s: while %s'%(c, chrom), 1, outlogfh)
						while chrm[c] == chrom and pos[c] < Gi:
							try:
								line[c] = fh[c].next()
								row = line[c][:-1].split('\t')
								 # [chrom, str(Gi), ref, str(LL[0]), str(depth), str(multiplicity), str(sum(totalmultiplicity)), str(udalpha)]
								Cchrom, Cpos, Cref, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:8]
								# if 'scaffold' in chrom: pipeit('    -in while: %s =? %s'%(Cpos,Gi), 1, outlogfh)
								
								# Cchrom, Cpos, Cref, ClikN, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:9]
								Clik = float(ClikY)
								chrm[c] = Cchrom
								pos[c] = int(Cpos)
								info[ci] = '%s,%s,%s,%s'%(Cdep,Cmul,Ctmul,Cd)
								tmul = int(Ctmul)
								Cd = float(Cd)
							
								if Cchrom != chrom: matched[c] = False
								else: matched[c] = True
							except ValueError: 
								pipeit('Error parsing likelihood line %s: %s'%(c,line[c]),0)
								pipeit('ERROR| Parsing likelihood line %s: %s'%(c,line[c]),0,outlogfh)
								
							except StopIteration: break
							# print '     while: chrm[c]=%s, pos[c]=%s'%(chrm[c], pos[c])
					else:
						# different chroms
						# if 'scaffold' in chrom: pipeit('  %s: unmatched chrom, have %s want %s'%(c, chrm[c], chrom), 1, outlogfh)
						if pos[c] == -1 or matched[c]:
							try:
								line[c] = fh[c].next()
								row = line[c][:-1].split('\t')
								 # [chrom, str(Gi), ref, str(LL[0]), str(depth), str(multiplicity), str(sum(totalmultiplicity)), str(udalpha)]
								Cchrom, Cpos, Cref, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:8]
								# print '- local', c, Cchrom, Cpos, '=?', Gi
								
								# Cchrom, Cpos, Cref, ClikN, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:9]
								Clik = float(ClikY)
								chrm[c] = Cchrom
								pos[c] = int(Cpos)
								info[ci] = '%s,%s,%s,%s'%(Cdep,Cmul,Ctmul,Cd)
								tmul = int(Ctmul)
								Cd = float(Cd)
								if Cchrom != chrom: matched[c] = False
								else: matched[c] = True
							except StopIteration: pass
							# print '     chrdiff: chrm[c]=%s, pos[c]=%s'%(chrm[c], pos[c])
						elif chrm[c] < chrom: # if this chrm[c] is < chrom, while loop to catch up
							# print 'NEW while', c, chrm, pos, Gi
							# if 'scaffold' in chrom: pipeit('  %s: while, have %s want %s...'%(c, chrm[c],chrom), 1, outlogfh)
							# catchcount = 1
							while chrm[c] < chrom:
								try:
									line[c] = fh[c].next()
									row = line[c][:-1].split('\t')
									 # [chrom, str(Gi), ref, str(LL[0]), str(depth), str(multiplicity), str(sum(totalmultiplicity)), str(udalpha)]
									Cchrom, Cpos, Cref, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:8]
									# print '- local', c, Cchrom, Cpos, '=?', Gi
									
									# Cchrom, Cpos, Cref, ClikN, ClikY, Cdep, Cmul, Ctmul, Cd = row[0:9]
									Clik = float(ClikY)
									chrm[c] = Cchrom
									pos[c] = int(Cpos)
									info[ci] = '%s,%s,%s,%s'%(Cdep,Cmul,Ctmul,Cd)
									tmul = int(Ctmul)
									Cd = float(Cd)
									
									if Cchrom != chrom: matched[c] = False
									else: matched[c] = True
									# catchcount += 1
								except ValueError: 
									pipeit('Error parsing likelihood line %s: %s'%(c,line[c]),0)
									pipeit('ERROR| Parsing likelihood line %s: %s'%(c,line[c]),0,outlogfh)
								except StopIteration: break
								# if catchcount % 10000 == 0: print '     catchup: chrm[c]=%s, pos[c]=%s'%(chrm[c], pos[c])
				except ValueError: 
					# pipeit('Parsing error: passing...',1)
					pipeit('Error (outer) parsing likelihood line %s: %s'%(c,line[c]),0)
					pipeit('ERROR|OUTER Parsing likelihood line %s: %s'%(c,line[c]),0,outlogfh)
				
				ref = Cref
				
				try:
					# store likelihood etc for this position-sample
					if chrm[c] == chrom and pos[c] == Gi and Clik != None:
						lik[ci] = Clik
						depth[c] = int(Cdep)
						totalmultiplicity[ci] = tmul
						dmean[ci] = Cd
					else:
						# no data for this position-sample
						lik[ci] = pseudolog
						depth[c] = pseudocount
						MISSINGDATAFLAG[ci] = 1
				except ValueError:
					lik[ci] = pseudolog
					depth[c] = pseudocount
					MISSINGDATAFLAG[ci] = 1
					
			# ALL data
			totaldepth = sum(depth.values())
			
			# only permit missing data on the input lane
			if (nucidx!=None and MISSINGDATAFLAG[nucidx]) or (ptmidx!=None and MISSINGDATAFLAG[ptmidx]): 
				if MDRunLength == 0: MDRunStart = Gi
				MDRunLength += 1
				continue
			elif MDRunLength > 0: 
				errstr = 'ERROR| MISSING DATA: SKIPPING LOCUS %s %s-%s (%s) Hypos=%s %s missing: inpidx=%s, nucidx=%s, ptmidx=%s'%(chrom, MDRunStart,Gi, MDRunLength, numHypotheses, MISSINGDATAFLAG, inpidx,nucidx,ptmidx)
				pipeit(errstr, 1, outlogfh)
				MDRunLength = 0
			
			# ... now we have the likelihood values
			
			
			# Correct likelihoods for read length, diameter, and coverage (within and between) differences
			# ============================================================================================
			# algebra:
			# lik[x] = log of the product of read likelihoods for a sample of a particular read length (from glimmrLikelihood.py)
			# assuming constant likelihood, P(r_i|s) = Z**R, where Z is the probability and R is read length
			#        = log(Prod(P(r_i|s)))
			#        = Sum(log(P(r_i|s)))
			#        = Sum(log(Z**R)) = Sum(R log(Z)) = R Sum(log(Z))
			# so we can correct across samples by the ratio S/R, where S is the reference read length
			# which we will define as the smallest read length (probably 36nt)
			
			# Also includes correction for diameter difference, which is analogous
			#        = Sum_reads (log(Z**R)) = R Sum_reads(log(Z)),
			#          where number of reads is proportional to diameter x expected covergage
			#          but expected coverage is equivalent across samples, so this simplifies
			#          to a ratio of diam1/diam2. 
			# that said we can compute empirical coverage correction
			
			# Correct likelihoods for coverage differences
			# Note that here the effect is exponential on the likelihood, just as above
			# this is because the likelihood scales as a product of coverage on natural scale
			# but linear with coverage on the log scale
			
			lik = map(lambda ci: lik[ci]*Kll[ci]*Kcov[ci], rlc)
			
			
			# Compute posterior probability
			# ----------------------------------------------------------
			
			# Complete the likelihoods by missing transitive probabilities
			# multiply:
			# L(R|inp) = L(Rinp|inp) x L(rnuc|inp)^d_nuc x L(rptm|inp)^d_ptm
			# L(R|nuc) = L(Rnuc|nuc) x L(rinp|nuc)^d_inp x L(rptm|nuc)^d_ptm
			# L(R|ptm) = L(Rptm|ptm) x L(rinp|ptm)^d_inp x L(rnuc|ptm)^d_nuc
			
			if modeltype == 1: # nuc vs inp
				# given state is x, what is prob of x emitting data from condition y
				dn = depth[A2S['nuc']]; di = depth[A2S['inp']]
				# we normalized the likelihoods, so now we need to normalize the coverage itself
				di = Kcov[inpidx]*di
				dn = Kcov[nucidx]*dn
				
				transLik = {'inp': (lik[nucidx]*(1./dn)+lPerr)*dn,
				            'nuc': (lik[inpidx]*(1./di)+lPerr)*di}
				# unified likelihood values given each model condition
				lik = map(lambda ci: lik[ci] + transLik[S2A[SID[ci]]], rlc)
				
			elif modeltype == 2: # ptm vs nuc
				# given state is x, what is prob of x emitting data from condition y
				dn = depth[A2S['nuc']]; dp = depth[A2S['ptm']]
				# we normalized the likelihoods, so now we need to normalize the coverage itself
				dn = Kcov[nucidx]*dn
				dp = Kcov[ptmidx]*dp
				
				transLik = {'nuc': (lik[ptmidx]*(1./dp)+lPerr)*dp,
				            'ptm': (lik[nucidx]*(1./dn)+lZptm)*dn}
				            # 'ptm': (lik[nucidx]*(1./dn)+lPerr+lZptm)*dn}
				# unified likelihood values given each model condition
				lik = map(lambda ci: lik[ci] + transLik[S2A[SID[ci]]], rlc)
				
			elif modeltype == 3: # ptm vs inp
				# given state is x, what is prob of x emitting data from condition y
				di = depth[A2S['inp']]; dp = depth[A2S['ptm']]
				# we normalized the likelihoods, so now we need to normalize the coverage itself
				di = Kcov[inpidx]*di
				dp = Kcov[ptmidx]*dp
				
				transLik = {'inp': (lik[ptmidx]*(1./dp)+lPerr)*dp,
				            'ptm': (lik[inpidx]*(1./di)+lPerr)*di}
				# unified likelihood values given each model condition
				lik = map(lambda ci: lik[ci] + transLik[S2A[SID[ci]]], rlc)
				
			elif modeltype == 4: # ptm vs nuc vs inp
				# given state is x, what is prob of x emitting data from condition y
				dn = depth[A2S['nuc']]; di = depth[A2S['inp']]; dp = depth[A2S['ptm']]
				
				# we normalized the likelihoods, so now we need to normalize the coverage itself
				di = Kcov[inpidx]*di
				dn = Kcov[nucidx]*dn
				dp = Kcov[ptmidx]*dp
				
				# P(nuc read | ptm is true state) contains Zptm, which upweights the ptm likelihood based on increased expected TP
				# from the relative occurrence of ptm within nuc
				# that is, Zptm says that some fraction of the nuc reads are actually modified by the ptm
				try:
					transLik = {'inp': (lik[nucidx]*(1./dn)+lPerr)*dn + (lik[ptmidx]*(1./dp)+lPerr)*dp,
					            'nuc': (lik[inpidx]*(1./di)+lPerr)*di + (lik[ptmidx]*(1./dp)+lPerr)*dp,
					            'ptm': (lik[inpidx]*(1./di)+lPerr)*di + (lik[nucidx]*(1./dn)+lZptm)*dn}
					            # 'ptm': (lik[inpidx]*(1./di)+lPerr)*di + (lik[nucidx]*(1./dn)+lPerr+lZptm)*dn}
				
					# unified likelihood values given each model condition
					lik = map(lambda ci: lik[ci] + transLik[S2A[SID[ci]]], rlc)
				except ZeroDivisionError: continue
			
			elif modeltype == 5: # atac vs baseline
				# take a real file and scramble the chrom/pos locations for atac control?
				sys.exit("modeltype 5 not implemented")
				
			elif modeltype == 6: # atac vs atacinp
				# given state is x, what is prob of x emitting data from condition y
				dn = depth[A2S['atac']]; di = depth[A2S['atacinp']]
				# we normalized the likelihoods, so now we need to normalize the coverage itself
				di = Kcov[inpidx]*di
				dn = Kcov[nucidx]*dn
				
				transLik = {'atacinp': (lik[nucidx]*(1./dn)+lPerr)*dn,
				            'atac': (lik[inpidx]*(1./di)+lPerr)*di}
				# unified likelihood values given each model condition
				lik = map(lambda ci: lik[ci] + transLik[S2A[SID[ci]]], rlc)

			elif modeltype == 7: # dnase vs dnaseinp
				# given state is x, what is prob of x emitting data from condition y
				dn = depth[A2S['dnase']]; di = depth[A2S['dnaseinp']]
				# we normalized the likelihoods, so now we need to normalize the coverage itself
				di = Kcov[inpidx]*di
				dn = Kcov[nucidx]*dn
				
				transLik = {'dnaseinp': (lik[nucidx]*(1./dn)+lPerr)*dn,
				            'dnase': (lik[inpidx]*(1./di)+lPerr)*di}
				# unified likelihood values given each model condition
				lik = map(lambda ci: lik[ci] + transLik[S2A[SID[ci]]], rlc)
				
			
			# Scale down the likelihoods so that natural scale values don't become 0
			# ------------------------------------------------------------------------------
			# we are currently on log scale so values are large negative exponents
			if (nchannels == 2 and lik[0]+lik[1] < -500) or (nchannels == 3 and lik[0]+lik[1]+lik[2] < -750):
				# to avoid numerical underflow due to exponentiation of large negative LL values
				# need to find a positive constant C such that math.e^LL+C is computable, eg < 1000
				
				# print >> sys.stderr, 'Warning: small log likelihoods:', lik
				
				# precursor step is to do one blind rescaling to bring values near computable range
				# if you don't do this, you will inevitably set many positions values to equality
				# whereas they could be quantified just fine
				maxll = max(lik)
				while maxll < -200:
					c = -maxll/4. # this can't become positive
					lik = map(lambda x: x+c, lik)
					# print >> sys.stderr, 'Step:', lik
					maxll = max(lik)
				
				badvalue = -100000000000000000000000000000000 # token indicator of too small value
				lik = map(lambda x: (x < -250 and badvalue) or x, lik)
				
				# confirm that not all have badvalue
				if len(filter(lambda x: x!=badvalue,lik)) == 0:
					pipeit('Warning: All likelihoods are 0 at position %s'%(Gi),1)
					pipeit('Warning: All likelihoods are 0 at position %s'%(Gi),1,outlogfh)
					# then all likelihoods have the same value
					lik = map(lambda x: -250, lik)
				else:
					idx = argmax(lik) # get the largest likelihood
					candC = -lik[idx] - 100 # makes the max possible log likelihood = -100: 10**-100 = 1e-100
					# all other liklihoods will be smaller than this
					lik = map(lambda x: (lik[x] == badvalue and -250) or (lik[x]+candC), rlc)
					
					
					# try: 
					# 	lik = map(lambda x: (lik[x] == badvalue and -250) or (lik[x]+candC), rlc)
					# 	# lik = map(lambda x: (lik[x] == badvalue and 1e-500) or 10**(lik[x]+candC), rlc)
					# except OverflowError:
					# 	# this should not except because we are exponentiating negative numbers >= 100
					# 	print '  ERROR:', 'C', candC, 'LL[%s]'%(idx), lik[idx], '=>', lik[idx]+candC#, '=>', math.e**(LL[idx]+C[idx])
					# 	print '  CANNOT EXPONENTIATE'
					# 	print '  LL', lik
					# 	# resolve by ignoring this genomic position
					# 	continue
					# 	
				# print >> sys.stderr, 'Corrected:', lik
			# else:
			# in this case likelihoods are fine, so just restore natural scale
			try: lik = map(lambda x: 10**x, lik)
			except OverflowError: 
				pipeit('ERROR| Stack overflow in exponentiating log likelihood at position %s: %s'%(Gi,lik),1,outlogfh)
				continue
			# DONE SMALL LIKELIHOOD CORRECTION
			
			# we are combining 2 or 3 harmonic means - use the standard average
			# udalpha = round(avg0(dmean), 3)
			
			try: dmean = map(floatna, dmean)
			except ValueError: continue
			
			# degeneracy with respect to key sample only - notably input degeneracy is likely higher
			udalpha = dmean[keyidx]
			
			try:
				ud = map(lambda x: x==nastr and 1 or 1+x, dmean)
			except TypeError:
				sys.stderr.write("Error: unexpected dmean: %s"%(dmean))
				pipeit("ERROR| Unexpected dmean: %s"%(dmean), 1, outlogfh)
				pipeit('ERROR| CRAZY HIGH LIKLIHOODS: SKIPPING LOCUS %s %s %s'%(chrom, Gi, lik), 1, outlogfh)
				
				continue
			
			# Posterior probabilities
			# ----------------------------------------------
			# we now have a length-10 vector of likelihood values 
			# for all reads at Gi given each genotype
			# calculate the posterior probability of each possible genotype
			# given the read set
			
			# we must calculate the PGxyR for all possible genomic string configurations C(SGxy(l))
			# assuming genotype at position l is xy
			# so sum up prob over the C(SGxy(l)) for each genotype xy
			
			# so basic assumption is G = GR so that C(SGxy(l)) is simply the k-bounded
			# set of the reference: C(SGxy(l)) = C(SGRxy(l))
			# this way there is no summing involved
			
			# Rescale the null likelihood based on the overall degeneracy of the locus
			# Essentialy we assume that if a protein/mark is present at this locus
			# then it cannot be present at the other homologous loci in the genome...
			# we compensate somewhat by capping the correction at 1/prior
			# thus if the number of paralogous loci exceeds the expectation for only 1 binding site
			# then we do not correct further, thus potentially allowing another binding event
			
			# this is a little conservative for total histone - we could adjust 
			altweight = udalpha*Pr[keyid]
			nulweight = max(1e-10, 1.-altweight) # if highly degenerate then set prior to minimum of 1e-10
			numer = map(lambda ci: ci==inpidx and (Pr[keyid]<0.33 or udalpha>3) and lik[ci]*nulweight + lik[keyidx]*altweight or lik[ci]*Pr[SID[ci]], rlc)
			# numer = map(lambda ci: ci==inpidx and lik[ci]*nulweight + lik[keyidx]*altweight or lik[ci]*Pr[SID[ci]], rlc)
			
			# standard computation
			# numer = map(lambda ci: lik[ci]*Pr[SID[ci]], rlc)
			
			# convert to decimal object for greater precision of PP value
			# numer = map(lambda x: decimal.Decimal(repr(x)), numer)
			
			denom = sum(numer)
			
			# loss of precision can occur here
			try: PP = map(lambda x: x/denom, numer)
			except ZeroDivisionError: 
				print >> sys.stderr, 'total numer mass is 0: L=%s, tL=%s'%(lik, transLik)
				continue
			
			# Determine MAP genotype and posterior
			# ------------------------------------
			# take argmax, but multiple genotypes may have similar or identical PP. 
			# This would happen if the pileup is highly degenerate. 
			idx = argmax(PP)
			model = sampleIDs[idx]
			
			# Record information for loci having desired model state
			if model == keyid or SAVEALLPOSTERIOR:
				
				# first try hardware/float precision
				MAP = PP[idx]
				try: 
					Q = -10*log(1-MAP)
				except ValueError:
					# recalculate with software precision
					# numer2 = map(lambda x: decimal.Decimal(repr(x)), numer)
					# MAP = numer2[idx]/sum(numer2)
					# Q = -10*(1-MAP).log10() # convert to log scale
					
					# must introduce minor truncation to maintain good runtime
					Q = infy; cprec = defaultprec
					try:
						while Q == infy and cprec <= 210:
							decimal.getcontext().prec = cprec
							numer2 = map(lambda x: decimal.Decimal(repr(x)), numer)
							MAP = numer2[idx]/sum(numer2)
							Q = -10*(1-MAP).log10()
							cprec = int(cprec * 1.3)
							# print Gi, 'trial', decimal.getcontext().prec, Q
					
						# if cprec > 190: print 'end of day', Gi, int(cprec/1.3), Q
					except Exception,e: continue
					
					if Q == infy: Q = lowestQ
					
					# caution should reset decimal precision if using again below
					# decimal.getcontext().prec = defaultprec # reset
				
				
				# carry over all of the degeneracy information from each channel
				likstr = map(lambda x: lik[x]==nastr and nastr or '%.3e'%(lik[x]), rlc)
				dstr = map(lambda x: '%s:%s'%(S2A[sampleIDs[x]],'%s,'%(likstr[x])+info[x]), rlc)
				rat1 = nastr; rat2 = nastr
				
				# ChIP-seq models
				if modeltype == 1:
					rat1 = 1
					rat2 = roundna(3)(divna(  divna(depth[nucsample],ud[nucidx]), divna(depth[inpsample],ud[inpidx])))
				elif modeltype == 2:
					rat1 = roundna(3)(divna( divna(depth[ptmsample],ud[ptmidx]), divna(depth[nucsample],ud[nucidx])))
					rat2 = 1
				elif modeltype == 3:
					rat1 = 1
					rat2 = roundna(3)(divna( divna(depth[ptmsample],ud[ptmidx]), divna(depth[inpsample],ud[inpidx])))
				elif modeltype == 4: 
					rat1 = roundna(3)(divna( divna(depth[ptmsample],ud[ptmidx]), divna(depth[nucsample],ud[nucidx])))
					rat2 = roundna(3)(divna( divna(depth[ptmsample],ud[ptmidx]), divna(depth[inpsample],ud[inpidx])))
				
				# ATAC-seq models
				elif modeltype == 6 or modeltype == 7:
					rat1 = 1
					rat2 = roundna(3)(divna(  divna(depth[nucsample],ud[nucidx]), divna(depth[inpsample],ud[inpidx])))
				
				
				info2 = [chrom, str(Gi), ref, str(S2A[model]), '%.1f'%(Q), str(totaldepth), '%.2f'%(rat1), '%.2f'%(rat2)]+dstr
				
				# if you want the total likelihood
				# info2 = [chrom, str(Gi), ref, str(S2A[model]), '%.1f'%(Q), str(totaldepth), '%.3f'%(rat1), '%.3f'%(rat2), str(-log(denom))]+dstr
				print >> outfh, '\t'.join(info2)
				numKeyHypotheses += 1
				
			numHypotheses += 1 # count up number of hypotheses we make
		
		if MDRunLength > 0:
			errstr = 'ERROR| MISSING DATA: SKIPPING LOCUS %s %s-%s (%s) Hypos=%s %s missing: inpidx=%s, nucidx=%s, ptmidx=%s'%(chrom, MDRunStart,Gi, MDRunLength, numHypotheses, MISSINGDATAFLAG, inpidx,nucidx,ptmidx)
			pipeit(errstr, 1, outlogfh)
		
	outfh.close()
	outlogfh.close()
else: print 'File %s found. DONE'%(postfile)
# DONE


# now filter out the significant calls for our ptmsample
# correct for multiple testing if this hasn't been done yet

# first check if we need to recompute anything
checkups = []
for minPvalue in PvalueRange:
	# multiple-test corrected snp variants
	keyfile = theptmdir+'%s.posterior'%(keyid)+'.Q%s%s.txt'%(minPvalue,boundslabel)
	if os.access(keyfile, os.F_OK) and not REDOFINALPEAKCALLING: continue
	# checkups += [(minPvalue,keyfile)]
	checkups += [minPvalue]
	
checkups.sort()
qccheck = lambda count: lambda line: len(line[:-1].split('\t'))>=count and True

if len(checkups)>0:
	pipeit('- Identifying significant loci from the %s/%s hypotheses...'%(numHypotheses, totalLoci),1)
	
	# first get total number of hypotheses from mother
	# and store data
	
	print '- Counting number of genomic positions evaluated from likelihood files...'
	
	# these files will be huge, so really need to do a wc
	# if not os.access(tmpdir+'motherwc.txt', os.F_OK):
	
	numHypos = {}
	for sid in sampleIDs:
		call = 'wc -l "%s" > "%s"'%(motherfiles[sid], tmpdir+sid+'.motherwc.txt')
		os.system(call)
		print call
		tmptmp = readList(tmpdir+sid+'.motherwc.txt')
		rem = tmptmp[0].split()
		if type(rem) == type([]): numHypos[sid] = int(rem[0])
		else: numHypos[sid] = int(rem)
		print '   %s: number of hypotheses:'%(sid), numHypos[sid]
	
	mnh = min(numHypos.values())
	mxh = max(numHypos.values())
	
	# if mother file contained no information
	if numHypos == 0: numHypos = 1
	print '- %s < final num hypotheses < %s...'%(mnh, mxh)
	
	nhstr = ', '.join(map(lambda x: '%s:%s'%(x,numHypos[x]), sampleIDs))
	
	numHypos = mxh # just use the max channel to report percentages
	
	# now we can post-filter the significant set 
	# based on the actual number of hypotheses we made
	for Pvalue in checkups:
		
		# generate list of significant genotypes
		keyfile = theptmdir+'%s.posterior.Q%s%s.txt'%(keyid,Pvalue,boundslabel)
		
		nsig = 0
		
		thefiletouse = theptmdir+'%s.posterior%s.txt'%(keyid,boundslabel)
		# print 'test', thefiletouse, '|', boundslabel
		if Pvalue != PvalueRange[0]:
			idx = PvalueRange.index(Pvalue)
			thefiletouse = theptmdir+'%s.posterior.Q%s%s.txt'%(keyid,PvalueRange[idx-1],boundslabel)
		infh = open(thefiletouse)
		
		outfh = open(keyfile, 'w')
		print >>outfh, '# Total positions tested: %s'%(nhstr)
		
		# for i in range(len(data)):
		for line in infh:
			if not qccheck(10)(line): continue
			row = line[:-1].split('\t')
			
			Q = (row[4]==nastr or row[4]=='NA' and nastr) or float(row[4])
			if Q != nastr and Q >= Pvalue:
				info = '\t'.join(row)
				print >>outfh, info
				nsig += 1
				
		infh.close()
		outfh.close()
		print '- There are %s/%s bound loci at Q > %s (P < %s).'%(nsig, numHypos, Pvalue, 10**(-Pvalue/10.))
	infh.close()
		
		# filter this file into components?
# end loop

# CLEANUP
os.system('rm -r "%s"'%(tmpdir)) # remove tmpdir

# finished processing, so write to an outfile so parent process knows to launch another
PID = ','.join(map(str,singlesampleid))
fin = open(completedir+'posterior_%s.done'%(PID), 'w'); fin.close()

print '\nResults saved to %s'%(theptmdir)
print 'Exiting glimmrHet.py...\n'

