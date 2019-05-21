#! /usr/bin/env python

import os, re, sys, math, time, stat, time, random, subprocess, decimal, glob
import glimmrAccessories as ga

VERSION = "1.3"

help = """Usage: %s [options]
version %s

Description: Run Kallisto on a directory of fastq files

Author: Daniel F. Simola
Email: dfsimola@fastmail.fm
Created: 16 May 2018


Version history:
v1.3  - added CLI parms for fragment length and sd
v1.25 - various small bug fixes
v1.2  - removed --kalisto-format flag; by default will save TPM, count, and FPKM matrices
v1.1: - changed default from outer join to simple concatenation which assumes input samples were processed using same index (--outerjoin to enable)
      - Added ability to source fastq file paths from metadata file under headers specified by user (--mapfastq R1,R2)
      - Added additional parallelization to specify threads per Kallisto call (-x or -xthreads) as well as parallel calls (-xjobs)

"""%(sys.argv[0], VERSION)

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
outdir = ''
REDO = False

ensgtfile = 'ENSEMBL.transcript.annotation.GRCh38.p10.txt' 
KAL_INDEX = 'GRCh38.kallisto.index'
KAL_PATH = 'kallisto_linux-v0.44.0/kallisto'

FLIPG2T = False

KAL_BOOTSTRAPS = 0
KAL_SINGLEEND = '--single'
KAL_FRAGLEN = 250
KAL_FRAGSD = 50
KAL_GCBIAS = False
KAL_CUSTOM = "" # custom string fed to kallisto call
nprocs = 1    # number of processors per single Kallisto execution
nthreads = 1  # number of parallel calls to Kallisto (with nprocs threads per job)

PE1token = '_R1'
PE2token = '_R2'

fname = 'kallisto' # user-settable prefix for resulting tables
# KAL_format = 'TPM' # or count or FPKM
KAL_length = 'eff_length' # for FPKM conversion

mapfile = '' # optional fastq to SampleID labeling of samples
mapkey = 'SampleID' # unique key to merge against column names of expression matrix
filenamekey = 'Filename' # basename of a fastq file, e.g. ABCD_DS123_ATTACTCG_L001_R1.fastq.gz has Filename of ABCD_DS123_ATTACTCG_L001
FASTQFROMMETADATA = None# = ['trimmedFastqR1','trimmedFastqR2']


ONLY_MERGE = False
OUTER_JOIN = False
TEST_FOR_ACCESS = False

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
			
	for job in running_jobs: job.wait() # halt code until all jobs are finished.
	return 1



# CLI ARGS
# -------------------------------------------------------------------
help = """
Usage: %s --in/-i <FASTQ_DIRECTORY/> --index <KALLISTO_INDEX_FILE> [--out/-o <OUTDIR> -xjobs parallel_jobs<INT> -x num_threads<INT> --map <FILENAME>]

Additional flags:
--label <STRING>                      # prefix of experiment for saving resulting files
--redo                                # Redo Kallisto analysis for previously processed files 
--kallisto-format <TPM,count,FPKM>    # How to report expression values
--kallisto-length <length|eff_length> # only needed if --kallisto-format FPKM
--kallisto-path </path/to/kallisto>   # path to executable
--samplekey <STRING=SampleID>         # Which field to parse SampleID (if not "SampleID")
--filekey <STRING=Filename>           # Which field to parse Filename (if not "Filename")
--mapfastq <STRING,STRING>            # Column name (for SE) or names (for PE) denoting hard paths to fastq_R1 [and fastq_R2] files, 
                                      # e.g. --mapfastq fastqR1,fastqR2
                                      # assumes mapfile provided via --map as source of fastq files
--gcbias                              # GC bias flag
--fraglen                             # fragment length mean (default: 250)
--fragsd                              # fragment length SD (default: 50)
--petokens <_R1,_R2>                  # string pattern to match PE fastq files (default '_R1,_R2')
--annotation <path/to/gene.ann/       # path to gene to name annotation file
--custom <STRING>                     # any custom flags provided to kallisto
--merge                               # skip alignment+quant, even if there are files that have not been processed
--outer_join                          # perform full outer join (needed if combining samples from different reference indices)
--flipannotation

Note, you can specify a subset of fastq files in a directory using wild card, e.g. -i "/path/to/dir/HYF3M_DS59787A_TCCGCGAA-AGGCTATA_L001_R*"

Mapfile is of the format:
SampleID\tFilename
sampleA_1\tHYF3M_DS59787A_TCCGCGAA-AGGCTATA_L001 (filename minus suffix and removing _R1/_R2)
...

For example, if real filename was /path/to/file_L001_R1.fastq.gz, enter simply file_L001 in the Filename column of the metadata file.

"""%(sys.argv[0])
nhelps = 0; helplimit = 0
args = sys.argv[:]
argstr = ''.join(args)
ai = 1
userformat = False

while ai < len(args):
	arg = args[ai].strip('-').strip('--')#.lower()
	try: val = args[ai+1]
	except IndexError: val = ''
	
	if re.match('^in$|^i$', arg): indirfile = val; nhelps += 1
	elif re.match('^out$|^o$', arg): outdir = slash(val)
	elif arg == 'label': fname = val
	
	elif arg in ['index', 'kalidx','kallisto-index']: KAL_INDEX = val
	elif arg == 'kallisto-path': KAL_PATH = val
	# elif arg == 'kallisto-format': KAL_format = val # 'TPM' or 'FPKM' or 'count'
	elif arg == 'kallisto-length': KAL_length = val # 'length' or 'eff_length'
	elif arg == 'gcbias': KAL_GCBIAS = True; ai-=1
	elif arg == 'fraglen': KAL_FRAGLEN = int(val)
	elif arg == 'fragsd': KAL_FRAGSD = int(val)
	
	elif arg == 'testaccess': TEST_FOR_ACCESS = True; ai-=1
	
	elif arg == 'annotation': ensgtfile = val
	elif arg == 'flipannotation'; FLIPG2T = True; ai-=1

	elif arg == 'map' or arg == 'metadata': mapfile = val
	elif arg == 'mapkey' or arg == 'key' or arg == 'samplekey': mapkey = val
	elif arg == 'filekey': filenamekey = val
	elif arg == 'mapfastq': FASTQFROMMETADATA = val.split(',')
	
	elif arg == 'PEtokens': PE1token,PE2token = val.split(',')
	
	elif arg == 'custom': KAL_CUSTOM = " %s "%(val)
	elif arg == 'redo': REDO = True; ai-=1
	elif arg == 'merge': ONLY_MERGE = True; ai-=1
	elif arg == 'outerjoin': outer_join = True; ai-=1
	elif arg == 'xjobs': nprocs = int(val)
	elif arg == 'xthreads' or arg == 'x': nthreads = int(val)

	else:
		help += "=> The specified argument \""+arg+"\" does not parse."
		print >> sys.stderr, help
		sys.exit()
	ai += 2
if nhelps < helplimit or re.match('.*-h ', argstr): 
	sys.exit(help)

ga.createdir(outdir)

ga.printList(['python ' + ' '.join(args)], file=outdir+'logfile.txt', ioMethod='a')

# call = '- Quantifying TPM using kallisto for files in directory: %s'%(indirfile)
if not outdir: 
	outdir = 'kallisto/'
	ga.createdir(outdir)

basecall = '%s quant -i "%s" -b %s -t %s'%(KAL_PATH, KAL_INDEX, KAL_BOOTSTRAPS, nthreads)
basecall += KAL_CUSTOM
if KAL_GCBIAS: basecall += ' --bias'

SEstr = '%s -l %s -s %s'%(KAL_SINGLEEND,KAL_FRAGLEN,KAL_FRAGSD)

# iterate over unique file names
fastqfiles = []
sampleGroups = {}

if FASTQFROMMETADATA:
	if not mapfile or not os.access(mapfile, os.F_OK):
		sys.exit('-mapfastq specified but cannot access mapfile "%s"'%(mapfile), 1)

	pipeit('- Loading mapfile: %s'%(mapfile),1)
	md,mr,mc = ga.readTable(mapfile,rownames=0)
	
	labels = ''
	try:
		samplekeycol = mc.index(mapkey)
		try:
			labels = ga.vslice(md,samplekeycol)
		except IndexError:
			pipeit('ERROR: cannot parse data from column %s for samplekey "%s".'%(samplekeycol, mapkey),1)
			ga.printFormattedTable(md[:5], mc, pipe=sys.stderr)
			sys.exit('\n\nEXITING.')
	except IndexError:
		sys.exit('ERROR: cannot parse samplekey "%s" from map header: %s.\n\nEXITING.'%(mapkey, mc))
		

	if len(FASTQFROMMETADATA) == 1:
		fkey = mc.index(FASTQFROMMETADATA[0])
		fastqfiles = ga.vslice(md,fkey)
		sampleGroups = dict(zip(labels,fastqfiles))
	
	elif len(FASTQFROMMETADATA) == 2:
		f1key = mc.index(FASTQFROMMETADATA[0])
		f2key = mc.index(FASTQFROMMETADATA[1])
		fpairs = zip(ga.vslice(md,f1key), ga.vslice(md,f2key))
		sampleGroups = dict(zip(labels,fpairs))
		fastqfiles = ga.flatten(fpairs)

	else:
		sys.exit('ERROR PARSING FASTQFROMMETADATA: %s'%(FASTQFROMMETADATA), 1)

else:
	if '*' in indirfile: fastqfiles = glob.glob(indirfile)
	else: fastqfiles = ga.getFiles(indirfile, include='fastq')


call = '- Output directory: %s'%(outdir)
pipeit(call,1)
pipeit('- Found %s files in input directory.'%(len(fastqfiles)),2)
pipeit('  First 10:', 1)
ga.printList(fastqfiles[:10])
if not len(fastqfiles): sys.exit('NO DATA. EXITING.')
pipeit('',1)


if TEST_FOR_ACCESS:
	yes = 0
	no = 0
	pipeit('- Confirming access to fastq files...',1)

	for f in fastqfiles:
		yesno = os.access(f, os.F_OK)
		if yesno == True: yes += 1
		else: no += 1
		pipeit('- ACCESS? %s => %s'%(f, yesno), 1)
	pipeit('\n\nCan access %s/%s'%(yes,len(fastqfiles)),1)
	pipeit('CANNOT access %s/%s'%(no,len(fastqfiles)),1)
	sys.exit('DONE')


# group into paired end
if not sampleGroups:
	for f in fastqfiles:
		lab = ga.getLabel(f).split('.fastq')[0]
		if PE1token in lab: lab = lab.split(PE1token)[0]
		if PE2token in lab: lab = lab.split(PE2token)[0]

		try: sampleGroups[lab] += [f]
		except KeyError: sampleGroups[lab] = [f]

count = 1
jobs = []
fhco = open(outdir+'call.log.txt', 'w')
for lab in sorted(sampleGroups.keys()):
	labfiles = sampleGroups[lab]
	
	call = basecall + ' -o "%s"'%(outdir+'%s/'%(lab))
	
	if len(labfiles) == 1:
		call += ' %s "%s"'%(SEstr, labfiles[0])
	elif len(labfiles) == 2:
		call += ' "%s" "%s"'%(labfiles[0],labfiles[1])
	else:
		pipeit('- kalliso error: more than 2 files for unique label %s: %s'%(lab, labfiles),1)
	
	pipeit('- Call %s: %s'%(count,call), nl=1, pipe=fhco)
	
	if not ONLY_MERGE and not os.access(outdir+'%s/abundance.tsv'%(lab), os.F_OK) and not REDO:
		pipeit('- Call %s: %s'%(count,call),1)
		jobs += [call]
		# os.system(call)
	else:
		pipeit('- Already completed %s.'%(lab),1)
	count += 1
fhco.close()

pipeit('',1)

# launch jobs
launchProcesses(jobs, nsubprocs=nprocs, sleep=10)


pipeit('\n- Merging data...',1)

if mapfile: 
	pipeit('- Merge kallisto files from mapfile: %s'%(mapfile), 1)
pipeit('- Input directory: %s'%(outdir), 1)

# if not fname: fname = 'kallisto'

enst2g = {}   # transcript to gene id
# enst2len = {} # transcript length for FPKM normalization

try:
	pipeit('- Loading transcript to gene annotation file...')
	d,r,c = ga.readTable(ensgtfile,rownames=0)
	for row in d:
		# print 'loading %s - %s'%(row[0], row[1])
		if FLIPG2T:
			enst2g[row[0]] = row[1]
		else:
			enst2g[row[1]] = row[0]
		
		# enst2len[row[1]] = int(row[2])
	pipeit('done.',1)
	pipeit('Note, expected format is geneid\ttranscriptid. Use -flipg2t to reverse',1)
except IOError:
	pipeit('WARNING: Cannot access -annotation file, i.e. the gene-transcript annotation file.', 1)

kaldirs = ga.getDirectories(outdir)

samples = map(ga.getDirectory,kaldirs)
# link file name to sampleID
N2S = dict(zip(samples,samples))

if mapfile and os.access(mapfile, os.F_OK) and not FASTQFROMMETADATA:
	N2S = {}
	pipeit('- Loading mapfile: %s'%(mapfile),1)
	md,mr,mc = ga.readTable(mapfile,rownames=0)
	sididx = mc.index(mapkey)
	fnidx = mc.index(filenamekey)
	
	for row in md:
		fnfull = row[fnidx]
		# parse just the file part
		lab = fnfull.split('/')[-1]
		# split if PE
		lab = lab.replace('R[12]', 'R1')
		lab = lab.split('.fastq')[0]
		if PE1token in lab: lab = lab.split(PE1token)[0]
		if PE2token in lab: lab = lab.split(PE2token)[0]
		# print 'hey', lab, samples
		match = filter(lambda x: lab in x, samples)
		for l in match: N2S[l] = row[sididx]

elif mapfile and not os.access(mapfile, os.F_OK):
	print 'Cannot access mapfile "%s"'%(mapfile)
else:
	# no mapfile - just use directory name as sampleID
	for sid in samples:
		N2S[sid] = sid

print '\nSample to Label mapping (10 entries):'
ga.printList(N2S.items()[:10])
print

# merge FPKM columns into single matrix
labels = []
tab_df = []

pairs = {}
pairs['count'] = []
pairs['TPM'] = []
pairs['FPKM'] = []

# CONTINUE HERE!
# mat = [[ga.nastr for l in labels] for g in genes]

genetab = []

# alignment stats - parsing json file
infodct = {}
countidx = 1
for tmpdir in kaldirs:
	try:
		d,r,c = ga.readTable(tmpdir+'abundance.tsv',rownames=0)
		pipeit('- working %s: %s'%(countidx, tmpdir+'abundance.tsv'),1)
		lenidx = c.index(KAL_length)

		malab = ga.getDirectory(tmpdir)
		labels += [malab]
		
		json = ga.readList(tmpdir+'run_info.json')[1:-1]
		json = map(lambda x: "".join(x[2:]).split('":'), json)
		json = dict(map(lambda y: [y[0], y[1].strip(' ').strip(',')], json))
		# print 'TEST', json
		infodct[ga.getDirectory(tmpdir)] = json
		
		
		idx = c.index('tpm')
		pairs['TPM'] += [[(row[0], row[idx]) for row in d]]
		pairs['FPKM'] = [(row[0], ga.divna(ga.multna(1000,row[idx]),row[lenidx])) for row in d]
		
		idx = c.index('est_counts')
		pairs['count'] += [[(row[0], row[idx]) for row in d]]
		
		tab_df += [[malab, N2S[malab]]+[row[0].split('|')[0] in enst2g and enst2g[row[0].split('|')[0]] or "NA", row[0].split('|')[0], row[0], row[c.index('tpm')], row[c.index('est_counts')]] for row in d]
		
		# if KAL_format == 'TPM':
		# 	idx = c.index('tpm')
		# 	pairs += [[(row[0], row[idx]) for row in d]]
		# 	# tab += [[malab, N2S[malab]]+[row[0].split('.')[0] in enst2g and enst2g[row[0].split('.')[0] or "NA"],row[0],row[idx]] for row in d]
		# elif KAL_format == 'count':
		# 	idx = c.index('est_counts')
		# 	pairs += [[(row[0], row[idx]) for row in d]]
		# 	# tab += [[malab, N2S[malab]]+[row[0].split('.')[0] in enst2g and enst2g[row[0].split('.')[0] or "NA"],row[0],row[idx]] for row in d]
		# elif KAL_format == 'FPKM':
		# 	idx = c.index('tpm')
		# 	fpkm = [(row[0], ga.divna(ga.multna(1000,row[idx]),row[lenidx])) for row in d]
		# 	pairs += [fpkm]
		# 	# Note conversion from TPM to FPKM here!
		# 	# tab += [ [malab, N2S[malab]]+[row[0].split('.')[0] in enst2g and enst2g[row[0].split('.')[0]] or "NA", row[0], ga.divna(ga.multna(1000,row[idx]),row[lenidx])] for row in d]
		# else:
		# 	sys.exit('Kallisto format must be one of TPM or FPKM.')
	except IOError: pass
	countidx += 1

samples = sorted(infodct.keys())
tab2 = []
header = sorted(ga.unique(ga.flatten([infodct[sid].keys() for sid in samples])))
for sid in samples:
	row = [N2S[sid], sid] + [infodct[sid][k] for k in header]
	tab2 += [row]
ga.printTable(tab2,['SampleID', 'Filename']+header,file=outdir+'run_info.table.txt')
pipeit('- Saved run info to %s'%(outdir+'run_info.table.txt'),1)


for KAL_format in ['TPM', 'count']:#, 'FPKM']:
	matname = outdir+'%s.%s.%s.matrix.txt'%(fname,'transcript',KAL_format)

	if OUTER_JOIN:
		# outer join - fill with NA
		pipeit('- Performing outer join to merge data for %s (may take a few minutes)...'%(KAL_format))
	
		names, mat = ga.joinSet(pairs[KAL_format],inner=0)
		# get rid of ancillary header
		names2 = map(lambda x: x.split('|')[0], names)
		ga.printTable(mat,['Transcript']+labels, names2, file=matname)
		
	else:
		# default, which is to simply concatenate the columns from each file assuming row order is invariant
		matT = []
		names = []
		# just load each file, transpose, and cat onto tab
		# each entry in pair is vector of names and values
		for pair in pairs[KAL_format]:
			Xvec,Yvec = ga.transpose(pair)
			matT += [Yvec]
			if not len(names): names = Xvec
			elif Xvec != names: sys.exit('Error concatenating results. Order of gene names incorrect!')

		mat = ga.transpose(matT)
		# get rid of ancillary header
		names2 = map(lambda x: x.split('|')[0], names)
		ga.printTable(mat,['Transcript']+labels, names2, file=matname)

	pipeit('\n- Saved %s'%(matname),1)
	
	# print 'names2', names2[:10]

	# create gene model matrix as the sum of the transcript values per gene per sample
	genematd = {} # keyed by gene
	nvalid = {}
	nerrors = {}
	for i in range(len(names2)):
		enst = names2[i]
		ensg = ga.nastr
		# print('HEY %s'%(enst))
		try:
			# print("TEST %s %s"%(enst, enst.split('.')[0] in enst2g))
			# print('enst2g: %s'%(enst2g[enst2g.keys()[0]]))
			ensg = enst2g[enst]#.split('.')[0]
			try: 
				genematd[ensg] += [mat[i]]
				nvalid[enst] = 1
			except KeyError: 
				genematd[ensg] = [mat[i]]
		except KeyError:
			# may be ENST point version issue. trim and try again
			try:
				enst2 = enst.split('.')[0]
				# print "bollocks. try again with %s"%(enst2)
				ensg = enst2g[enst2]#.split('.')[0]
				try:
					genematd[ensg] += [mat[i]]
					nvalid[enst] = 1
				except KeyError: 
					genematd[ensg] = [mat[i]]
				# print('version error MATCH %s/%s'%(enst2, ensg))
			except KeyError: 
				# print('BAD BAD BAD enst2g: %s'%(enst2g[enst2g.keys()[0]]))
				nerrors[enst] = 1
				pass
	del mat
	print '%s valid and %s errors mapping ENST to ENSG'%(len(nvalid), len(nerrors))

	# gene level estimates: take the sum of transcript FPKMs (analogous to sum of TPMs divided by length but less accurate)
	genemat = []
	genenames = sorted(genematd.keys())
	for g in genenames:
		genemat += [map(ga.sumna, ga.transpose(genematd[g]))]
	
	matname = outdir+'%s.%s.%s.matrix.txt'%(fname,'gene',KAL_format)
	ga.printTable(genemat,['Gene']+labels, genenames, file=matname)
	pipeit('- Saved %s'%(matname),1)
	del genematd,genemat

# should also save data frames for this....
# incorporate map file into data frame
pipeit('- Saving data as data frame...',1)
tabname = outdir+'%s.%s.table.txt'%(fname,'transcript')
ga.printTable(tab_df,['SampleID','Label','Gene','Transcript','Info','TPM','Count'],file=tabname)
pipeit('- Saved %s'%(tabname),1)
del tab_df
pipeit('- Finished merging kallisto files.',1)

pipeit('DONE. EXITING.',1)




def main():
	pass

if __name__ == "__main__": main()