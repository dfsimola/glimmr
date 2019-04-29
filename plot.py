#! /usr/bin/env python

from utilities import *

plotWidth = 2048
plotHeight = 1536

keyLabel = ''
corMethod  = 'pearson'
togglepdf = False
UNDERSCORE = False

import Gnuplot as gp
gnuplotPath = '/GWD/bioinfo/projects/cb04/simola/bin.working/gnuplot-5.0.5/bin/gnuplot'
gp.GnuplotOpts.gnuplot_command = gnuplotPath

def setPath(path):
	global gnuplotPath
	gnuplotPath = path
	gp.GnuplotOpts.gnuplot_command = path

def getPath():
	global gnuplotPath
	return gnuplotPath


def setKey(spacing=1,slen=4):
	global keyLabel
	keyLabel = 'set key samplen %s spacing %s width 0 height 0'%(slen,spacing)

def getKey():
	return keyLabel

def setSize(width,height):
	global plotWidth, plotHeight
	
	plotWidth = width
	plotHeight = height
	return 1

def getSize():
	return (plotWidth,plotHeight)

def setCorMethod(arg): 
	global corMethod
	corMethod = arg
	return

def getCorMethod():
	return corMethod


def setUnderscore(value):
	global UNDERSCORE
	UNDERSCORE = value
	return UNDERSCORE

def getUnderscore():
	global UNDERSCORE
	return UNDERSCORE


# add cumulate to this
def hist(lst=[], nbins=100, categorical=False, bins=[], histdat=[], yfreq=False, logscale=False, shift=0, dtype=None, title='', xlabel=' ', xlabels=[], ylabel='Count', file='', legend='', legends=[], fontface='Helvetica', fontsize=8, fontsize_x=-1, fontsize_y=-1, reverseOrder=False, logbase=10, tittype='float', custom='', showstats=True, strict=1, xtics='', rotateLabels=-45, stagger=True, spacing=.9, ytotal=None, xticrounding=roundna(3), col=0, colors=[], labCoef=1., style='histogram', lw=1, smooth=None, xlabelspacing=5, pophead=0, grid='', error=[], ecolors='black', errorwidth=.05, ratio=1):
	"""This is a work in progres..."""
	
	# styles: histogram, outline, lines, linespoints
	
	if not len(lst) and not len(histdat): 
		print >> sys.stderr, 'plot.hist() error: no input data provided.'
		return -1
	
	legends = map(str,legends)
	
	# print 'underscore', getUnderscore()
	if getUnderscore()==False: legends = [l.replace('_','-') for l in legends]
	if title and getUnderscore()==False: title = title.replace('_','-')
	
	if len(lst) and type(lst[0]) != type([]): lst = [lst]
	if len(histdat) and (type(histdat[0]) != type([]) or len(histdat[0])==2): histdat = [histdat]
	if len(error) and type(error[0]) != type([]): error = [error]
	
	if categorical: xlabelspacing = 0
	else:
		if len(bins) > 0 and len(bins) < 20: xlabelspacing = 0
		if len(bins) > 100: xlabelspacing = 10
	
	
	try: import Gnuplot as gp
	except ImportError: return -1
	g = gp.Gnuplot()
	
	face = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	if file: fsize /= 2.0
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	
	g('set terminal pdf enhanced font "'+face+','+str(fsize)+'"')
	# g('set terminal postscript enhanced font "'+face+','+str(fsize)+'"')
	# g('set terminal pdf enhanced')
	format = GLOB_FORMAT
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file.split('/')[-1])
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - %s"'%(file))
		else:
			g('set terminal %s' % (format))
			g('set out "'+file+'"')
	
	if title: g.title(title)
	if yfreq and ylabel == 'Count': ylabel = 'Frequency'
	g.xlabel(xlabel)
	g.ylabel(ylabel)
	
	if logscale == 'x' or logscale == True: g('set logscale x')
	elif logscale == 'y': g('set logscale y')
	elif logscale == 'xy' or logscale: g('set logscale')
	elif logscale == True: g('set logscale y')
	
	
	lfd = len(lst)
	# scale font size if need be
	fticsize = minna([fsize, divna(fsize*labCoef,ln(float(lfd)))])
	# print 'lfd', lfd, 'labCoef', labCoef, 'fsize', fsize, fticsize, divna(fsize*labCoef,ln(float(lfd)))
	
	# trial
	# fticsize = minna([fsize, divna(fsize*3.3,ln(float(lfd))), divna(fsize*1.,ln(float(len(lst[0]))))])
	
	if fontsize_x == -1: fontsize_x = fticsize
	if fontsize_y == -1: fontsize_y = fticsize

	g('set xtics font "%s,%s"'%(face,fontsize_x))
	g('set ytics font "%s,%s"'%(face,fontsize_y))
	g('set ylabel font "%s,%s"'%(face,fsize))

	
	doreverse = ''
	if reverseOrder: doreverse = 'reverse'
	
	lst = map(lambda X: filterstr(nastr)(X), lst)
	
	hists = [None for i in range(len(lst))]
	if len(lst):
		for j in range(len(lst)):
			
			hists[j] = histogram(lst[j], nbins=nbins, bins=bins, yfreq=False, logscale=False, shift=shift, dtype=dtype, logbase=logbase, strict=strict, categorical=categorical)
	else:
		for i in range(len(histdat)): hists += [histdat[i]]
	
	
	# print 'hists'
	# print hists
	
	# pop some of the entries if users so desires
	if pophead > 0:
		for zz in range(pophead):
			for zj in range(len(hists)):
				hists[zj].pop(0)
	#
	
	
	if yfreq:
		for i in range(len(hists)):
			total = sum(vslice(hists[i],1))
			# total = len(lst[i])
			if ytotal: total = float(ytotal)
			if total > 0:
				hists[i] = map(lambda x: (x[0],x[1]/float(total)), hists[i])
	
	# set a margin for yrange
	if 1: #not yfreq:
		mnc = GLOB_INF; mxc = -GLOB_INF
		for his in hists:
			vsl = vslice(his,1)
			if not len(vsl): continue
			mn = min(vsl)
			mx = max(vsl)
			if mx > mxc: mxc = mx
			if mn < mnc: mnc = mn
		margin = .07*abs(mxc-mnc)
		
		if mnc == 0: mnc = 1e-6
		if mnc == GLOB_INF: mnc = '*'
		g('set yrange [%f:%f]' % (mnc,mxc+margin))
	# else:
		# g('set yrange [1e-6:1]')
	
	if legend and not legends: legends = [legend]
	_legends = map(str,legends[:])
	if len(legends) != len(hists):
		_legends = ['Dataset %d'%(i) for i in range(1,len(hists)+1)]
	
	dorotate = ''
	if rotateLabels!=None: dorotate = 'rotate by %s'%(rotateLabels)
	
	userlabels = False
	if xlabels: 
		userlabels = True
		
		# get xvalues for the data
		xlabelrange = map(str, range(len(xlabels)))
		# xlabelrange = []
		# if not categorical:
		# 	xlabelrange = map(str, xlabels)
		# else:
		# 	xlabelrange = map(str, range(len(xlabels)))
		
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), xlabelrange)
		
		# tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, xlabels))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		
		g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, xticstr, fontface, fsize))
	
	if not xlabels and not categorical:
		xlabels = map(lambda x: str(xticrounding(x)),sorted(vslice(hists[0],0)))
		# if bins:
		# 	xlabels = map(lambda x: str(x),sorted(vslice(hists[0],0)))
		# else:
		# 	xlabels = map(lambda x: str(xticrounding(x)),sorted(vslice(hists[0],0)))
	
	elif not xlabels:
		# get xlabels as the category keys
		xlabels = map(lambda x: str(x),vslice(hists[0],0))
		
	
	if len(xlabels) > 15 and not userlabels and not bins and not categorical:
		modfactor = nbins/float(50)
		xlabels = map(lambda x,y: x%modfactor==0. and y or '', range(len(xlabels)), xlabels)
	
	
	
	# re-coordinatize the bins themselves
	for j in range(len(hists)): hists[j] = map(lambda i: (i,hists[j][i][1]), range(len(hists[j])))
	if len(hists): bins = sorted(vslice(hists[0],0))
	
	if getUnderscore()==False: xlabels = [str(l).replace('_','-') for l in xlabels]
	
	if not userlabels:
		tmp = None
		if xticrounding and not categorical:
			try:
				tmp = zip(map(lambda x: '"'+str(xticrounding(x))+'"', xlabels), map(str, range(len(xlabels))))
			except ValueError:
				tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, range(len(xlabels))))
		else:
			tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, range(len(xlabels))))
			
			
		if xlabelspacing:
			tmp = map(lambda x: float(x[1]) % xlabelspacing == 0 and x or ('" "', x[1]), tmp)
			
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
			
		g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, xticstr, face, fticsize))
	
	# elif not xlabels:
	# 	g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, '', face, fticsize))
	
	
	#
	if stagger:
		# place hist bars adjacent to each other
		# width = (10/float(len(bins)))/float(len(hists))
		width = divna(1,float(len(hists)))
		boxwidth = spacing*width
		g('set boxwidth %s'%(boxwidth))
		
		fold = 0
		for i in range(1,len(hists)):
			fold += width
			hists[i] = map(lambda x: (x[0]+(spacing*fold),x[1]),hists[i])
			
			# add a little extra space after first column of each bin
			# for j in range(len(hists[i])):
			# 	if j == 0: fold += .05*width
			# 	hists[i][j] = (hists[i][j][0]+fold, hists[i][j][1])
				
			
	
	lsb = listRange(bins)
	binmargin = .05*lsb
	# binmargin = .1*lsb
	if lsb < binmargin: binmargin = 0
	bmn = min(bins)-binmargin
	# if bmn == 0: bmn = 1e-6
	if bmn == 0: bmn = -1
	
	binmarginhigh = binmargin+.61
	# binmarginhigh += sum([.13 for foobar in hists])
	g('set xrange [%f:%f] %s' % (bmn, max(bins)+binmarginhigh, doreverse))
	
	if grid != '':
		if grid == 'xy': g('set grid')
		else: g('set grid %s'%(grid))
	
	default = 'set size ratio %s; '%(ratio)
	g(default+custom)
	g(keyLabel)
	
	
	#
	# handle error bars
	ocol = col
	if len(colors) < len(lst) or colors == []:
		colors = [ocol+tmp+1 for tmp in range(len(lst))]
	# set styles once
	elwfactor = 1.
	if ecolors=='black':
		g('set style arrow 1 nohead lw %s lt %s'%(lw/elwfactor, -1))
		# print '- no colerr', 'set style arrow 1 nohead lw %s lt %s'%(lw/1.5, -1)
	elif ecolors == None:
		# print 'setting up ecolors'
		for li in xrange(len(hists)):
			g('set style arrow %d nohead lw %s lt %s'%(li+1, lw/elwfactor, colors[li]))
			
	elif type(ecolors) == type([]):
		for li in xrange(len(hists)):
			g('set style arrow %d nohead lw %s lt %s'%(li+1, lw/elwfactor, ecolors[li]))
			# print '- colerr', 'set style arrow 1 nohead lw %s lt %s'%(lw/1.5, -1)
	# --------------------------------------------
	
	
	# assume err is a list of ([low,high] pairs) or (low) values
	if len(error): 
		for li in range(len(hists)):
			dat = hists[li]
			err = error[li]
			te = type(err[0])
			DOUBLEUP = False
			if te==type(1) or te==type(1.) or te==type(nastr): 
				DOUBLEUP = True
				err = zip(err,err)
			
			# if logscale then make sure none are negative
			err = map(lambda x: (x[0]<0 and 1e-10 or x[0], x[1]<0 and 1e-10 or x[1]), err)
			
			# draw them as arrows with no heads at each x-value
			# dat should be a list of (x,y) pairs
			
			# x-values
			ds = vslice(dat,0)
			# compute average delta
			ds0 = ds[0]
			dx = avg(map(lambda x: abs(ds[x]-ds[x-1]), range(1,len(ds))))
			
			ass = 1
			if ecolors=='black':
				ass = 1
			elif ecolors == None:
				ass = li+1
			elif type(ecolors) == type([]):
				ass = li+1
			
			logfunc = lambda x: log(x)+.0001
			for zi in range(len(dat)):
				xlow = dat[zi][0] - dx*float(errorwidth)
				xhi = dat[zi][0] + dx*float(errorwidth)
				try:
					if logscale == 'x' or logscale=='xy':
						pass
					
					# not implemented
					if logscale == 'x' or logscale=='xy':
						# draw the verical line
						g('set arrow from %s,%s to %s,%s as %d'%(logfunc(dat[zi][0]), dat[zi][1]-err[zi][0], logfunc(dat[zi][0]), dat[zi][1]+err[zi][1],ass))
						# draw the bars
						g('set arrow from %s,%s to %s,%s as %d'%(logfunc(xlow), dat[zi][1]-err[zi][0], logfunc(xhi), dat[zi][1]-err[zi][0],ass))
						g('set arrow from %s,%s to %s,%s as %d'%(logfunc(xlow), dat[zi][1]+err[zi][1], logfunc(xhi), dat[zi][1]+err[zi][1],ass))
					
					else:
						
						if DOUBLEUP:
							# draw the verical line
							g('set arrow from %s,%s to %s,%s as %d'%(dat[zi][0], dat[zi][1]-err[zi][0], dat[zi][0], dat[zi][1]+err[zi][1],ass))
							# draw the bars
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, dat[zi][1]-err[zi][0], xhi, dat[zi][1]-err[zi][0],ass))
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, dat[zi][1]+err[zi][1], xhi, dat[zi][1]+err[zi][1],ass))
						else:
							# don't offset the error bars
							# draw the verical line
							g('set arrow from %s,%s to %s,%s as %d'%(dat[zi][0], err[zi][0], dat[zi][0], err[zi][1],ass))
							# draw the bars
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, err[zi][0], xhi, err[zi][0],ass))
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, err[zi][1], xhi, err[zi][1],ass))
							
						
				except TypeError: pass
	
	
	
	
	# # add a blank entry so last bin isn't cut off
	# if categorical and len(hists)>1:
	# 	# print 'test', topbin, hists[:10]
	# 	for i in range(len(hists)):
	# 		topbin = max([x[0] for x in hists[i]])
	# 		# print 'i', i, hists[i]
	# 		hists[i] += [(topbin+1,0)]
	# 		print hists[i]
	
	# add stats
	for i in range(len(hists)):
		av = med = var = nastr
		nfilt = 0
		
		try:
			av = avg(lst[i]); med = median(lst[i]); var = variance(lst[i]); nfilt = len(lst[i])
		except IndexError: pass
		except ValueError: pass # if categorical
		
		tit = ''
		if showstats:
			try:
				if tittype == 'float': 
					tit = 'Avg=%.2f Median=%.2f Var=%.2f (n=%d)' % (av, med, var, nfilt)
				else:
					tit = 'Avg=%.4e Median=%.4e Var=%.4e (n=%d)' % (av, med, var, nfilt)
			except TypeError:
				dequote = lambda x: str(x).strip('"')
				tit = 'Avg=%s Median=%s Var=%s (n=%d)' % (dequote(r4(av)), dequote(r4(med)), dequote(r4(var)), nfilt)
			_legends[i] = (len(_legends[i]) and tit+', '+_legends[i]) or tit
	
	if not colors or len(colors)<len(hists):
		colors = range(1,len(hists)+1)
	
	
	if len(style.split(','))>1:
		jc = 0
		for sty in style.split(','):
			if sty == 'histogram':
				g('set style data boxes')
				g('set style fill solid border -1')
			elif sty == 'outline':
				pass
			# else:
			# 	g('unset style data')
			# 	g('unset style fill')
			
			dat = []
			for i in range(len(hists)):
				if sty == 'histogram' or sty == 'outline':
					dat += [ gp.Data(hists[i], with_='boxes lt %s lw %s'%(colors[i],lw)) ]
				elif sty == 'lines' or sty == 'linespoints':
					tmp = hists[i]
					if smooth: tmp = map(lambda x,y: [x[0],y], tmp, kernelsmooth([x[1] for x in hists[i]], w=smooth))
					dat += [ gp.Data(tmp, with_='%s lt %s lw %s'%(sty, colors[i]+1, lw+2),title=_legends[i]) ]
			
			if jc == 0: g.plot(*dat)
			else: 
				g('set out')
				g('set terminal x11 enhanced font "'+face+','+str(fsize)+'"')
				format = GLOB_FORMAT
				if file:
					crap,format,crap = re.split('^.+\.(.+)$', file.split('/')[-1])
					if format == 'ps': format = 'postscript'
					if togglepdf:
						g('set terminal postscript enhanced')
						g('set out "| ps2pdf - %s"'%(file))
					else:
						g('set terminal %s enhanced' % (format))
						g('set out "'+file+'"')
				g.replot(*dat)
			jc += 1
			
	else:
		if style == 'histogram':
			g('set style data boxes')
			g('set style fill solid border -1')
			
		dat = []
		for i in range(len(hists)):
			
			tmp = hists[i]
			if smooth: tmp = map(lambda x,y: [x[0],y], tmp, kernelsmooth([x[1] for x in hists[i]], w=smooth))
			
			if style == 'histogram' or style == 'outline':
				dat += [ gp.Data(tmp, with_='boxes lt %s lw %s'%(colors[i],lw),title=_legends[i]) ]
			elif style == 'lines' or style == 'linespoints':
				dat += [ gp.Data(tmp, with_='%s lt %s lw %s'%(style, colors[i], lw),title=_legends[i]) ]
		g.plot(*dat)
	g('set out')
	


def scatterOne(X,Y=[], xlabel='x-axis',ylabel='y-axis', title='', logscale='', style='points', file='', reverseOrder=False, ps=.5, lw=2, col=1, rgbcol='', pt=1, legend='', key=True, xlabels=[], ylabels=[], rotateLabels=False, inkey=corMethod,custom='', pointLabel=None, pointLabels=[], colors=[], fontface='Helvetica', fontsize=0, blanklegend=False):
	try: import Gnuplot as gp
	except ImportError: return -1
	
	dat = None
	if not len(Y): dat = X
	else: dat = zip(X,Y)
	
	fontface = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	if file: fsize /= 3.0
	
	g = gp.Gnuplot()
	g('set terminal pdf enhanced font "'+fontface+','+str(fsize)+'"')
	format = GLOB_FORMAT
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file)
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - \"%s\""'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
	g.ylabel(ylabel)
	g.xlabel(xlabel)
	if title: g.title(title)
	
	ret = corWithPvalue(dat, method=inkey, tails='2', p='approx')
	if ret['P'] == nastr: ret['P'] = 1
	r = ret['r']
	d2 = pairedfilterstr(nastr)(dat[:])
	nfilt = len(d2)
	
	if logscale == 'x': g('set logscale x')
	elif logscale == 'y': g('set logscale y')
	elif logscale == 'xy' or logscale: g('set logscale')
	
	if reverseOrder: g('set xrange [] reverse')
	
	
	legstrchar = 'r'
	if inkey == 'spearman': legstrchar = 's'
	
	dorotate = ''
	if rotateLabels: dorotate = 'rotate by -75'
	if xlabels:
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, xlabels))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		
		g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, xticstr, fontface, fsize))
	if ylabels:
		g('set ytics nomirror '+'('+','.join(ylabels)+') font "%s,%s"'%(fontface, fsize))
	
	g(custom)
	
	if pointLabel:
		# apply same label to all points
		text = pointLabel; color = rgbcol
		if list(pointLabel) == pointLabel:
			text,color = pointLabel[0],pointLabel[1]
		for li in range(len(d2)):
			xco = d2[li][0]
			yco = d2[li][1]
			lab = 'set label "%s" at %s,%s center front font "%s,%s"'%(text,xco,yco,fontface,fsize)
			if len(color): lab += 'textcolor rgb \"%s\"'%(color)
			g(lab)
	
	if pointLabels:
		# draw a text label at each point
		for li in range(len(d2)):
			text = pointLabels[li]
			xco = d2[li][0]
			yco = d2[li][1]
			g('set label "%s" at %s,%s center front font "%s,%s"'%(text,xco,yco,fontface,fsize))
	
	if colors:
		# draw a color label at each point
		for li in range(len(d2)):
			acol = colors[li]
			xco = d2[li][0]
			yco = d2[li][1]
			
			if xco == 0: xco += 1e10
			if yco == 0: yco += 1e10
			
			g('set label "%s" at %s,%s textcolor lt %s center front font "%s,%s"'%(pt,xco,yco,acol,fontface,fsize))
	
	
	legstr = ''
	if blanklegend == True: 
		legstr = ''
		if legend: legstr = legend
	else:
		try:
			legstr = '{/Helvetica-Oblique %s} = %.2f, {/Helvetica-Oblique P} < %.2e (n=%d)'%(legstrchar,ret['r'],ret['P'],nfilt)
		except TypeError:
			legstr = '(n=%d)'% (nfilt)
		if legend: legstr += ', '+legend
	
	try:
		if rgbcol:
			if key:
				g.plot(gp.Data(d2, with_=style+" lw %d pt '%s' ps %f rgb \"%s\"" % (lw, pt, ps,rgbcol), title=legstr))
			else:
				g.plot(gp.Data(d2, with_=style+" lw %d pt '%s' ps %f rgb \"%s\"" % (lw, pt, ps, rgbcol)))
		else:
			if key:
				
				g.plot(gp.Data(d2, with_=style+" lt %d lw %d pt '%s' ps %f" % (col, lw, pt, ps), title=legstr))
			else:
				g.plot(gp.Data(d2, with_=style+" lt %d lw %d pt '%s' ps %f" % (col, lw, pt, ps)))
	except AssertionError:
		print 'utilities.scatterplot: There was an error in plotting this data.'
	if file: g('set out')
	return g


def scatterRE(dat, replot, style='points', ps=.5, lw=2, col=1, rgbcol='', pt=1, legend='', inkey=corMethod,key=True, file=''):
	
	assert replot, "Sorry, scatterREplot needs a Gnuplot object as an argument."
	import Gnuplot as gp
	g = replot
	
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file)
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - \"%s\""'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
	
	ret = corWithPvalue(dat, method=inkey, tails='2', p='approx')
	if ret['P'] == nastr: ret['P'] = 1
	r = ret['r']
	d2 = pairedfilterstr(nastr)(dat[:])
	nfilt = len(d2)
	
	legstrchar = 'r'
	if inkey == 'spearman': legstrchar = 's'
	
	try:
		if rgbcol:
			if key:
				legstr = '{/Helvetica-Oblique %s} = %.2f, {/Helvetica-Oblique P} < %.2e (n=%d)'%(legstrchar,ret['r'],ret['P'],nfilt)
				# legstr = "R=%.4f (n=%d)" % (r, nfilt)
				if legend: legstr += ', '+legend
				g.replot(gp.Data(d2, with_=style+" lw %d pt '%s' ps %f rgb \"%s\"" % (lw, pt, ps,rgbcol), title=legstr))
			else:
				g.replot(gp.Data(d2, with_=style+" lw %d pt '%s' ps %f rgb \"%s\"" % (lw, pt, ps, rgbcol)))
		else:
			if key:
				# if r == nastr: r = -0
				# legstr = "r=%.4f,  (n=%d)" % (r, nfilt)
				try:
					legstr = '{/Helvetica-Oblique %s} = %.2f, {/Helvetica-Oblique P} < %.2e (n=%d)'%(legstrchar,ret['r'],ret['P'],nfilt)
				except TypeError:
					legstr = '(n=%d)'% (nfilt)
				if legend: legstr += ', '+legend
				g.replot(gp.Data(d2, with_=style+" lw %d pt '%s' ps %f" % (lw, pt, ps), title=legstr))
			else:
				g.replot(gp.Data(d2, with_=style+" lw %d pt '%s' ps %f" % (lw, pt, ps)))
	except AssertionError:
		print 'utilities.scatterplot: There was an error in plotting this data.'
	if file: g('set out')
	return g



def scatter(lst, xlabel=' ', ylabel=' ', title='', logscale='', style='points', styles=[], col=1, rgbcol='', colors=[], reverseOrder=False, ps=1, pointSizes=[], lineWeights=[], lw=2,  pt=.9, pointTypes=[], legend=None, legends=[], xlabels=[], ylabels=[], rotateLabels=None, key=True, inkey=None, replot=None, custom='', fontface='Helvetica', fontsize=0, legendfontsize=0, showstats=True, apre=False, file='', error=[], ls='', pointLabel=None, pointLabels=[], labelOffset=[0,0], labelidx=0, errorwidth=.05, labelsize=None, grid='', xlabelspacing=None, noleg=[], ecolors=None, errorbarwidth=None, regress=False, showslope=True, smooth=0, ratio=1):
	"""Difference between this and scatter plot is this takes a list, containing 1 or more lists of x,y pairs."""
	try: import Gnuplot as gp
	except ImportError: 
		print >> sys.stderr, 'Error: Cannot import Gnuplot.py'
		return -1
	ocol = col
	
	global nastr

	if getUnderscore()==False and len(xlabels): 
		xlabels = [str(l).replace('_','-') for l in xlabels]
	if getUnderscore()==False and len(legends): 
		legends = [str(l).replace('_','-') for l in legends]
	if title and getUnderscore()==False: title = title.replace('_','-')
	
	
	if len(lst):
		if type(lst[0]) != type([]): lst = [lst]
		elif type(lst[0]) == type([]) and len(lst[0]) and type(lst[0][0]) != type([]) and type(lst[0][0]) != type(()): lst = [lst]
	
	if len(error) and type(error[0]) != type([]): error = [error]
	
	
	# delete empty data sets
	lst2 = []
	for tmp in lst:
		if len(tmp)>0:
			lst2 += [tmp]
	lst = lst2
	
	if ps == 1 and len(lst[0]) > 2500:
		ps = .75
		
	
	if type(pt) == type('asdf'):
		pt = '"'+pt+'"'
	ptdef = pt
	psdef = ps
	
	fontface = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	flabelsize = GLOB_FONT_SIZE
	
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	if file: fsize /= 3.0
	if legendfontsize == 0: legendfontsize = fsize
	if labelsize: flabelsize = labelsize
	if file: flabelsize/=3.0
	
	g = gp.Gnuplot()
	if replot: g = replot
	g('set terminal pdf enhanced font "'+fontface+','+str(fsize)+'"')
	format = GLOB_FORMAT
	if file and not replot:
		format = 'pdf'
		
		# does file specify format in extension
		format = getSuffix(file)
		
		# try: crap,format,crap = re.split('^.+\.(.+)$', file)
		# except ValueError: file = file+'.%s'%(format)
			
		if format == 'ps': format = 'postscript'
		
		if format == 'png':
			g('set terminal %s enhanced size %s, %s' % (format, plotWidth, plotHeight))
			g('set out "'+file+'"')
		elif togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - %s"'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
	
	lfd = len(lst)
	
	# scale font size if need be
	fticsize = minna([fsize, divna(fsize*3.3,ln(float(lfd)))])
	
	g('set xlabel font "%s,%s"'%(face,fticsize))
	g('set xtics font "%s,%s"'%(face,fticsize))
	g('set ytics font "%s,%s"'%(face,fticsize))
	g('set ylabel font "%s,%s"'%(face,fticsize))
	g.ylabel(ylabel)
	g.xlabel(xlabel)
	if title: g.title(title)
	
	if logscale == 'x': g('set logscale x')
	elif logscale == 'y': g('set logscale y')
	elif logscale == 'xy' or logscale: g('set logscale')
	
	if reverseOrder: g('set xrange [] reverse')
	
	
	
	
	# --------------------------------------------
	if len(colors) < len(lst) or colors == []:
		colors = [ocol+tmp+1 for tmp in range(len(lst))]
	
	# set styles once
	elwfactor = 5.
	if ecolors=='black':
		g('set style arrow 1 nohead lw %s lt %s'%(lw/elwfactor, 1))
		# print '- no colerr', 'set style arrow 1 nohead lw %s lt %s'%(lw/1.5, -1)
	elif ecolors == None:
		# print 'setting up ecolors'
		for li in xrange(len(lst)):
			g('set style arrow %d nohead lw %s lt %s'%(li+1, lw/elwfactor, colors[li]))
			
	elif type(ecolors) == type([]):
		for li in xrange(len(lst)):
			g('set style arrow %d nohead lw %s lt %s'%(li+1, lw/elwfactor, ecolors[li]))
			# print '- colerr', 'set style arrow 1 nohead lw %s lt %s'%(lw/1.5, -1)
	# --------------------------------------------
	
	
	dorotate = ''
	if rotateLabels != None: dorotate = 'rotate by %s'%(rotateLabels)
	if xlabels != None:
		# get xvalues for the data
		xlabelrange = []
		vls = vslice(lst[0],0)
		
		if not xlabelspacing:
			xlabelrange = map(str, vls)
		else:
			try: xlabelrange = map(str, xrange(0,max(vls),xlabelspacing))
			except ValueError: xlabelrange = map(str, vls)
			
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), xlabelrange)
		# tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, xlabels))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		g('set xtics nomirror %s %s font "%s,%s";'%(dorotate, xticstr, fontface, fsize))
	else:
		g('set xtics (-5,-5)')
		# g('set xtics ""')
	
	
	if ylabels:
		tmp = zip(map(lambda x: '"'+str(x)+'"', ylabels), map(str, ylabels))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		yticstr = '('+','.join(map(str, tmp2))+')'
		# print xticstr
		g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, yticstr, fontface, fsize))
		
	if grid != '':
		if grid == 'xy': g('set grid')
		else: g('set grid %s'%(grid))
	
	if len(styles) < len(lst):
		styles = [style for tmp in range(len(lst))]
	
	if len(pointSizes) < len(lst):
		# print >> sys.stderr, 'Warning: |pointSizes| < |data|', len(pointSizes), len(lst)
		pointSizes = [ps for tmp in range(len(lst))]
	
	#
	if len(lineWeights) < len(lst):
		lineWeights = [lw for tmp in range(len(lst))]
	
	if legend and not len(legends):
		legends = [legend]
	
	
	
	# for li in reversed(range(len(lst))):
	
	# print 'list', len(lst), len(pointTypes), len(pointSizes)
	# regressdata = []
	
	# get min and max x data
	
	# print 'HEY'
	# for x in lst:
	# 	print 'test', vslice(x,0)
	# 	map(floatna,flatten(vslice(x,0)))
	
	xdat = map(floatnacrap, flatten([vslice(x,0) for x in lst]))
	minx = minna(xdat)
	maxx = maxna(xdat)
	# print 'hey', minx, maxx
	
	for li in xrange(len(lst)):
		
		style = styles[li]
		col = colors[li]
		lw = lineWeights[li]
		try: pt = pointTypes[li]
		except IndexError: pt = psdef
		try: ps = type(pointSizes[li])==type('a') and '"'+pointSizes[li]+'"' or str(pointSizes[li])
		except IndexError: ps = ptdef
		
		# print li, 'style:', style, col, ps, pt, lw
		
		if li > 0:
			if not replot: 
				g('set out')
				if file: 
					g('set out "'+file+'"')
		
		dat = lst[li]
		
		# handle error bars
		# assume err is a list of ([low,high] pairs) or (low) values
		err = None
		if len(error): 
			err = error[li]
			te = type(err[0])
			DOUBLEUP = False
			if te==type(1) or te==type(1.) or te==type(nastr): 
				DOUBLEUP = True
				err = zip(err,err)
			
			# if logscale then make sure none are negative
			err = map(lambda x: (x[0]<0 and 1e-10 or x[0], x[1]<0 and 1e-10 or x[1]), err)
			
			# draw them as arrows with no heads at each x-value
			# dat should be a list of (x,y) pairs
			
			# x-values
			ds = vslice(dat,0)
			# compute average delta
			ds0 = ds[0]
			if not errorbarwidth:
				dx = avg(map(lambda x: abs(ds[x]-ds[x-1]), range(1,len(ds))))
			else:
				dx = errorbarwidth
			if dx == nastr: dx = 0
			
			ass = 1
			if ecolors=='black':
				ass = 1
			elif ecolors == None:
				ass = li+1
			elif type(ecolors) == type([]):
				ass = li+1
			
			logfunc = lambda x: log(x)+.0001
			for zi in range(len(dat)):
				xlow = dat[zi][0] - dx*float(errorwidth)
				xhi = dat[zi][0] + dx*float(errorwidth)
				try:
					if logscale == 'x' or logscale=='xy':
						pass
					
					# not implemented
					if logscale == 'x' or logscale=='xy':
						# draw the verical line
						g('set arrow from %s,%s to %s,%s as %d'%(logfunc(dat[zi][0]), dat[zi][1]-err[zi][0], logfunc(dat[zi][0]), dat[zi][1]+err[zi][1],ass))
						# draw the bars
						g('set arrow from %s,%s to %s,%s as %d'%(logfunc(xlow), dat[zi][1]-err[zi][0], logfunc(xhi), dat[zi][1]-err[zi][0],ass))
						g('set arrow from %s,%s to %s,%s as %d'%(logfunc(xlow), dat[zi][1]+err[zi][1], logfunc(xhi), dat[zi][1]+err[zi][1],ass))
					
					else:
						
						if DOUBLEUP:
							# draw the verical line
							g('set arrow from %s,%s to %s,%s as %d'%(dat[zi][0], dat[zi][1]-err[zi][0], dat[zi][0], dat[zi][1]+err[zi][1],ass))
							# draw the bars
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, dat[zi][1]-err[zi][0], xhi, dat[zi][1]-err[zi][0],ass))
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, dat[zi][1]+err[zi][1], xhi, dat[zi][1]+err[zi][1],ass))
						else:
							# don't offset the error bars
							# draw the verical line
							g('set arrow from %s,%s to %s,%s as %d'%(dat[zi][0], err[zi][0], dat[zi][0], err[zi][1],ass))
							# draw the bars
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, err[zi][0], xhi, err[zi][0],ass))
							g('set arrow from %s,%s to %s,%s as %d'%(xlow, err[zi][1], xhi, err[zi][1],ass))
							
						
				except TypeError: pass
				
				# try:
				# 	# draw the verical line
				# 	g('set arrow from %s,%s to %s,%s as %s'%(dat[zi][0], dat[zi][1]-err[zi][0], dat[zi][0], dat[zi][1]+err[zi][1]),ass)
				# 	# draw the bars
				# 	g('set arrow from %s,%s to %s,%s as %d'%(xlow, dat[zi][1]-err[zi][0], xhi, dat[zi][1]-err[zi][0]),ass)
				# 	g('set arrow from %s,%s to %s,%s as %d'%(xlow, dat[zi][1]+err[zi][1], xhi, dat[zi][1]+err[zi][1]),ass)
				# except TypeError: pass
		
		
		#################
		if pointLabel:
			# apply same label to all points
			text = pointLabel; color = rgbcol
			if list(pointLabel) == pointLabel:
				text,color = pointLabel[0],pointLabel[1]
			for li in range(len(dat)):
				xco = dat[li][0]
				yco = dat[li][1]
				lab = 'set label "%s" at %s,%s center front font "%s,%s"'%(text,xco,yco,fontface,flabelsize)
				if len(color): lab += 'textcolor rgb \"%s\"'%(color)
				g(lab)
				
		if pointLabels:
			if type(pointLabels[li]) == type([]):
				for li2 in range(len(dat)):
					text = pointLabels[li][li2]
					xco = dat[li2][0]
					yco = dat[li2][1]
					g('set label "%s" at %s,%s center front font "%s,%s"'%(text,xco+labelOffset[0],yco+labelOffset[1],fontface,flabelsize))
				
			elif li == labelidx:
				# draw a text label at each point
				for li2 in range(len(dat)):
					text = pointLabels[li2]
					xco = dat[li2][0]
					yco = dat[li2][1]
					g('set label "%s" at %s,%s center front font "%s,%s"'%(text,xco+labelOffset[0],yco+labelOffset[1],fontface,flabelsize))
			# else:
				# g('set label')
				
		# if colors:
		# 	# draw a color label at each point
		# 	for li in range(len(dat)):
		# 		acol = colors[li]
		# 		xco = dat[li][0]
		# 		yco = dat[li][1]
		# 
		# 		if xco == 0: xco += 1e10
		# 		if yco == 0: yco += 1e10
		# 
		# 		g('set label "%s" at %s,%s textcolor lt %s center front font "%s,%s"'%(pt,xco,yco,acol,fontface,fsize))
		
		
		
		#################
		
		legstrchar = 'r'
		if inkey == 'spearman' or getCorMethod()=='spearman': legstrchar = 's'
		
		d2 = pairedfilterstr(nastr)(dat[:])
		nfilt = len(d2)
		
		# regression slope
		m = b = 0
		
		if len(legends) > 0 and len(legends) < len(lst): 
			# add blanks
			while len(legends) < len(lst): legends += ['']
			# print 'Error: too few legends for provided data sets.'
		
		if showstats == False or (showstats > 1 and li+1 > showstats):
			legstr = ''
			if legends and legends[li]: legstr = legends[li]
			g('set key')
		elif showstats == 1 or (showstats > 1 and li+1 <= showstats):
			ret = corWithPvalue(dat, method=inkey and inkey or getCorMethod(), tails='2', p='approx')
			if ret['P'] == nastr: ret['P'] = 1
			r = ret['r']
			
			
			try:
				res = linearRegression(dat)
				m = res['m']; b = res['b']
			except ZeroDivisionError: pass
			except TypeError: m = nastr; b = nastr
			# regressdata += [[b,m,col]]
			
			try:
				if regress and showslope:
					legstr = '{/Helvetica-Oblique m} = %s, {/Helvetica-Oblique b} = %s; {/Helvetica-Oblique %s} = %.2f, {/Helvetica-Oblique P} < %.2e (n=%d)'%(r2(m),r2(b),legstrchar,ret['r'],ret['P'],nfilt)
				else:
					legstr = '{/Helvetica-Oblique %s} = %.2f, {/Helvetica-Oblique P} < %.2e (n=%d)'%(legstrchar,ret['r'],ret['P'],nfilt)
				
			except TypeError:
				legstr = '(n=%d)'% (nfilt)
			
			if legends and legends[li]: legstr += ', '+str(legends[li])
		
		legstr = "{/%s=%s %s}"%(fontface,legendfontsize,legstr)
		
		if noleg and noleg[li]==1: legstr = ''
		
		if style == 'linespoints' or style == 'points':
			style += ' pt %s ps %s'%(pt, ps)
		
		g(keyLabel)
		
		default = 'set size ratio %s; '%(ratio)
		if li == len(lst)-1: g(default+custom) # only apply this once at the end
		
		if regress:
			if m == 0 and b == 0:
				res = linearRegression(dat)
				m = res['m']; b = res['b']
			# add a regression line here
			# g.replot(gp.Data(d3), with_='lines lt %d lw %d'%(col, lw))
			# print 'm', m, minx
			yl = b + m*minx
			yh = b + m*maxx
			g('set arrow nohead from %s,%s to %s,%s front lt %d'%(minx,yl,maxx,yh, col))
			# print 'REGRESS', 'm',m,'b',b, minx, maxx, 'set arrow nohead from %s,%s to %s,%s front lt %d'%(minx,yl,maxx,yh, col)
		
		try:
			if rgbcol:
				if key:
					if li > 0: g.replot(gp.Data(d2, with_=style+" lw %d rgb \"%s\"" % (lw, rgbcol), title=legstr))
					else:      g.plot(gp.Data(d2, with_=style+" lw %d rgb \"%s\"" % (lw, rgbcol), title=legstr))
				else:
					if li > 0: g.replot(gp.Data(d2, with_=style+" lw %d rgb \"%s\"" % (lw, rgbcol)))
					else:      g.plot(gp.Data(d2, with_=style+" lw %d rgb \"%s\"" % (lw, rgbcol)))
			else:
				if key:
					if li > 0: g.replot(gp.Data(d2, with_=style+" lt %d lw %d" % (col, lw), title=legstr))
					else:      g.plot(gp.Data(d2, with_=style+" lt %d lw %d" % (col, lw), title=legstr))
				else:
					if li > 0: g.replot(gp.Data(d2, with_=style+" lt %d lw %d" % (col, lw)))
					else:      g.plot(gp.Data(d2, with_=style+" lt %d lw %d" % (col, lw)))
			
		except AssertionError:
			print 'plot.scatter(): There was an error in plotting this data.'
	
	if file and not replot and apre != True: g('set out')
	return g


multiscatter = scatter

def heatmap(dat, title='', xlabel=' ', ylabel=' ', labels='', xlabels='', ylabels='', file='', xcells=0, ycells=0, z=4, color=-1, contour=0, pt=.1, type='2D', cblabel='Value', formula=[], ratio=1, na=True, invertXY=True, imgsize=[640,480], cbrange=[], custom='', fontface='Helvetica', fontsize=0, trimylabel=False):
	"""Data must be provided in the form of triples (x,y,z), where x and y are grid coordinates and z is the value.
	
	how to get a good heatmap
	
	ut.heatmap(trips,labels=r,custom='set ytics rotate by -60 font "Helvetica,12"; set xtics rotate by -30 font "Helvetica,12";unset border', cblabel="Spearman's s")
	
	ut.heatmap(trips,labels=r,custom='set ytics norotate font "Helvetica,12"; set xtics norotate font "Helvetica,12"; unset border; set margin 15', cblabel="Spearman's s")
	
	
	"""
	
	try: import Gnuplot as gp
	except ImportError: return -1
	face = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	g = gp.Gnuplot()
	g('set terminal pdf enhanced font "'+face+','+str(fsize)+'"')
	format = GLOB_FORMAT
	
	
	crap,format,crap = re.split('^.+\.(.+)$', file)
	if format == 'ps': format = 'postscript'
	
	if file:
		g('set terminal %s enhanced' % (format))
		if format == 'png':
			face = 'Verdana'
			g('set terminal png nocrop enhanced size %d,%d' % (imgsize[0], imgsize[1]))
			
		g('set out "'+file+'"')
	else:
		# custom = 'set ytics rotate font "Helvetica,12"; set xtics norotate font "Helvetica,12"; unset border'
		margin = 5
		custom += '; set bmargin %d; set tmargin %d; set lmargin %d; set rmargin %d' % (margin, margin, margin, margin)
		
	# change these to helvetica, fsize
	if ylabel: g.ylabel(ylabel)
	if xlabel: g.xlabel(xlabel)
	if title: g.title(title)
	
	if fontface: face = fontface
	
	if file and format == 'pdf': fsize /= 3.0
	if fontsize: fsize = fontsize
	
	if ratio == 1: g('set size square')
	else: g('set size ratio %f' % (ratio))
	
	if labels:
		tmp = zip(map(lambda x: '"'+str(x)+'"', labels), map(str, range(len(labels))))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		
		tmplab = labels[:]
		if trimylabel: tmplab[0] = '     '
		tmp = zip(map(lambda x: '"'+str(x)+'"', tmplab), map(str, range(len(tmplab))))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		yticstr = '('+','.join(map(str, tmp2))+')'
	
	if xlabels:
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, range(len(xlabels))))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
	
	if ylabels:
		tmplab = ylabels[:]
		tmp = zip(map(lambda x: '"'+str(x)+'"', tmplab), map(str, range(len(tmplab))))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		yticstr = '('+','.join(map(str, tmp2))+')'
	
	if labels or xlabels or ylabels:
		if file and format == 'pdf': 
			g('set xtics border in scale 0,0 nomirror rotate by -30 offset character -1, 0, 0 '+xticstr+' font "%s,%s"'%(face, fsize))
			g('set ytics border in scale 0,0 nomirror rotate by -38 offset character 0, 0, 0 '+yticstr+' font "%s,%s"'%(face, fsize))
		elif not file:
			g('set xtics border out scale 0,0 nomirror rotate by -30 offset character 1.3, -.7, 0 '+xticstr+' font "%s,%s"'%(face, fsize))
			shifty = 0
			if trimylabel: shifty = -1
			g('set ytics border out scale 0,0 nomirror rotate by -38 offset character -2.6, %d, 0 '%(shifty)+yticstr+' font "%s,%s"'%(face, fsize))
		else:
			g('set xtics border out scale 0,0 nomirror rotate by -30 offset character 1.3, -.7, 0 '+xticstr)
			shifty = 0
			if trimylabel: shifty = -1
			g('set ytics border out scale 0,0 nomirror rotate by -38 offset character -2.6, %d, 0 '%(shifty)+yticstr)
		g('set ztics border out scale 0,0 nomirror norotate  offset character 0, 0, 0')
		
	if na: 
		dat = triplefilterstr(nastr)(dat)
	
	xflat = [dat[i][0] for i in range(len(dat))]
	yflat = [dat[i][1] for i in range(len(dat))]
	zflat = [dat[i][2] for i in range(len(dat))]
	
	# import utilities as ut; a = [[0,0,1],[0,1,3],[1,0,3],[1,1,2.5]]; ut.heatmap(a, type='2D')
	g("""unset key
	set view map
	set style data linespoints
	set nocbtics
	set tics scale 0
	""")
	
	if len(formula) == 3:
		g('set palette rgbformula %d,%d,%d' % (formula[0], formula[1], formula[2]))
	
	# g('set grid back')
	g('set border 4096 front lw 1')
	
	if type == '3D':
		g('set pm3d')
		if xcells == 0: xcells = len(unique(xflat))
		if ycells == 0: ycells = len(unique(yflat))
		g("set dgrid3d %d,%d,%d" % (xcells, ycells, z))
		if contour: g('set contour surface')
		g(custom)
		g.splot(gp.Data(dat,with_='points lt %d pt %f' % (color, pt)))
	else:
		g('set urange [ * : * ] noreverse nowriteback')
		g('set vrange [ * : * ] noreverse nowriteback')
		if invertXY:
			g('set xrange [%d:%d] noreverse nowriteback' % (0, len(unique(yflat))-.5))
			g('set yrange [%d:%d] noreverse nowriteback' % (0, len(unique(xflat))-.5))
		else:
			g('set xrange [%d:%d] noreverse nowriteback' % (0, len(unique(xflat))-.5))
			g('set yrange [%d:%d] noreverse nowriteback' % (0, len(unique(yflat))-.5))
		g('set zrange [ * : * ] noreverse nowriteback')
		
		# need a set cb position
		
		if format == 'png' and file:
			g('set cblabel "%s" offset -2' % (cblabel))
		else:
			g('set cblabel "%s" offset -2 font "%s,%d"' % (cblabel, face, fsize))
		
		if not file: g('set colorbox user origin .75,.245 size .03, .635')
		if len(cbrange) == 2: g('set cbrange [%f:%f]' % (cbrange[0], cbrange[1]))
		else: g('set cbrange [*:*]')
		# g('set cbrange [%d:%d]' % (min(zflat), max(zflat)))
		g(custom)
		g.splot(gp.Data(dat, with_='image'))
	
	return g


def Rheatmap(mat=None, rows=[], cols=[], trimrow=False, tmpfile='', file='myheatmap.pdf', col='brewer.pal(256,"Blues")', outdir='/tmp/', scale='"none"', rmtmp=False, title='', xlabel='', ylabel='', fontsize=None, colfontsize=1, rowfontsize=.5, key=True, version='builtin', cluster=None, width=2048, height=1536, custom="", invert=True, clusterMethod='average', Rversion='R',scriptname='rscpt.heatmap.R'):
	# library(RColorBrewer)
	# pal = brewer.pal(256,"Blues")
	# pal = brewer.pal(1024,"Oranges")
	# # pal = brewer.pal(256,"BuGn")
	# # pal = brewer.pal(256,"BrBG")
	# # pal = brewer.pal(256,"Accent")
	# # pal = brewer.pal(1024,"PiYG")
	# # pal = brewer.pal(256,"PuOr")
	
	if fontsize: 
		colfontsize = rowfontsize = fontsize
		
	# cluster by correlation
	if cluster == 'default' or cluster == 'euclidean':
		cluster = "Rowv=TRUE,Colv=TRUE"
		# invert = False
	elif cluster == 'row' or cluster == 'x.euclidean':
		cluster = "Rowv=TRUE,Colv=NA"
	elif cluster == 'col' or cluster == 'y.euclidean':
		cluster = "Rowv=NA,Colv=TRUE"
	elif cluster == None or cluster.lower() == 'none':
		cluster = "Rowv=NA,Colv=NA"
		# invert = False
	elif cluster == 'correlation':
		# cluster = "Rowv=as.dendrogram(hclust(dist(1-cor(t(AMD))),method='%s')),Colv=as.dendrogram(hclust(dist(1-cor(AMD)),method='%s'))"%(clusterMethod,clusterMethod)
		cluster = "Rowv=TRUE,Colv=as.dendrogram(hclust(dist(1-cor(AMD)),method='%s'))"%(clusterMethod)
	elif cluster == 'abs.correlation':
		# cluster = "Rowv=as.dendrogram(hclust(dist(1-cor(t(AMD))),method='%s')),Colv=as.dendrogram(hclust(dist(1-cor(AMD)),method='%s'))"%(clusterMethod,clusterMethod)
		cluster = "Rowv=TRUE,Colv=as.dendrogram(hclust(dist(1-abs(cor(AMD))),method='%s'))"%(clusterMethod)
	elif cluster == 'sq.correlation':
		# cluster = "Rowv=as.dendrogram(hclust(dist(1-cor(t(AMD))),method='%s')),Colv=as.dendrogram(hclust(dist(1-cor(AMD)),method='%s'))"%(clusterMethod,clusterMethod)
		cluster = "Rowv=TRUE,Colv=as.dendrogram(hclust(dist(1-cor(AMD)^2),method='%s'))"%(clusterMethod)
	
	elif cluster == 'y.correlation':
		cluster = "Rowv=NA,Colv=as.dendrogram(hclust(dist(1-cor(AMD)),method='%s'))"%(clusterMethod)
	elif cluster == 'x.correlation':
		cluster = "Rowv=as.dendrogram(hclust(dist(1-cor(t(AMD))),method='%s')),Colv=NA"%(clusterMethod)
	elif cluster == 'y.abs.correlation':
		cluster = "Rowv=NA,Colv=as.dendrogram(hclust(dist(1-abs(cor(AMD))),method='%s'))"%(clusterMethod)
	elif cluster == 'x.abs.correlation':
		cluster = "Rowv=as.dendrogram(hclust(dist(1-abs(cor(t(AMD)))),method='%s')),Colv=NA"%(clusterMethod)
	elif cluster == 'y.sq.correlation':
		cluster = "Rowv=NA,Colv=as.dendrogram(hclust(dist(1-cor(AMD)^2),method='%s'))"%(clusterMethod)
	elif cluster == 'x.sq.correlation':
		cluster = "Rowv=as.dendrogram(hclust(dist(1-cor(t(AMD))^2),method='%s')),Colv=NA"%(clusterMethod)
	
	elif cluster == 'matrix.correlation':
		cluster = "Rowv=as.dendrogram(hclust(as.dist(1-AMD),method='%s')),Colv=as.dendrogram(hclust(as.dist(1-AMD),method='%s'))"%(clusterMethod,clusterMethod)
		invert=False
	elif cluster == 'matrix.euclidean':
		cluster = "Rowv=as.dendrogram(hclust(as.dist(AMD),method='%s')),Colv=as.dendrogram(hclust(as.dist(AMD),method='%s'))"%(clusterMethod,clusterMethod)
		invert=False
		
	# no clustering
	# cluster='Rowv=NA,Colv=NA'
	
	# redscale(N)
	# greenscale(N)
	# bluescale(N)
	# blueyellow(N)
	# redgreen(N)
	# jetColors(N)
	
	import sio as io
	# scale = "none", "row", "column"
	if not tmpfile: tmpfile = outdir+'Rheatmaptmp.txt'
	
	scriptname = outdir+scriptname
	# file = outdir+file
	
	# reverse order of mat
	if invert:
		mat2 = []; rows2 = []
		for i in range(len(mat)):
			mat2 = [mat[i]]+mat2
			if len(rows)==len(mat): rows2 = [rows[i]]+rows2
		mat = mat2
		rows = rows2
	
	if mat:
		mat2 = []
		for row in mat:
			newrow = []
			for x in row:
				if x == nastr: newrow += [0]
				else: newrow += [x]
			mat2 += [newrow]
		mat = mat2
		
		# mat = map(lambda x: map(lambda y: y==nastr and random.uniform(0,0.001) or y, x), mat)
		io.printTable(mat,cols,rows,file=tmpfile)
	
	if not file: return -1
	
	if len(cols):
		hasheader = 'TRUE'
	else:
		hasheader = 'FALSE'
	
	
	rowstr = 'NULL'
	if len(rows): rowstr = '1'
	
	colstr = 'FALSE'
	if len(cols): colstr = 'TRUE'
	
	extra = ''
	if len(title): extra += ', main="%s"'%(title)
	if len(xlabel): extra += ', xlab="%s"'%(xlabel)
	if len(ylabel): extra += ', ylab="%s"'%(ylabel)
	
	rc = col
	
	if len(cluster): cluster = ', '+cluster+', '
	
	keystr = 'TRUE'
	if key == False: keystr = 'FALSE'
	
	parts = file.split('/')
	format = parts[-1].split('.')[-1]
	# print 'FORMAT', format
	
	formatstr = '%s("%s")'%(format,file)
	if format == 'png':
		if width and height:
			formatstr = '%s("%s", width = %s, height = %s)'%(format,file, width, height)
		else:
			formatstr = '%s("%s")'%(format,file)
	# library(ClassDiscovery)
	Rscpt = """
	
	library(RColorBrewer)
	dat = read.delim("%s",header=%s)
	"""%(tmpfile,hasheader)
	if trimrow: Rscpt += 'dat = dat[,2:length(dat)]\n'
	Rscpt += 'AMD = as.matrix(dat)\n'
	Rscpt += "\n"
	Rscpt += custom
	Rscpt += "\n\n"
	
	if version != 'gplots':
		Rscpt += """
		
		# AMD = ut.cormat()
		
		%s
		heatmap(AMD, scale=%s, col=%s%s%s, cexRow=%s, cexCol=%s)
		dev.off()\n"""%(formatstr, scale, col, extra, cluster, rowfontsize, colfontsize)
	
	else:
		Rscpt += """library(gplots)
		library(ClassDiscovery)
		dat = read.delim("%s", row.names=%s, header=%s)
		%s
		heatmap.2(as.matrix(dat), Rowv=FALSE,Colv=FALSE, scale=%s, col=%s%s, cexRow=%s, cexCol=%s, key=%s, density.info="none", keysize=1, dendrogram="none")
		dev.off()\n"""%(tmpfile, rowstr, colstr, formatstr, scale, col, extra, rowfontsize, colfontsize, keystr)
	
	fh = open(scriptname, 'w')
	print >> fh, Rscpt
	fh.close()
	os.system("%s --vanilla < %s > %s"%(Rversion,scriptname, scriptname+'.out'))
	# os.system('rm -f "%s"'%(outdir+'rscpt.heatmap.R'))
	# os.system('rm -f "%s"'%(outdir+'rscpt.heatmap.R.out'))
	# if rmtmp:
		# os.system('rm -f "%s"'%(tmpfile))
	return 1



def pairsToMatrix(dat, bins=[], ntics=8, logscale='', dtype=r1, categorical=True, freq=False):
	# first make 2D histograms of the (X,Y) data
	
	if logscale: 
		Y = map(dtype, map(logfunc, vslice(dat,0)))
		X = map(dtype, map(logfunc, vslice(dat,1)))
	else:
		Y = map(dtype, vslice(dat,0))
		X = map(dtype, vslice(dat,1))
	
	zipXY = zip(X,Y)
	zipXY = pairedfilterstr(nastr)(zipXY)
	uXY = unique(zipXY)
	
	if not bins: bins=uXY
	# print 'bins', bins
	Hxy = histogram(zipXY, bins=bins, categorical=categorical, dtype=None)
	# print Hxy
	
	XY,Z = unzip(Hxy)
	X = sorted(vslice(XY,0))
	Y = sorted(vslice(XY,1))
	uX = unique(X); uY = unique(Y)
	
	# hash the XY pairs for fast lookup
	XYD = {}
	# print 'confirm XY', len(XY)
	for i in range(len(XY)): XYD[XY[i]] = i
	
	# aggregate results
	dct = {}
	for j in range(len(X)):
		for i in range(len(Y)):
			# try: dct[(X[j],Y[i])] += Z[XY.index((X[j],Y[i]))]
			try: dct[(X[j],Y[i])] += Z[XYD[(X[j],Y[i])]]
			except KeyError: dct[(X[j],Y[i])] = 0
			except ValueError: dct[(X[j],Y[i])] = 0
		
	Z2 = [[0 for j in range(len(uY))] for i in range(len(uX))]
	if not freq:
		for j in range(len(uX)):
			for i in range(len(uY)):
				Z2[j][i] = dct[(uX[j],uY[i])]
	else:
		denom = float(sum(dct.values()))
		for j in range(len(uX)):
			for i in range(len(uY)):
				Z2[j][i] = dct[(uX[j],uY[i])]/denom
	
	
	# if bins:
	# 	rngeystr = bins#makeRange(min(uX),max(uX),n=ntics)
	# 	rngexstr = bins#makeRange(min(uY),max(uY),n=ntics)
	# else:
	# 	rngey = makeRange(min(range(len(uX))),max(range(len(uX))),n=ntics)
	# 	rngex = makeRange(min(range(len(uY))),max(range(len(uY))),n=ntics)
	# 	
	# 	rngeystr = makeRange(min(uX),max(uX),n=ntics)
	# 	rngexstr = makeRange(min(uY),max(uY),n=ntics)
	
	return Z2#, rngex, rngey


def pylabHeatmap(dat, file, type='pairs', bins=[], freq=False, dtype=r1, xlabel='', ylabel='', stat=avg, ntics=8, shading='flat', logscale=10, cbarLabel='', ticrounding=r1, title=''):
	import pylab as pil, plot, numpy as np, sys
	
	Z2 = None
	
	logfunc = None
	if logscale>1: logfunc = logb(logscale)
	
	if type == 'pairs':
		# first make 2D histograms of the (X,Y) data
		if logscale: 
			Y = map(dtype, map(logfunc, vslice(dat,0)))
			X = map(dtype, map(logfunc, vslice(dat,1)))
		else:
			Y = map(dtype, vslice(dat,0))
			X = map(dtype, vslice(dat,1))
			
		zipXY = zip(X,Y)
		uXY = unique(zipXY)
		
		if not bins: bins=uXY
		Hxy = histogram(zipXY, bins=bins, categorical=True, dtype=None)
		
		XY,Z = unzip(Hxy)
		X = sorted(vslice(XY,0))
		Y = sorted(vslice(XY,1))
		uX = unique(X); uY = unique(Y)
		
		# aggregate results
		dct = {}
		for j in range(len(X)):
			for i in range(len(Y)):
				try: dct[(X[j],Y[i])] += Z[XY.index((X[j],Y[i]))]
				except KeyError: dct[(X[j],Y[i])] = 0
				except ValueError: dct[(X[j],Y[i])] = 0
		
		Z2 = [[0 for j in range(len(uY))] for i in range(len(uX))]
		if not freq:
			for j in range(len(uX)):
				for i in range(len(uY)):
					Z2[j][i] = dct[(uX[j],uY[i])]
		else:
			denom = float(sum(dct.values()))
			for j in range(len(uX)):
				for i in range(len(uY)):
					Z2[j][i] = dct[(uX[j],uY[i])]/denom
				
		rngey = makeRange(min(range(len(uX))),max(range(len(uX))),n=ntics)
		rngex = makeRange(min(range(len(uY))),max(range(len(uY))),n=ntics)
		
		rngeystr = makeRange(min(uX),max(uX),n=ntics)
		rngexstr = makeRange(min(uY),max(uY),n=ntics)
		
	elif type == 'triples':
		# we are provided with a list of X,Y,Z triples
		# convert this to a 2d matrix
		print
		
		# pile it up using a dictionary
		D = {}
		for i in range(len(dat)):
			y,x,z = dat[i]
			if not logscale:
				try:
					D[(dtype(x),dtype(y))] += [z]
				except KeyError:
					D[(dtype(x),dtype(y))] = [z]
			else:
				try:
					D[(dtype(logfunc(x)),dtype(logfunc(y)))] += [z]
				except KeyError:
					D[(dtype(logfunc(x)),dtype(logfunc(y)))] = [z]
				
				
		# get array dimensions
		xs = sorted(unique(vslice(D.keys(), 0)))
		ys = sorted(unique(vslice(D.keys(), 1)))
		Z2 = [[0 for j in range(len(ys))] for i in range(len(xs))]
		# fill in
		for x,y in D.keys():
			ix = xs.index(x)
			iy = ys.index(y)
			Z2[ix][iy] = stat(D[(x,y)])
		
		rngey = map(int,makeRange(min(range(len(xs))),max(range(len(xs))),n=ntics))
		rngex = map(int,makeRange(min(range(len(ys))),max(range(len(ys))),n=ntics))
		
		rngeystr = makeRange(min(xs),max(xs),n=ntics)
		rngexstr = makeRange(min(ys),max(ys),n=ntics)
		
	elif type == 'matrix':
		# we are provided with a 2D matrix
		print
	
	if not cbarLabel:
		cbarstring = 'Count'
		if freq: cbarstring = 'Frequency'
	else: cbarstring = cbarLabel
		
	pil.pcolormesh(np.array(Z2), shading=shading)
	
	if logscale:
		rngexstr = map(lambda x: logscale**x, rngexstr)
		rngeystr = map(lambda x: logscale**x, rngeystr)
	
	pil.xticks(rngex, map(ticrounding,rngexstr))
	pil.yticks(rngey, map(ticrounding,rngeystr))
	pil.ylabel(ylabel)
	pil.xlabel(xlabel)
	if title: pil.title(title)
	cbar = pil.colorbar(); cbar.set_label(cbarstring)
	pil.savefig(file)
	pil.clf()
	return 1


def chart(dat, title='', xlabel='', ylabel='', file='', color=1, xlabels=[], boxwidth=1, error=[], errorwidth=.05, legend=None, custom='', fontface=None, fontsize=None, rotateLabels=-45, style='filled', bottombar=False, basey=None,sumstat=avg, errorstat=stderr, logscale=False, logscaleymin=0.01, ratio=1):
	try: import Gnuplot as gp
	except ImportError: return -1
	g = gp.Gnuplot()
	
	if len(dat) > 0 and type(dat[0]) == type([]):
		if not len(error): error = map(errorstat,dat)
		dat = map(sumstat,dat)
    
	if getUnderscore()==False: xlabels = [str(l).replace('_','-') for l in xlabels]
	if title and getUnderscore()==False: title = title.replace('_','-')
	
	if not len(dat): 
		print >> sys.stderr, 'Error: no input data provided.'
		return -1
	
	format = GLOB_FORMAT
	face = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	if file: fsize /= 2.0
	
	dat = map(lambda x: x!= nastr and x or 0, dat)
	error = map(lambda x: x!= nastr and x or 0, error)
	
	g('set terminal pdf enhanced font "'+face+','+str(fsize)+'"')
	
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file)
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - \"%s\""'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
	g.ylabel(ylabel)
	g.xlabel(xlabel)
	if title: g.title(title)
	
	g('set font "%s,%s"'%(face,fsize))
	g('set style data boxes')
	if style == 'filled':
		g('set style fill solid border -1')
	
	g('set boxwidth %s'%(boxwidth))
	
	lfd = len(dat)
	
	# if you provide dat = [1,2,3] make this into dat = [(0,1), (1,2), (2,3)]
	if lfd and type(dat[0]) != type([]) and type(dat[0]) != type(()):
		dat = map(lambda x: (x,dat[x]), range(lfd))
	
	
	# scale font size if need be
	fticsize = minna([fsize, divna(fsize*3.3,ln(float(lfd)))])
	
	g('set ytics font "%s,%s"'%(face,fticsize))
	g('set ylabel font "%s,%s"'%(face,fticsize))
	
	if xlabels:
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, range(len(xlabels))))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		
		g('set xtics nomirror rotate by %s '%(rotateLabels)+xticstr+' font "%s,%s"'%(face, fticsize))
	
	#
	if logscale == 'x' or logscale == True: g('set logscale x')
	elif logscale == 'y': g('set logscale y')
	elif logscale == 'xy' or logscale: g('set logscale')
	
	for i in range(len(dat)):
		if dat[i][1] == 0 and logscale != False and 'y' in logscale:
			dat[i] = (dat[i][0], dat[i][1]+logscaleymin)
	
	
	# handle error bars
	# assume err is a list of ([low,high] pairs) or (low) values
	top = 0
	bot = GLOB_INF
	use = 0
	if len(error):
		err = error
		for i in range(len(err)):
			er = err[i]
			if er == nastr or er == 0: er = [0,0]
			elif type(er) == type(1.): er = [er,er]
			else: er = list(er)
			# else: print 'ISSUE', er
			err[i] = er
		
		# draw them as arrows with no heads at each x-value
		xs = vslice(dat,0)
		ys = vslice(dat,1)
		
		# compute average delta
		if len(xs) > 1:
			dx = 4*avg(map(lambda x: abs(xs[x]-xs[x-1]), range(1,len(xs))))
		else:
			dx = 2
			
		
		
		for zi in range(len(xs)):
			xlow = xs[zi] - dx*float(errorwidth)
			xhi = xs[zi] + dx*float(errorwidth)
			
			ex = ys[zi]+err[zi][1]
			if ex > top: top = ex
			ex = ys[zi]-err[zi][1]
			if bot == None or ex < bot: bot = ex
			
			# and draw the verical line
			g('set arrow nohead from %s,%s to %s,%s front'%(xs[zi], ys[zi], xs[zi], ys[zi]+err[zi][1]))
			if bottombar: g('set arrow nohead from %s,%s to %s,%s front'%(xs[zi], ys[zi], xs[zi], ys[zi]-err[zi][1]))
			# draw the bars
			# g('set arrow nohead from %s,%s to %s,%s as 1'%(xlow, ys[zi], xhi, ys[zi]))
			g('set arrow nohead from %s,%s to %s,%s front'%(xlow, ys[zi]+err[zi][1], xhi, ys[zi]+err[zi][1]))
			if bottombar: g('set arrow nohead from %s,%s to %s,%s front'%(xlow, ys[zi]-err[zi][1], xhi, ys[zi]-err[zi][1]))
			
		# account for error bars
		use = ys[0]
		if len(ys)>1:
			use = max(map(lambda x: abs(ys[x]-ys[x-1]), range(1,len(ys))))
			# print 'use', use, top
		
	
	# print 'ok', basey, top, use, .05*use
	if basey != None:
		g('set yrange [%s:%s]'%(basey,top + .05*use))
	elif bottombar:
		# print 'test', bot,.05*use
		g('set yrange [%s:%s]'%(bot-.05*use, top+.05*use))
	else:
		g('set yrange [*:%s]'%(top + .05*use))
	
	
	
	default = 'set size ratio %s; '%(ratio)
	g(default+custom)
	
	if legend: g.plot(gp.Data(dat,with_='boxes lt %d' % (color), title=legend))
	else: g.plot(gp.Data(dat,with_='boxes lt %d' % (color)))
	g('set out')


def stackedChart(dat, header, errordat=[], error=[], title='', xlabel='', ylabel='y-axis', file='', color=-1, tmpfile='tmp.dat', outdir='', key=True, custom='', rmtmp=False, useonly=0, xlabels=[], rounding=None, rotateLabels=-45, fontface='Helvetica', fontsize=0, verbose=0, errorwidth=.05, bottombar=False, grid='', ratio=1):
	
	try: import Gnuplot as gp
	except ImportError: return -1
	g = gp.Gnuplot()
	
	face = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	if file: fsize /= 3.0
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	
	if getUnderscore()==False: xlabels = [l.replace('_','-') for l in xlabels]
	
	g('set terminal pdf enhanced font "'+face+','+str(fsize)+'"')
	format = GLOB_FORMAT
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file)
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - \"%s\""'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
	
	g.ylabel(ylabel)
	g.xlabel(xlabel)
	if title: g.title(title)
	
	dat = map(list, dat)
	
	header = ['Key']+header
	if xlabels:
		# need to append to beginning of each row
		if not header: header = range(len(dat[0]))
		if not rounding:
			for i in range(len(dat)): dat[i] = [xlabels[i]]+dat[i]
		else:
			for i in range(len(dat)): dat[i] = [xlabels[i]]+map(rounding, dat[i])
	else:
		if not rounding:
			for i in range(len(dat)): dat[i] = ['%s'%(i+1)]+dat[i]
		else:
			for i in range(len(dat)): dat[i] = ['%s'%(i+1)]+map(rounding, dat[i])
	
	dorotate = ''
	if rotateLabels!=None: dorotate = 'rotate by %s'%(rotateLabels)
	
	if not outdir: outdir = os.getcwd()+'/'
	tmpfile = outdir+tmpfile
	import sio as io
	if verbose:
		print 'Printing stacked chart to', tmpfile
		io.printFormattedTable(dat,header=header)
	# print 'dat', dat
	io.printFormattedTable(dat,header=header,file=tmpfile)
	
	g("""set boxwidth 0.75 absolute
	set style fill solid 1.00 border -1
	set style histogram rowstacked
	set style data histograms""")
	g('unset key')
	
	
	if grid != '':
		if grid == 'xy': g('set grid')
		else: g('set grid %s'%(grid))
	
	
	if key:
		g("set key outside right top vertical Left reverse enhanced autotitles columnhead nobox")
		g("set key invert samplen 4 spacing 1 width 0 height 0")
	
	xlabelspacing = None
	if xlabels: 
		userlabels = True
		
		# get xvalues for the data
		xlabelrange = []
		if not xlabelspacing:
			xlabelrange = map(str, xlabels)
			# xlabels2 = map(lambda x: round(x,1), xlabels)
			# xlabelrange = map(str, map(lambda x: x in xlabels2 and x or ' ', bins))
		else:
			xlabelrange = map(str, xlabels)
		
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), xlabelrange)
		# tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, xlabels))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, xticstr, fontface, fsize))
	
	
	useonly = len(dat[0])-1
	# if useonly == 0: useonly = 3+len(dat[0])-1
	# stng = ''.join([", '' using %d" % (i) for i in range(3,3+len(dat[0])-1)])
	stng = ''.join([", '' using %d" % (i) for i in range(3,useonly+2)])
	
	# handle error bars
	# assume err is a list of ([low,high] pairs) or (low) values
	if len(error):
		err = error
		for i in range(len(err)):
			er = err[i]
			if er == nastr or er == 0: er = [0,0]
			elif type(er) == type(1.): er = [er,er]
			else: er = list(er)
			err[i] = er
		
		# draw them as arrows with no heads at each x-value
		xs = range(len(errordat))
		ys = errordat
		
		# compute average delta
		if len(xs) > 1:
			dx = 4*avg(map(lambda x: abs(xs[x]-xs[x-1]), range(1,len(xs))))
		else:
			dx = 2
			
		top = None
		
		for zi in range(len(xs)):
			xlow = xs[zi] - dx*float(errorwidth)
			xhi = xs[zi] + dx*float(errorwidth)
			
			ex = ys[zi]+err[zi][1]
			if ex > top: top = ex
			
			# and draw the verical line
			g('set arrow nohead from %s,%s to %s,%s front'%(xs[zi], ys[zi], xs[zi], ys[zi]+err[zi][1]))
			if bottombar: g('set arrow nohead from %s,%s to %s,%s front'%(xs[zi], ys[zi], xs[zi], ys[zi]-err[zi][1]))
			# draw the bars
			# g('set arrow nohead from %s,%s to %s,%s as 1'%(xlow, ys[zi], xhi, ys[zi]))
			g('set arrow nohead from %s,%s to %s,%s front'%(xlow, ys[zi]+err[zi][1], xhi, ys[zi]+err[zi][1]))
			if bottombar: g('set arrow nohead from %s,%s to %s,%s front'%(xlow, ys[zi]-err[zi][1], xhi, ys[zi]-err[zi][1]))
			
		# account for error bars
		use = ys[0]
		if len(ys)>1:
			use = max(map(lambda x: abs(ys[x]-ys[x-1]), range(1,len(ys))))
		g('set yrange [*:%s]'%(top + .05*use))
	
	
	default = 'set size ratio %s; '%(ratio)
	g(default+custom)
	
	g("plot '"+tmpfile+"' using 2:xtic(1)"+stng)
	g('set out')
	
	# if tmpfile != 'tmp.dat': rmtmp = False
	if rmtmp: os.system('rm -f '+tmpfile)


def box(dat, style='group', xvals=[], xlabel=' ', ylabel='Value', file='', logscale=False, plot=True, title='', key='', whiskers='full', tails=[.025,0.975], showtails=True, margin=.05, inkey='', xlabels=[], statistic='', logtype='gnuplot', apre=False, custom='', boxwidth=.75, rotateLabels=-45, alpha=0.05, fontface=None, fontsize=None, ratio=1):
	"""
	Takes a list of [x,y] value pairs and gives a set of box plots at each x-value.
	The style argument specifies if boxplots should be generated by a list of data points
	at each index of dat (group), or as lists of [x,y] dat pairs at specified xvalues.
	
	"""
	
	if getUnderscore()==False: xlabels = [str(l).replace('_','-') for l in xlabels]
	if title and getUnderscore()==False: title = title.replace('_','-')
	
	dat = castemat(dat,floatna)
	
	def logtmp(x):
	 	if x==nastr or float(x)<=0.0: return 1e-18
		try: return math.log(float(x),10)
		except OverflowError: return 1e-18
	
	face = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	if file: fsize /= 2.0
	
	
	try: import Gnuplot as gp
	except ImportError: return -1
	
	xlabels = map(str,xlabels)
	if getUnderscore()==False: xlabels = [l.replace('_','-') for l in xlabels]
	
	tfgroups = {}
	if style=='group': 
		assert (len(xvals)==0), "xvals should not be specified for this style."
		
		xvals = range(len(dat))
		tfgroups = dat
		if not xlabels and len(dat) > 1:
			xlabels = xvals[:]
		if xvals == []: xvals = [0]
	else:
		assert (len(xvals)), \
			"xvals list must be specified for this type of plot."
		
		# group gene variances by number of tfs
		for tfn in range(int(maxna(xvals)+1)): tfgroups[tfn] = []
		for tfn,gvar in dat: tfgroups[int(tfn)] += [gvar]
		
	yvals = []; allvarvals = []; pairs = []
	boxes = []; medboxes = []
	varz = []; medvarz = []
	lowvals = []; highvals = []
	# print 'xvals',xvals
	supermax = '*'
	for tfn in range(int(maxna(xvals)+1)):
		# if logscale and logtype != 'gnuplot': 
		# 	gvars = map(log, tfgroups[tfn])
		# else:
		
		gvars = filterstr(nastr)(tfgroups[tfn])
		
		sm = maxna(gvars)
		if supermax == '*' or (sm > supermax and sm != nastr): supermax = sm
		
		# print xlabels[tfn], 'gvars', sorted(gvars)
		if logscale and logtype != 'gnuplot': 
			p00 = logtmp(percentile(gvars, tails[0]))
			p25 = logtmp(percentile(gvars, .25))
			p50 = logtmp(percentile(gvars, .5))
			p75 = logtmp(percentile(gvars, .75))
			p10 = logtmp(percentile(gvars, tails[1]))
			
			# get extreme values
			plrem = [logtmp(zz) for zz in gvars if logtmp(zz) < p00]
			phrem = [logtmp(zz) for zz in gvars if logtmp(zz) > p10]
			
		else:
			p00 = percentile(gvars, tails[0])
			p25 = percentile(gvars, .25)
			p50 = percentile(gvars, .5)
			p75 = percentile(gvars, .75)
			p10 = percentile(gvars, tails[1])
			
			# get extreme values
			plrem = [zz for zz in gvars if zz < p00]
			phrem = [zz for zz in gvars if zz > p10]
			if len(plrem) == 0: plrem = [0]
			if len(phrem) == 0: phrem = [0]
			
			# print 'cutoff', p00, 'rem', plrem
			# print 'cutoff', p10, 'rem', phrem
			# print
			
			
			
		yvals += [p10,p00]
		# print xlabels[tfn], p25, p00, p10, p75
		boxes += [ [tfn, p25, p00, p10, p75] ]
		medboxes += [ [tfn, p50, p50, p50, p50] ]
		if showtails:
			lowvals += map(lambda x: (tfn,x), plrem)
			highvals += map(lambda x: (tfn,x), phrem)
		
		stat = nastr; L1 = nastr; L2 = nastr
		# if not statistic: whiskers = 'mean'
		if statistic != '': whiskers = 'statistic'
		
		if statistic == 'mean': stat, L1, L2 = meanWithConfidence(gvars, alpha=alpha)
		elif statistic == 'variance': stat, L1, L2 = varianceWithConfidence(gvars, alpha=alpha, method='equal tails')
		elif len(statistic): raise TypeError('statistic must be one of {\'mean\', \'variance\'}')
		
		allvarvals += [stat,L1,L2]; pairs += [[tfn,stat]]
		if L1 != nastr and L2 != nastr and stat != nastr:
			varz += [ [tfn, stat, L1, L2, stat] ]
		if stat != nastr:
			medvarz += [ [tfn, stat, stat, stat, stat] ]
		
	if inkey:
		ret = corWithPvalue(pairs, method=inkey, tails='2', p='approx')
		if ret['P'] == nastr: ret['P'] = 1
		key = '{/Helvetica-Oblique s} = %.2f, {/Helvetica-Oblique P} < %.2e'%(ret['r'],ret['P'])
		# print 'boxplot!', ret['r'], ret['p-value']
		# print 'varvals', pairs
	
	# plot the results
	g = gp.Gnuplot()
	format = GLOB_FORMAT
	g('set terminal pdf enhanced font "'+face+','+str(fsize)+'"')
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file)
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - \"%s\""'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
		
	
	if not plot:
		g('set terminal %s enhanced' % (format))
		g('set out')
	
	
	if logscale and logtype == 'gnuplot': g('set logscale y')
	if title: g('set title "%s"' % (title))
	g('set xlabel "%s"' % (xlabel)); g('set ylabel "%s"' % (ylabel))
	
	mnx = minna(xvals); mxx = maxna(xvals); xrnge = listRange(xvals)
	mny = minna(yvals); mxy = maxna(yvals); yrnge = listRange(yvals)
	if whiskers != 'full':
		mny = minna(allvarvals); mxy = maxna(allvarvals); yrnge = listRange(allvarvals)
	
	# print 'min/max', mny,mxy,yrnge
	
	# if mnx == ut.nastr: mnx = '*'
	dx = max(margin*xrnge,0); dy = max(margin*yrnge,0)
	if mnx == nastr: 
		mnx = 0
		#xlow = '*'
	xlow = mnx-dx
	
	# if mny == nastr: ylow = '*'
	# else: ylow = mny-dy
	# print 'mxx', mxx, 'dx', dx, 'mxy', mxy, 'dy', dy
	xhi = mxx+dx
	# yhi = mxy+dy
	
	g('set boxwidth %s'%(boxwidth))
	
	# if logscale:
	# 	yhi = max(xhi,1e-6)
	# 	# ylow = max(ylow,1e-6)
	
	rng = listRange([xlow,xhi])
	delt = boxwidth/2.#*rng
	# print 'test', xlow, boxwidth*rng
	g('set xrange [%s:%s]'%(xlow-delt,xhi+delt))
	
	if len(xvals) == 1: g('set xrange [%f:%f]'%(-.5,.5))
	
	# print 'supermax', supermax
	try: g('set yrange [*:%f]'%(supermax+.025*supermax))
	except TypeError: g('set yrange [*:*]')
	
	dorotate = ''
	if rotateLabels != None:  dorotate = 'rotate by %s'%(rotateLabels)
	
	fsize = GLOB_FONT_SIZE
	if file: fsize /= 3.0
	
	if xlabels:
		tmp = zip(map(lambda x: '"'+str(x)+'"', xlabels), map(str, range(len(xlabels))))
		# tmp is [ ["lab1",0], ["lab2",1],... ]
		tmp2 = map(lambda x: ' '.join(x), tmp)
		xticstr = '('+','.join(map(str, tmp2))+')'
		
		g('set xtics nomirror %s '%(dorotate)+xticstr+' font "%s,%s"'%(GLOB_FONT_FACE, fsize))
	else:
		g('set xtics nomirror %s %s font "%s,%s"'%(dorotate, '', GLOB_FONT_FACE, fsize))
	
	default = 'set size ratio %s; '%(ratio)
	g(default+custom)
	
	# print 'tails:', lowvals, highvals
	
	if whiskers == 'full':
		if showtails:
			if key: 
				g.plot(gp.Data(boxes, with_='candlesticks lt -1 lw 1.5 whiskerbars'),
				       gp.Data(boxes, with_='candlesticks lt 0 lw 1'),
				       gp.Data(medboxes, with_='candlesticks lt 1 lw 2.5', title=key),
				       gp.Data(lowvals, with_='points lt 1 ps 1 pt 1'),
				       gp.Data(highvals, with_='points lt 1 ps 1 pt 1'),
				)
			else:
				g.plot(gp.Data(boxes, with_='candlesticks lt -1 lw 1.5 whiskerbars'),
				       gp.Data(boxes, with_='candlesticks lt 0 lw 1'),
				       gp.Data(medboxes, with_='candlesticks lt 1 lw 2.5'),
				       gp.Data(lowvals, with_='points lt 1 pt 5 ps .5'),
				       gp.Data(highvals, with_='points lt 1 pt 5 ps .5'),
				)
		else:
			if key: 
				g.plot(gp.Data(boxes, with_='candlesticks lt -1 lw 1.5 whiskerbars'),
				       gp.Data(boxes, with_='candlesticks lt 0 lw 1'),
				       gp.Data(medboxes, with_='candlesticks lt 1 lw 2.5', title=key)
				)
			else:
				g.plot(gp.Data(boxes, with_='candlesticks lt -1 lw 1.5 whiskerbars'),
				       gp.Data(boxes, with_='candlesticks lt 0 lw 1'),
				       gp.Data(medboxes, with_='candlesticks lt 1 lw 2.5')
				)
			
	else:
		if key: 
			g.plot(gp.Data(varz, with_='candlesticks lt -1 lw 2 whiskerbars'),
			       gp.Data(varz, with_='candlesticks lt 0 lw 0'),
			       gp.Data(medvarz, with_='points pt "." lt 1 lw 8',title=key)
			)
		else:
			g.plot(gp.Data(varz, with_='candlesticks lt -1 lw 1 whiskerbars'),
			       gp.Data(varz, with_='candlesticks lt 0 lw 1'),
			       gp.Data(medvarz, with_='points pt "." lt 1 lw 8')
			)
	if file and apre != True: g('set out')
	return g













def multiStackedChart(dat, header, title='', xlabel='x-axis', ylabel='y-axis', file='', colors=[], tmpfile='tmp.dat', outdir='', key=True, custom='', rmtmp=False, useonly=0, xlabels=[], rounding=None, rotateLabels=-45, fontface='Helvetica', fontsize=0, boxwidth=1, spacing=.9, legends=[]):
	
	try: import Gnuplot as gp
	except ImportError: return -1
	g = gp.Gnuplot()
	
	face = GLOB_FONT_FACE
	fsize = GLOB_FONT_SIZE
	if file: fsize /= 3.0
	if fontface: face = fontface
	if fontsize: fsize = fontsize
	
	g('set terminal pdf enhanced font "'+face+','+str(fsize)+'"')
	format = GLOB_FORMAT
	if file:
		crap,format,crap = re.split('^.+\.(.+)$', file)
		if format == 'ps': format = 'postscript'
		
		if togglepdf:
			g('set terminal postscript enhanced')
			g('set out "| ps2pdf - \"%s\""'%(file))
		else:
			g('set terminal %s enhanced' % (format))
			g('set out "'+file+'"')
		
	
	g.ylabel(ylabel)
	g.xlabel(xlabel)
	if title: g.title(title)
	
	
	metadat = map(lambda d: map(list, d), dat)
	
	if not header: header = range(len(metadat[0][0]))
	header = ['Key']+header
	if not xlabels and not legends:
		xlabels = ['%s'%(i+1) for i in range(len(metadat[0]))]
	elif not xlabels:
		xlabels = ['' for i in range(len(metadat[0]))]
	
	if not legends:
		legends = [str(i) for i in range(len(metadat))]
	
	
	
	# interleave the data
	ndata = len(metadat)
	dat = []
	
	blankdata = ['0' for i in range(len(metadat[0][0]))]
	
	# need to append to beginning of each row
	if not rounding:
		# for each sample
		for i in range(len(metadat[0])):
			# for each data set
			for zz in range(ndata):
				astr = '%s%s'%(legends[zz], xlabels[i])
				dat += [ [astr]+metadat[zz][i] ]
			# add a blank between
			dat += [ ['.']+blankdata ]
	else:
		for i in range(len(metadat[0])):
			for zz in range(ndata):
				astr = '%s%s'%(legends[zz], xlabels[i])
				dat += [ [astr]+map(rounding, metadat[zz][i]) ]
			# add a blank between
			dat += [ ['.']+blankdata ]
	
	# remove final blank
	dat = dat[:-1]
	
	if not outdir: outdir = os.getcwd()+'/'
	tmpfile = outdir+tmpfile
	# print 'printing stacked chart to', tmpfile
	import sio
	# io.printFormattedTable(dat,header=header)
	io.printFormattedTable(dat,header=header,file=tmpfile)
	
	g('set boxwidth %s'%(boxwidth))
	g("""set style fill pattern 1 border -1
	set style histogram rowstacked
	set style data histograms""")
	g('unset key')
	
	dorotate = ''
	if rotateLabels: dorotate = 'rotate by %s'%(rotateLabels)
	g('set xtics nomirror %s font "%s,%s"'%(dorotate, face, fsize))#*(3/float(ndata))))
	
	
	if custom: g(custom)
	
	if key:
		g("set key outside right top vertical Left reverse enhanced autotitles columnhead nobox")
		g("set key invert samplen 4 spacing 1 width 0 height 0")
	
	
	useonly = len(dat[0])-1
	# if useonly == 0: useonly = 3+len(dat[0])-1
	stng = ''.join([", '' using %d" % (i) for i in range(3,useonly+2)])
	g("plot '"+tmpfile+"' using 2:xtic(1)"+stng)
	g('set out')
	
	# if tmpfile != 'tmp.dat': rmtmp = False
	if rmtmp: os.system('rm -f '+tmpfile)

