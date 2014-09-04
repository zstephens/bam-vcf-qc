#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

python genReport.py bamMetrics.data allErrors outReportDir

************************************************** """

import sys
import os
import re
import time
import numpy as np
import matplotlib.pyplot as mpl

# sort list of strings according a number at the end of element
def sortByNumericSuffix(l):
	ol = sorted(l)
	regions  = []
	sortInds = []
	startRI = 0
	for i in xrange(len(ol)):
		ind = len(ol[i])-1
		while ol[i][ind].isdigit() and ind >= 0:
			ind -= 1
		ind += 1
		if i == 0:
			prevRP = ol[0][:ind]
		if ol[i][:ind] != prevRP:
			regions.append((startRI,i))
			startRI = i
			prevRP = ol[i][:ind]
		if not ol[i][-1].isdigit():
			sortInds.append((0,i))
		else:
			sortInds.append((int(ol[i][ind:]),i))
	regions.append((startRI,len(ol)))

	ool = []
	for r in regions:
		si = sorted(sortInds[r[0]:r[1]])
		for n in si:
			ool.append(ol[n[1]])
	return ool

def mpl_xticks_ints(eps=0.00001):
	locs, labels = mpl.xticks()
	xlab = []
	for n in locs:
		if abs(n-int(n)) < eps:
			xlab.append(str(int(n)))
		else:
			xlab.append('')
	mpl.xticks(locs,xlab)

MAX_READLEN = 10000

def main():

	tt = time.time()

	if len(sys.argv) != 4:
		print '\npython genReport.py bamMetrics.data allErrors outReport.html\n'
		exit(1)
	else:
		inData = sys.argv[1]
		inErrs = sys.argv[2]
		outRep = sys.argv[3]
		if outRep[-1] != '/':
			outRep += '/'

	sys.stdout.write('generating output report...')
	sys.stdout.flush()

	# eventually just hardcode all the reportHeader stuff into a single string, will make my life easier.
	inf = open('reportHeader.html','r')
	REPORT_HEADER = inf.read()
	inf.close()

	"""*************************************************"""
	"""		PARSE DATA FILE PRODUCED BY BAMQC APP		"""
	"""*************************************************"""

	myAligner = None
	myVersion = None

	nRecords        = None
	nUniquelyMapped = None
	nMultiMappedIDs = None
	nUnmapped       = None
	totalUniqueIDs  = None

	refLen_byRef        = {}
	readCount_byRef     = {}
	baseCount_byRef     = {}
	nReads_byReadLength = {}
	myNFreq             = {}
	myIMD               = {}
	covHist_byRef       = {}
	mapQHist_byRef      = {}

	f = open(inData,'r')
	for line in f:
		splt = line[:-1].split('\t')

		if splt[0][:2] == '>>':
			read_rf = False
			read_rc = False
			read_rl = False
			read_nf = False
			read_id = False
			read_mh = False
			read_ch = False

		if (read_mh or read_ch) and splt[0][0] == '>':
			currentRefName = splt[0][1:]
			if read_mh:
				mapQHist_byRef[currentRefName] = {}
			elif read_ch:
				covHist_byRef[currentRefName] = {}

		if splt[0] == '>>ALIGNER':
			myAligner = splt[1]
			myVersion = splt[2]
		elif splt[0] == '>>RECORDS_READ':
			nRecords = int(splt[1])
		elif splt[0] == '>>UNIQUELY_MAPPED':
			nUniquelyMapped = int(splt[1])
		elif splt[0] == '>>MULTI_MAPPED_IDS':
			nMultiMappedIDs = int(splt[1])
		elif splt[0] == '>>UNMAPPED':
			nUnmapped = int(splt[1])
		elif splt[0] == '>>TOTAL_UNIQUE_IDS':
			totalUniqueIDs = int(splt[1])

		elif splt[0] == '>>REF_LENGTHS':
			read_rf = True
		elif splt[0] == '>>READ_COUNT_PER_REF':
			read_rc = True
		elif splt[0] == '>>READ_LENGTH':
			read_rl = True
		elif splt[0] == '>>N_FREQUENCY':
			read_nf = True
		elif splt[0] == '>>IMD_HISTOGRAM':
			read_id = True
		elif splt[0] == '>>MAPQ_HISTOGRAM':
			read_mh = True
		elif splt[0] == '>>COVERAGE_HISTOGRAM':
			read_ch = True

		else:
			if read_rf:
				refLen_byRef[splt[0]] = int(splt[1])
			elif read_rc:
				readCount_byRef[splt[0]] = int(splt[1])
				baseCount_byRef[splt[0]] = int(splt[2])
			elif read_rl:
				nReads_byReadLength[int(splt[0])] = int(splt[1])
			elif read_nf:
				myNFreq[int(splt[0])] = int(splt[1])
			elif read_id:
				myIMD[int(splt[0])] = int(splt[1])
			elif read_mh and splt[0][0] != '>':
				mapQHist_byRef[currentRefName][int(splt[0])] = int(splt[1])
			elif read_ch and splt[0][0] != '>':
				covHist_byRef[currentRefName][int(splt[0])] = int(splt[1])
	f.close()

	allErrorData = []
	f = open(inErrs,'r')
	for line in f:
		allErrorData.append(line[:-1].split('\t'))
	f.close()

	refList = sortByNumericSuffix(readCount_byRef.keys())


	"""*********************************"""
	"""		CREATE OUTPUT REPORT		"""
	"""*********************************"""


	# will overwrite existing report if it exists
	if not os.path.isdir(outRep):
		os.system('mkdir '+outRep)
	if not os.path.isdir(outRep+'images/'):
		os.system('mkdir '+outRep+'images/')
	if not os.path.isdir(outRep+'logs/'):
		os.system('mkdir '+outRep+'logs/')

	of = open(outRep+'bamQC_report.html','w')

	REPORT_HEADER = re.sub("<title>-","<title>"+"BamQC Report",REPORT_HEADER)
	of.write(REPORT_HEADER)

	OUT_HTML = '<h2><font color="#444444">'+'BamQC Report'+'</font></h2>\n\n'

	#
	#	basic stats table
	#
	OUT_HTML += '<h2>Basic Stats:</h2>\n'
	OUT_HTML += '<div class="datagrid">\n'
	OUT_HTML += '<table>\n'
	OUT_HTML += '<tbody>\n'
	OUT_HTML += '<tr class="alt"><td><b>'+'Aligner:'+'</b></td><td>'+myAligner+' '+myVersion+'</td></tr>\n'
	OUT_HTML += '<tr><td><b>'+'# BAM Records:'+'</b></td><td>'+'{0:,}'.format(nRecords)+'</td></tr>\n'
	OUT_HTML += '<tr class="alt"><td><b>'+'# Total Unique Reads:'+'</b></td><td>'+'{0:,}'.format(totalUniqueIDs)+'</td></tr>\n'
	OUT_HTML += '<tr><td><b>'+'# Uniquely-Mapped Reads:'+'</b></td><td>'+'{0:,}'.format(nUniquelyMapped)+'</td></tr>\n'
	OUT_HTML += '<tr class="alt"><td><b>'+'# Multi-Mapped Reads:'+'</b></td><td>'+'{0:,}'.format(nMultiMappedIDs)+'</td></tr>\n'
	OUT_HTML += '<tr><td><b>'+'# Unmapped Reads:'+'</b></td><td>'+'{0:,}'.format(nUnmapped)+'</td></tr>\n'
	OUT_HTML += '</tbody>\n'
	OUT_HTML += '</table>\n'
	OUT_HTML += '</div>\n\n'

	#
	#	warnings table
	#
	OUT_HTML += '<h2>Quality Control Warnings/Failures:</h2>\n'

	if len(allErrorData) > 0:
		OUT_HTML += '<div class="datagrid">\n'
		OUT_HTML += '<table>\n'
		OUT_HTML += '<thead><tr><th>Program</th><th>Error</th><th>Severity</th><th># Records</th></tr></thead>\n'
		OUT_HTML += '<tbody>\n'

		isEven = False
		for n in allErrorData:
			if isEven:
				OUT_HTML += '<tr>'
			else:
				OUT_HTML += '<tr class="alt">'
			isEven = not(isEven)
			OUT_HTML += '<td>'+n[0]+'</td>\n'
			OUT_HTML += '<td>'+n[1]+'</td>\n'
			OUT_HTML += '<td><span class="'+n[2]+'Text">'+n[2]+'</span></td>\n'
			OUT_HTML += '<td>'+n[3]+'</td>\n'

		OUT_HTML += '</tbody>\n'
		OUT_HTML += '</table>\n'
		OUT_HTML += '</div>\n<br>\n\n'
	else:
		OUT_HTML += 'No errors.<br>\n'

	#
	#	read-length distribution
	#
	OUT_HTML += '<h2>Read Statistics:</h2>\n'
	OUT_HTML += '<a href="images/readLength_distribution.png"><img width="300" height="200" src="images/readLength_distribution.png" /></a>\n'

	#
	#	overall N-freq distribution
	#
	OUT_HTML += '<a href="images/nFreq_distribution.png"><img width="300" height="200" src="images/nFreq_distribution.png" /></a>\n'
	OUT_HTML += '<br>\n\n'

	#
	#	per-ref breakdown
	#
	COV_THRESH = 1

	OUT_HTML += '<h2>Mapping Statistics:</h2>\n'
	OUT_HTML += '<div class="datagrid">\n'
	OUT_HTML += '<table>\n'

	OUT_HTML += '<thead><tr><th>Reference</th>'
	OUT_HTML += '<th>Length (bp)</th>'
	OUT_HTML += '<th># Mapped (% Total)</th>'
	OUT_HTML += '<th>Bases >= '+str(COV_THRESH)+'x</th>'
	OUT_HTML += '<th>Avg Coverage</th>'
	OUT_HTML += '<th>Additional Info</th>'
	OUT_HTML += '</tr></thead>\n'

	OUT_HTML += '<tbody>\n'
	isEven = False
	refInd = 0
	for k in refList:
		if isEven:
			OUT_HTML += '<tr>'
		else:
			OUT_HTML += '<tr class="alt">'
		isEven = not(isEven)
		OUT_HTML += '<td><b>'+k+'</b></td>'
		OUT_HTML += '<td align="center">'+'{0:,}'.format(refLen_byRef[k])+'</td>'
		OUT_HTML += '<td align="center">'+'{0:,}'.format(readCount_byRef[k])+' <b>({0:.2f}%)</b>'.format(100.*readCount_byRef[k]/float(nRecords - nUnmapped))+'</td>'
		basesCov = sum([covHist_byRef[k][m] for m in covHist_byRef[k].keys() if m >= COV_THRESH])
		OUT_HTML += '<td align="center">'+'<b>{0:.2f}%</b>'.format(100.*basesCov/float(refLen_byRef[k]))+'</td>'
		OUT_HTML += '<td align="center">'+'{0:.2f}'.format(baseCount_byRef[k]/float(refLen_byRef[k]))+'</td>'
		lnk1 = '<div class="ZoomIt"><ul><li><a href="images/coverage'+str(refInd)+'.png"><img class="zit" src="images/coverage'+str(refInd)+'.png" /></a></li>'
		lnk2 = '<li><a href="images/mapQ'+str(refInd)+'.png"><img class="zit2" src="images/mapQ'+str(refInd)+'.png" /></a></li></ul></div>'
		OUT_HTML += '<td>'+lnk1+lnk2+'</td>'
		OUT_HTML += '</tr>'
		refInd   += 1
	OUT_HTML += '</tbody>\n'
	OUT_HTML += '</table>\n'
	OUT_HTML += '</div>\n\n'

	#
	#	overall coverage distribution
	#
	OUT_HTML += '<br>\n'
	OUT_HTML += '<a href="images/totalCoverage.png"><img width="300" height="200" src="images/totalCoverage.png" /></a>\n'

	#
	#	total MapQ distribution
	#
	OUT_HTML += '<a href="images/totalMapQ.png"><img width="300" height="200" src="images/totalMapQ.png" /></a>\n'
	OUT_HTML += '<br>\n'

	#
	#	insert-size distribution
	#
	OUT_HTML += '<a href="images/imd_distribution.png"><img width="300" height="200" src="images/imd_distribution.png" /></a>\n'
	OUT_HTML += '<br>\n\n'

	#
	#	links to log files
	#
	OUT_HTML += '<h2>Logs:</h2>\n'
	OUT_HTML += '<br><a href="logs/bamQC_output.log">bamQC_output.log</a>\n'
	OUT_HTML += '<br><a href="logs/picardValidateSam_output.log">picardValidateSam_output.log</a>\n'
	OUT_HTML += '<br><a href="logs/raw_bamQC_metrics.dat">raw_bamQC_metrics.dat</a>\n'

	#
	#	write output
	#
	OUT_HTML += '<br><br><br>this page was generated on:  <b>'+time.strftime("%Y-%m-%d %H:%M:%S")+'</b><br><br>\n\n'
	of.write(OUT_HTML)
	of.close()


	"""*********************************"""
	"""			GENERATE FIGURES		"""
	"""*********************************"""


	#
	#	various plot-related parameters
	#
	iDir   = outRep+'images/'
	MY_DPI = 80.
	mpl.rcParams.update({'font.size': 18, 'font.weight':'bold', 'lines.linewidth': 2})
	mpl.rcParams['savefig.facecolor'] = '0.92'
	mpl.rcParams['savefig.edgecolor'] = '0.92'
	N_COLORS = 10
	colors = []
	colormap = mpl.cm.Paired(range(256))
	incr = np.linspace(0,255,N_COLORS)
	for i in xrange(N_COLORS):
		colors.append(colormap[incr[i]])
	MAIN_COLOR = colors[1]

	#
	#	readLength_distribution (600x400)
	#
	mpl.figure(0,figsize=(600/MY_DPI,400/MY_DPI))
	x = sorted(nReads_byReadLength.keys())
	y = [nReads_byReadLength[n] for n in x]
	markerline, stemlines, baseline = mpl.stem(x,y)
	mpl.setp(stemlines, linewidth=2, color=MAIN_COLOR)
	mpl.setp(markerline, 'markerfacecolor', MAIN_COLOR)
	mpl_xticks_ints()
	mpl.grid()
	mpl.title('Read Length Distribution')
	mpl.ylabel('# Records')
	mpl.xlabel('Read Length (bp)')
	mpl.tight_layout()
	mpl.savefig(iDir+'readLength_distribution.png')

	#
	#	nFreq_distribution (600x400)
	#
	mpl.clf()
	x = sorted(myNFreq.keys())
	y = [myNFreq[n] for n in x]
	if len(x) <= 1:
		mpl.savefig(iDir+'nFreq_distribution.png')
	else:
		mpl.semilogy(x,y,color=MAIN_COLOR)
		mpl_xticks_ints()
		mpl.grid()
		mpl.title("'N' basecall Frequency")
		mpl.ylabel('# Records')
		mpl.xlabel("# of 'N' Calls in Read")
		mpl.tight_layout()
		mpl.savefig(iDir+'nFreq_distribution.png')

	#
	#	imd_distribution	(600x400)
	#
	mpl.clf()
	x = sorted(myIMD.keys())
	y = [myIMD[n] for n in x]
	if len(x) <= 1:
		mpl.savefig(iDir+'imd_distribution.png')
	else:
		mpl.semilogy(x,y,color=MAIN_COLOR)
		mpl.title('Inner Mate Distance Distribution')
		mpl.ylabel('# Read-Pairs')
		mpl.xlabel('Inner Mate Distance (bp)')
		mpl.grid()
		mpl.tight_layout()
		mpl.savefig(iDir+'imd_distribution.png')

	#
	#	totalCoverage	(600x400)
	#
	totalCov = {}
	for k in refList:
		for k2 in covHist_byRef[k].keys():
			if not k2 in totalCov:
				totalCov[k2] = 0
			totalCov[k2] += covHist_byRef[k][k2]
	mpl.clf()
	x = sorted(totalCov.keys())
	y = [totalCov[n] for n in x]
	if len(x) <= 1:
		mpl.savefig(iDir+'totalCoverage.png')
	else:
		mpl.loglog(x,y,color=MAIN_COLOR)
		mpl.grid()
		mpl.title("Total Per-Base Coverage")
		mpl.ylabel('# Ref Bases')
		mpl.xlabel("x Coverage")
		mpl.tight_layout()
		mpl.savefig(iDir+'totalCoverage.png')

	#
	#	totalMapQ	(600x400)
	#
	totalMapQ = {}
	for k in refList:
		for k2 in mapQHist_byRef[k].keys():
			if not k2 in totalMapQ:
				totalMapQ[k2] = 0
			totalMapQ[k2] += mapQHist_byRef[k][k2]
	mpl.clf()
	x = sorted(totalMapQ.keys())
	y = [totalMapQ[n] for n in x]
	if len(x) <= 1:
		mpl.savefig(iDir+'totalMapQ.png')
	else:
		mpl.semilogy(x,y,color=MAIN_COLOR)
		mpl.grid()
		mpl.title("Total MAPQ Distribution")
		mpl.ylabel('# Records')
		mpl.xlabel("Mapping Quality")
		mpl.tight_layout()
		mpl.savefig(iDir+'totalMapQ.png')

	#
	#	coverage#	(500x400)
	#
	mpl.rcParams['font.size'] = 14
	mpl.figure(1,figsize=(500/MY_DPI,400/MY_DPI))
	refInd = 0
	for k in refList:
		mpl.clf()
		x = sorted(covHist_byRef[k].keys())
		y = [covHist_byRef[k][n] for n in x]
		if len(x) <= 1:
			mpl.savefig(iDir+'coverage'+str(refInd)+'.png')
		else:
			mpl.loglog(x,y,color=MAIN_COLOR)
			mpl.grid()
			mpl.title("Per-Base Coverage Distribution")
			mpl.ylabel('# Ref Bases')
			mpl.xlabel("x Coverage")
			mpl.tight_layout()
			mpl.savefig(iDir+'coverage'+str(refInd)+'.png')
		refInd += 1

	#
	#	mapQ#	(500x400)
	#
	refInd = 0
	for k in refList:
		mpl.clf()
		x = sorted(mapQHist_byRef[k].keys())
		y = [mapQHist_byRef[k][n] for n in x]
		if len(x) <= 1:
			mpl.savefig(iDir+'mapQ'+str(refInd)+'.png')
		else:
			markerline, stemlines, baseline = mpl.stem(x,y)
			mpl.setp(stemlines, linewidth=2, color=MAIN_COLOR)
			mpl.setp(markerline, 'markerfacecolor', MAIN_COLOR)
			mpl.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			mpl.grid()
			mpl.title("Mapping Quality Distribution")
			mpl.ylabel('# Records')
			mpl.xlabel("Mapping Quality Score")
			mpl.tight_layout()
			mpl.savefig(iDir+'mapQ'+str(refInd)+'.png')
		refInd += 1


	print ' done.  '#,time.time()-tt,'(sec)'


if __name__ == '__main__':
	main()


