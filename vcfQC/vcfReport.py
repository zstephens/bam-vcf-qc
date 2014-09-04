#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

python vcfReport.py in.vcf vcfTools_validate.log vcfTools_stats.dat outReport/

************************************************** """

import os
import sys
import re
import time
import numpy as np
import matplotlib.pyplot as mpl

# if a line from the vcftools-validate log contains these phrases, it's indicative of a fatal error
FAIL_MESSAGE_PREFIXES = ["Could not parse the line:",
						 "Could not parse the line, wrong number of columns:",
						 "Sorry, empty lines not allowed.",
						 "Multiple header blocks ",
						 "Wrong number of fields",
						 "Could not parse the header line:",
						 "Broken VCF header, no column names?",
						 "Could not parse the mandatory columns",
						 "Expected key=value pair in the header:",
						 "Missing the ID tag in ",
						 "Could not parse ALT ",
						 "Could not parse the allele",
						 "cannot be used, it is a reserved word.",
						 "Unable to parse the GT field"]

NUCL_IND = {'A':0, 'C':1, 'G':2, 'T':3}

# 0: transition, 1: transversion
TS_TV = [[9,1,0,1],
		 [1,9,1,0],
		 [0,1,9,1],
		 [1,0,1,9]]

def main():

	if len(sys.argv) != 5:
		print '\npython vcfReport.py in.vcf vcfTools_validate.log vcfTools_stats.dat outReport/\n'
		exit(1)
	else:
		inVcf = sys.argv[1]
		vlog = sys.argv[2]
		slog = sys.argv[3]
		outRep = sys.argv[4]

	inf = open('reportHeaderVcf.html','r')
	REPORT_HEADER = inf.read()
	inf.close()

	allErrors = []
	f = open(vlog,'r')
	for line in f:
		failWorthy = 0
		for fmp in FAIL_MESSAGE_PREFIXES:
			if len(line) >= len(fmp) and fmp in line:
				failWorthy = 1
		allErrors.append([line,failWorthy])
	f.close()

	nFail = sum([n[1] for n in allErrors])

	f = open(slog,'r')
	fr = f.read()
	f.close()

	noStats = False
	# if we don't have stats...
	if len(fr) == 0:
		noStats = True
	else:
		nTotal = int(re.findall(r"'count' => (.*?)(?=,)",fr)[0])
		nSNP   = int(re.findall(r"'snp_count' => (.*?)(?=,)",fr)[0])
		nIndel = int(re.findall(r"'indel_count' => (.*?)(?=,)",fr)[0])
		nOther = int(re.findall(r"'other_count' => (.*?)(?=,)",fr)[0])

		pat = re.compile(r"'indel' => {(.*?)(?=})", re.MULTILINE | re.DOTALL)
		indels = re.findall(pat,fr)[0].split('\n')[1:-1]
		indLen = {}
		for n in indels:
			splt = n.split("'")
			indLen[int(splt[1])] = int(splt[2].strip(',')[4:])

		pat = re.compile(r"'snp' => {(.*?)(?=})", re.MULTILINE | re.DOTALL)
		snps = re.findall(pat,fr)[0].split('\n')[1:-1]
		snpMap = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
		ts = 0
		tv = 0
		for n in snps:
			splt = n.split("'")
			i1 = NUCL_IND[splt[1][0]]
			i2 = NUCL_IND[splt[1][2]]
			nSnp = int(splt[2].strip(',')[4:])
			snpMap[i1][i2] += nSnp
			if TS_TV[i1][i2]:
				tv += nSnp
			else:
				ts += nSnp
		tstv = float(ts)/float(tv)


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

	of = open(outRep+'vcfQC_report.html','w')

	REPORT_HEADER = re.sub("<title>-","<title>"+"VcfQC Report",REPORT_HEADER)
	of.write(REPORT_HEADER)

	OUT_HTML = '<h2><font color="#444444">'+inVcf+'</font></h2>\n\n'

	#
	#	QC warnings
	#
	OUT_HTML += '<h2>Quality Control Warnings/Failures:</h2>\n'
	nStr = str(nFail)
	if nFail > 0:
		nStr = '<span class="failText">'+nStr+'</span>'

	OUT_HTML += '<div class="datagrid">\n<table>\n<tbody>\n'
	OUT_HTML += '<tr class="alt"><td>Fatal Errors:</td><td>'+nStr+'</td></tr>\n'
	OUT_HTML += '</tbody>\n</table>\n</div>\n\n'

	OUT_HTML += '<br><div class="datagrid">\n<table>\n<tbody>\n'
	isEven = True
	for n in allErrors:
		if n[1]:
			if isEven:
				OUT_HTML += '<tr class="alt"><td>'+n[0]+'</td></tr>\n'
			else:
				OUT_HTML += '<tr><td>'+n[0]+'</td></tr>\n'
			isEven = not(isEven)
	OUT_HTML += '</tbody>\n</table>\n</div>\n\n'

	OUT_HTML += '<br><a href="logs/vcfTools_validate.log">vcfTools_validate.log</a>\n'
	OUT_HTML += '<br><a href="logs/vcftools_stats.errors">vcftools_stats.errors</a>\n'

	if not noStats:
		#
		#	basic stats table
		#
		OUT_HTML += '<h2>Basic Stats:</h2>\n'
		OUT_HTML += '<div class="datagrid">\n'
		OUT_HTML += '<table>\n'
		OUT_HTML += '<tbody>\n'
		OUT_HTML += '<tr class="alt"><td><b>'+'Total Variants:'+'</b></td><td>'+'{0:,}'.format(nTotal)+'</td></tr>\n'
		OUT_HTML += '<tr><td><b>'+'Total SNPs:'+'</b></td><td>'+'{0:,}'.format(nSNP)+'</td></tr>\n'
		OUT_HTML += '<tr class="alt"><td><b>'+'Total Indels:'+'</b></td><td>'+'{0:,}'.format(nIndel)+'</td></tr>\n'
		OUT_HTML += '<tr><td><b>'+'Total Other:'+'</b></td><td>'+'{0:,}'.format(nOther)+'</td></tr>\n'
		OUT_HTML += '<tr class="alt"><td><b>'+'Ts/Tv:'+'</b></td><td>'+'{0:.2f}'.format(tstv)+'</td></tr>\n'
		OUT_HTML += '</tbody>\n'
		OUT_HTML += '</table>\n'
		OUT_HTML += '</div>\n\n'

		#
		#	fun little plots
		#
		OUT_HTML += '<br><a href="images/snpMap.png"><img width="300" height="200" src="images/snpMap.png" /></a>\n'
		OUT_HTML += '<a href="images/indels.png"><img width="300" height="200" src="images/indels.png" /></a>\n'
		OUT_HTML += '<br>\n\n'

	#
	#	write output
	#
	OUT_HTML += '<br><br><br>this page was generated on:  <b>'+time.strftime("%Y-%m-%d %H:%M:%S")+'</b><br><br>\n\n'
	of.write(OUT_HTML)
	of.close()


	"""*********************************"""
	"""			GENERATE FIGURES		"""
	"""*********************************"""


	if not noStats:
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
		#	snpMap (600x400)
		#
		for i in xrange(len(snpMap)):
			snpMap[i] = [float(n)/sum(snpMap[i]) for n in snpMap[i]]
		snpMap = np.array(snpMap)
		mpl.figure(0,figsize=(600/MY_DPI,400/MY_DPI))
		Z = snpMap
		X, Y = np.meshgrid( range(0,len(Z[0])+1), range(0,len(Z)+1) )
		mpl.pcolormesh(X,Y,Z[::-1],vmin=0.,vmax=0.5)
		mpl.axis([0,len(Z[0]),0,len(Z)])
		mpl.yticks(np.arange(0,len(Z))+0.5,['A','C','G','T'][::-1])
		mpl.xticks(np.arange(0,len(Z[0]))+0.5,['A','C','G','T'])
		mpl.ylabel('n1')
		mpl.xlabel('n2')
		mpl.title('SNP transitions: n1 -> n2')
		cbar = mpl.colorbar()
		cbar.set_ticks([0,.1,.2,.3,.4,.5])
		cbar.set_ticklabels(['0.0','0.1','0.2','0.3','0.4','0.5'])
		mpl.savefig(iDir+'snpMap.png')

		#
		#	indels (600x400)
		#
		mpl.clf()
		x = sorted(indLen.keys())
		y = [indLen[n] for n in x]
		mpl.plot(x,y)
		mpl.ylabel('# indels')
		mpl.xlabel('<-- del length | ins length -->')
		mpl.title('Indel Length Distribution')
		mpl.savefig(iDir+'indels.png')


		if nFail > 0:
			exit(1)


if __name__ == '__main__':
	main()


