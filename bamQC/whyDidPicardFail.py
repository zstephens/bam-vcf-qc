#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

python whyDidPicardFail.py picardErrorLog picardErrorLog picardErrorOutput

************************************************** """

import sys

def main():

	if len(sys.argv) != 3:
		print '\npython whyDidPicardFail.py picardErrorLog picardErrorOutput\n'
		exit(1)
	else:
		inF = sys.argv[1]
		ouF = sys.argv[2]

	#
	#	scan through the picard error log, if we find a reason we know it failed, make note of it for the BamQC report
	#
	f = open(inF,'r')
	outputErrors = []
	for line in f:

		#
		#	did it fail because duplicate reads present? (i.e. multiple alignments present in bam, or merged bam from multiple runs)
		#
		if "PicardException: Value was put into PairInfoMap more than once." in line:
			outputErrors.append('ERROR:Picard cannot handle BAM with duplicate reads, sorry.\tn/a\n')

	f.close()

	#if len(outputErrors) == 0:
	#	outputErrors.append('ERROR:Picard failed for an unknown reason, check log.\tn/a\n')
	
	outF = open(ouF,'a')
	for n in outputErrors:
		outF.write(n)
	outF.close()


if __name__ == '__main__':
	main()


