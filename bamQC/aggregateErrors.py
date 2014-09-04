#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

python aggregateErrors.py bamQCData bamQCErrors picardErrors output

************************************************** """

import sys
import re

# OUTPUT FORMAT: [program performing check, check name, severity, # records affected]

# (severity, affects_all_reads, description)
BAMQC_ERROR_CODES = {0:  ('fail',True,'Invalid Magic String'),
					 1:  ('warn',True,'Unknown Header Line'),
					 2:  ('fail',True,'Invalid @HD Line'),
					 3:  ('fail',True,'Invalid @HD Line'),
					 4:  ('fail',True,'Invalid @HD Line'),
					 5:  ('fail',True,'Invalid @SQ Line'),
					 6:  ('fail',True,'Invalid @SQ Line'),
					 7:  ('fail',True,'Unknown Header Tag'),
					 8:  ('fail',True,'Invalid Ref Length'),
					 9:  ('fail',True,'Unknown Reference'),
					 10: ('fail',True,'Unknown Reference'),
					 11: ('fail',True,'Unknown Tag Datatype'),
					 12: ('fail',True,'Unknown Tag Datatype'),
					 13: ('fail',True,'Invalid Ref ID'),
					 14: ('warn',True,'Ambiguous Multiplicity (Adjust your faith in multi-mapped read counts accordingly!)'),
					 15: ('warn',False,'Invalid Mapped Pos'),
					 16: ('fail',True,'Sorting Error')}

PICARD_ERRORS = {'INVALID_QUALITY_FORMAT' : 'warn',
				 'INVALID_FLAG_PROPER_PAIR' : 'info',
				 'INVALID_FLAG_MATE_UNMAPPED' : 'info',
				 'MISMATCH_FLAG_MATE_UNMAPPED' : 'info',
				 'INVALID_FLAG_MATE_NEG_STRAND' : 'info',
				 'MISMATCH_FLAG_MATE_NEG_STRAND' : 'info',
				 'INVALID_FLAG_FIRST_OF_PAIR' : 'info',
				 'INVALID_FLAG_SECOND_OF_PAIR' : 'info',
				 'PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND' : 'warn',
				 'INVALID_FLAG_NOT_PRIM_ALIGNMENT' : 'info',
				 'INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT' : 'info',
				 'INVALID_FLAG_READ_UNMAPPED' : 'info',
				 'INVALID_INSERT_SIZE' : 'info',
				 'INVALID_MAPPING_QUALITY' : 'warn',
				 'INVALID_CIGAR' : 'warn',
				 'ADJACENT_INDEL_IN_CIGAR' : 'warn',
				 'INVALID_MATE_REF_INDEX' : 'fail',
				 'MISMATCH_MATE_REF_INDEX' : 'warn',
				 'INVALID_REFERENCE_INDEX' : 'fail',
				 'INVALID_ALIGNMENT_START' : 'fail',
				 'MISMATCH_MATE_ALIGNMENT_START' : 'warn',
				 'MATE_FIELD_MISMATCH' : 'warn',
				 'INVALID_TAG_NM' : 'info', 
				 'MISSING_TAG_NM' : 'info',
				 'MISSING_HEADER' : 'fail',
				 'MISSING_SEQUENCE_DICTIONARY' : 'fail',
				 'MISSING_READ_GROUP' : 'warn',
				 'RECORD_OUT_OF_ORDER' : 'fail',
				 'READ_GROUP_NOT_FOUND' : 'warn',
				 'RECORD_MISSING_READ_GROUP' : 'warn',
				 'INVALID_INDEXING_BIN' : 'warn',
				 'MISSING_VERSION_NUMBER' : 'warn',
				 'INVALID_VERSION_NUMBER' : 'warn',
				 'TRUNCATED_FILE' : 'fail',
				 'MISMATCH_READ_LENGTH_AND_QUALS_LENGTH' : 'fail',
				 'EMPTY_READ' : 'info',
				 'CIGAR_MAPS_OFF_REFERENCE' : 'warn',
				 'MISMATCH_READ_LENGTH_AND_E2_LENGTH' : 'info',
				 'MISMATCH_READ_LENGTH_AND_U2_LENGTH' : 'info',
				 'E2_BASE_EQUALS_PRIMARY_BASE' : 'info',
				 'BAM_FILE_MISSING_TERMINATOR_BLOCK' : 'fail',
				 'UNRECOGNIZED_HEADER_TYPE' : 'fail',
				 'POORLY_FORMATTED_HEADER_TAG' : 'warn',
				 'HEADER_TAG_MULTIPLY_DEFINED' : 'warn',
				 'HEADER_RECORD_MISSING_REQUIRED_TAG' : 'fail',
				 'INVALID_DATE_STRING' : 'warn',
				 'TAG_VALUE_TOO_LARGE' : 'warn',
				 'INVALID_INDEX_FILE_POINTER' : 'fail',
				 'INVALID_PREDICTED_MEDIAN_INSERT_SIZE' : 'info',
				 'DUPLICATE_READ_GROUP_ID' : 'warn',
				 'MISSING_PLATFORM_VALUE' : 'warn',
				 'INVALID_PLATFORM_VALUE' : 'warn',
				 'DUPLICATE_PROGRAM_GROUP_ID' : 'info',
				 'MATE_NOT_FOUND' : 'warn',
				 'MATES_ARE_SAME_END' : 'warn',
				 'MISMATCH_MATE_CIGAR_STRING' : 'warn',
				 'MATE_CIGAR_STRING_INVALID_PRESENCE' : 'info',

				 'Picard cannot handle BAM with duplicate reads, sorry.': 'skip',
				 'Picard failed for an unknown reason, check log.': 'skip'}

def weighted_N_content(Nfreq):	#	give warning if above ~0.25?
	nDist = 0.
	for i in xrange(len(Nfreq)):
		nDist += i*Nfreq[i]
	return nDist

N_THRESH  = 0.25
N_WARNING = ('BamQC','High N Content','warn','n/a')

def main():

	if len(sys.argv) != 5:
		print '\npython aggregateErrors.py bamQCData bamQCErrors picardErrors output\n'
		exit(1)
	else:
		inBD = sys.argv[1]
		inBE = sys.argv[2]
		inPE = sys.argv[3]
		outF = sys.argv[4]

	allErrors    = []
	bamQC_errors = {}
	f = open(inBE,'r')
	for line in f:
		errorCode = int(re.findall(r"\[\d+\]",line)[0][1:-1])
		if errorCode not in bamQC_errors:
			res = BAMQC_ERROR_CODES[errorCode]
			bamQC_errors[errorCode] = ['BamQC',res[2],res[0],0]
			if res[1]:
				bamQC_errors[errorCode][3] = 'n/a'
		if bamQC_errors[errorCode][3] != 'n/a':
			bamQC_errors[errorCode][3] += 1
	f.close()
	for k in sorted(bamQC_errors.keys()):
		allErrors.append(bamQC_errors[k])

	f = open(inPE,'r')
	for line in f:
		if len(line) >= 6:
			if line[:6] == 'ERROR:' or line[:8] == 'WARNING:':
				splt = line[:-1].split(':')[1]
				splt = splt.split('\t')
				allErrors.append(['Picard ValidateSam',splt[0],PICARD_ERRORS[splt[0]],splt[1]])
	f.close()


	""" *****************************************************************
						PERFORM SOME ADDITIONAL CHECKS
	***************************************************************** """


	#
	#	N-frequency check
	#
	myNFreq = {}
	f = open(inBD,'r')
	for line in f:
		splt = line[:-1].split('\t')
		if splt[0][:2] == '>>':
			read_nf = False

		if splt[0] == '>>TOTAL_UNIQUE_IDS':
			totalUniqueIDs = int(splt[1])
		elif splt[0] == '>>N_FREQUENCY':
			read_nf = True

		elif splt[0] != '>':
			if read_nf:
				myNFreq[int(splt[0])] = int(splt[1])/float(totalUniqueIDs)
	f.close()
	nFreqList = [0 for n in xrange(max(myNFreq.keys())+1)]
	for k in myNFreq.keys():
		nFreqList[k] = myNFreq[k]
	if weighted_N_content(nFreqList) > N_THRESH:
		allErrors.insert(0,N_WARNING)


	weFailed = 0
	f = open(outF,'w')
	for n in allErrors:
		if n[2] == 'fail':
			weFailed = 1
		f.write(str(n[0])+'\t'+str(n[1])+'\t'+str(n[2])+'\t'+str(n[3])+'\n')
	f.close()


	if weFailed:
		exit(1)


if __name__ == '__main__':
	main()


