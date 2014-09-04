#!/usr/bin/env python
# encoding: utf-8

""" **************************************************

python EOF_stream.py

************************************************** """

import sys

CORRECT_EOF = [31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0]

READ_CHUNK = 4096

def main():
	
	eofLen  = len(CORRECT_EOF)

	print 'checking EOF...',
	prevChunk = ''
	count = 0
	while True:
		ch = sys.stdin.read(READ_CHUNK)
		if len(ch) < READ_CHUNK: # EOF
			break

		prevChunk = ch
		count += 1
		#sys.stdout.write(ch)

	myEOF = [ord(n) for n in (prevChunk+ch)[-eofLen:]]

	if myEOF == CORRECT_EOF:
		print 'pass.'
		exit(0)
	else:
		print 'fail.'
		exit(1)

	#sys.stdout.write(ch)
	#sys.stdout.close()


if __name__ == '__main__':
	main()