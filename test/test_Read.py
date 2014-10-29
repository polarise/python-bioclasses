#!/usr/bin/env python
from __future__ import division
import sys
from Read import *

def main( fq_file, n ):
	pos = 0
	with open( fq_file ) as f:
		for row in f:
			l = row.strip()
			if pos % 4 == 0:
				R = Read( l, trim_by=n )
			elif pos % 4 == 1:
				R.bases = l
			elif pos % 4 == 2:
				R.name_ = l
			elif pos % 4 == 3:
				R.quals = l
				R.remove_polyA()
				try:
					print R
				except TypeError:
					pass
			pos += 1

if __name__ == "__main__":
	try:
		fq_file = sys.argv[1]
		n = int( sys.argv[2] )
	except IndexError:
		raise IOError( "Missing FASTQ input file or INT trim." )
	main( fq_file, n )