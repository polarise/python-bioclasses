#!/usr/bin/env python
from __future__ import division
import sys
import gzip
from Sequence import *

def main():
	with open( "/home/paul/bioinf/Translational_Frameshifting/2_Euplotes_data/ALLCDS.fa" ) as f:
	#with open( "/home/paul/bioinf/Translational_Frameshifting/1_human_data/some_genes.fasta" ) as f:
	#with gzip.open( "/home/paul/Downloads/mart_export.txt.gz" ) as f:
		c = 0
		read_header = False
		for row in f:
			if c > 10:
				break
			if row[0] == "C":
				print row.strip( "\n" )
				continue
			sequence = ""
			while row[0] != "C":
				sequence += row.strip( "\n" )
				row = f.next()
			min_len = 999
			print sequence
			if len( sequence ) > min_len:
				for m in xrange( len( sequence )//min_len + 1 ):
					t = sequence[m*min_len:(m+1)*min_len]
					s = BiologicalSequence( t )
					print s.info()
					s.build_tree()
					print "There are %d paths." % len( s.paths )
					for frame in xrange( 3 ):
						print s.binary_frame( frame, sep="" )
					for sp in s.sorted_frame_paths:
						print sp, s.sorted_frame_paths[sp]
					print
			else:
				s = BiologicalSequence( sequence )
				print s.info()
				s.build_tree()
				print "There are %d paths." % len( s.paths )
				for frame in xrange( 3 ):
					print s.binary_frame( frame, sep="" )
				for sp in s.sorted_frame_paths:
					print sp, s.sorted_frame_paths[sp]
				print
			print
			print row.strip( "\n" )
			c += 1
	
	#print len( s.frame_paths[0] )
	#print len( s.frame_paths[1] )
	#print len( s.frame_paths[2] )
	
	#for l in s.sorted_frame_paths:
		#print "Frameshifts of cardinality %d:" % l
		#for paths in s.sorted_frame_paths[l]:
			#print paths
		#print

if __name__ == "__main__":
	main()