#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import sys
import BioClasses
import cPickle
import pysam
import scipy

def main():
	# load data from the GTF pic
	with open( "/home/paul/Resources/H_sapiens/test.hg19.chr.pic" ) as f:
		genes = cPickle.load( f )
	
	# iterate over a few Transcript objects
	c = 0
	for gene_id,G in genes.iteritems():
		for transcript_id,T in G.transcripts.iteritems():
			if c > 100:
				break
			
			print T,T.strand
			# iterate through the exons
			exons = T.exons.values()
			exons.sort()
			i = 0
			for exon_id,E in T.exons.iteritems():
				print exon_id,E,exons[i]
				i += 1
						
		c += 1

if __name__ == "__main__":
	main()