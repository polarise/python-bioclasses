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
		
	# load the indexed FASTA file
	fastafile = pysam.Fastafile( "/home/paul/Resources/H_sapiens/hg19.ens.chr.fa" )
	
	# iterate over a few Transcript objects
	c = 0
	for gene_id,G in genes.iteritems():
		for transcript_id,T in G.transcripts.iteritems():
			if c > 100:
				break
			
			print T
			print T.UTR
			
			if len( T.UTR ) > 0:
				print T.compute_UTR_GC_content( fastafile )
			print
						
		c += 1
	
	# close the fastafile
	fastafile.close()

if __name__ == "__main__":
	main()