#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import cPickle
import BioClasses

def main( fn ):
	with open( fn ) as f:
		genes = cPickle.load( f )
	
	c = 0
	for gene_id,G in genes.iteritems():
		for transcript_id,T in G.transcripts.iteritems():
			if c > 99:
				break
			if T.is_protein_coding():
				print "\t".join( map( str, [ gene_id, transcript_id, T.region_str(), T.strand, T.length_of_spliced_5utr(), T.get_cds_length()]))
			c += 0

if __name__ == "__main__":
	main( sys.argv[1] ) # hg19.chr.pic