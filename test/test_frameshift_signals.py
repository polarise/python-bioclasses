#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
from Sequence import *
from TransitionMatrix import *
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO

def main( fn ):
	TM = TransitionMatrix()
	TM.read( "euplotid_transition_matrix.pic" )
	pdf = PdfPages( "likelihood_profiles_test.pdf" )
	b = 0 # count the ones that pass
	c = 0 # count all
	for seq_record in SeqIO.parse( fn, "fasta" ):
		if c > 1000:
			break
		sequence = str( seq_record.seq )
		seq_name = seq_record.id
		s = Sequence( sequence=sequence, name=seq_name )
		s.truncate( effect_truncation=True, verbose=False )
		no_of_leaves = s.count_leaves()
		if no_of_leaves > 1000:
			print >> sys.stderr, "Complex tree with %s leaves...omitting." % no_of_leaves
			continue
		s.set_transition_matrix( TM )
		s.build_tree()
		s.get_frameshift_signals()
		s.estimate_likelihood()
		s.estimate_frameshift_likelihood()
		s.get_most_likely_frameshift()
		if s.most_likely_frameshift is not None:
			if 1 < len( s.most_likely_frameshift.path ) < 4:
				#s.plot_differential_graded_likelihood( outfile=pdf, show_path_str=True )
				s.plot_differential_graded_likelihood()
				b += 1
		c += 1
	pdf.close()
	
	print >> sys.stderr, "Processed %d (of %d) sequences [%.2f%%]." % ( b, c, b/c*100 )
	
if __name__ == "__main__":
	try:
		fn = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "./script.py <fasta-file>"
		sys.exit( 1 )
	main( fn )
