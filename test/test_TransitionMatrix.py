#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from TransitionMatrix import *
from Sequence import *
import numpy
import matplotlib.pyplot as plt
import math
from Bio import SeqIO
from LeafCounter import *
from matplotlib.backends.backend_pdf import PdfPages

def f( x ):
	return -(x - 1)*math.log( 64 )/3

def main( fn ):
	# initialise the TransitionMatrix
	TM = TransitionMatrix()
	# T.build( fn )
	# T.write( "euplotid_transition_matrix.pic" )
	TM.read( "euplotid_transition_matrix.pic" )
	
	pdf = PdfPages( "likelihood_profiles_1000fs.pdf" )
	
	c = 0
	for seq_record in SeqIO.parse( fn, "fasta" ):
		if c > 10:
			break
		sequence = str( seq_record.seq )
		seq_name = seq_record.id
		
		s = Sequence( sequence )
		s.truncate()
		s.get_stop_sequence()
		s.sanitise_stop_sequence()
		nodes = [ Node( *d ) for d in s.unique_stop_sequence ]
		L = LeafCounter()
		for n in nodes[:-3]:
			L.add_node( n )

		if L.leaf_count() > 1000:
			print >> sys.stderr, "Skipping complex sequence %s with %d leaves..." % ( seq_name, L.leaf_count())
			continue
		
		# set params
		s.set_transition_matrix( TM )
		s.build_tree()
		s.estimate_likelihood()
		s.estimate_frameshift_likelihood()
		s.get_most_likely_frameshift()
	
		if s.most_likely_frameshift is None:
			print >> sys.stderr, "%s admits no frameshifts..." % seq_name
			continue
		
		# original sequence
		y = numpy.linspace( 1, s.length, len( s.graded_likelihood ))
		plt.plot( y, map( lambda w: w[1] - f( w[0] ), zip( y, s.graded_likelihood )), ':', label="Original (fr 0)" )
		
		for fs in s.frameshift_sequences:
			F = s.frameshift_sequences[fs]
			y = numpy.linspace( 1, F.length, len( F.graded_likelihood ))
			# plt.figure( figsize=( 10, 5 ))
			if F == s.most_likely_frameshift:
				plt.plot( y, map( lambda w: w[1] - f( w[0] ), zip( y, F.graded_likelihood )), label="ML (fr %s)" % F.path[0][0] )
			else:
				plt.plot( y, map( lambda w: w[1] - f( w[0] ), zip( y, F.graded_likelihood )))
		
		plt.title( "%s\n(best of %d frameshift sequences)" % ( seq_name, len( s.frameshift_sequences )))
		
		# draw all the frameshift sites
		most_likely_signals = s.most_likely_frameshift.signals
		ymin, ymax = plt.ylim()
		i = 0
		for frame,position in s.most_likely_frameshift.path:
			if position >= 0:
				plt.axvline( x=position ) 
				plt.annotate( most_likely_signals[i], xy=( position+1, ymax - 4 ), rotation=90,  )
				i += 1

		plt.legend( loc="best", prop={'size': "small"} )
		# plt.grid()
		pdf.savefig()
		plt.close()
		
		c += 1
	
	pdf.close()

	# sequence = sequence.replace( "\n", "" )
	# print sequence
	# s = Sequence( sequence )
	# s.truncate()
	# print s
	# s.set_transition_matrix( TM )
	# s.build_tree()
	# print s.tree
	# print
	# print s.unique_stop_sequence
	# print
	# s.estimate_likelihood()
	# s.estimate_frameshift_likelihood()
	# for fs in s.frameshift_sequences:
	# 	print fs, s.frameshift_sequences[fs].likelihood
	#
	#
	#
	# print
	# s.get_most_likely_frameshift()
	# print "Most likely frameshift:\n%s" % s.most_likely_frameshift
	# print
	# print "Least likely frameshift:\n%s" % s.least_likely_frameshift
	# print
	#
	# likelihood_values = [ s.frameshift_sequences[fs].likelihood for fs in s.frameshift_sequences ]
	# print "Likelihood range: %s to %s" % ( min( likelihood_values ), max( likelihood_values ))
	# print
	# # print "All frameshifts:"
	# # for fs in s.frameshift_sequences:
	# # 	print s.frameshift_sequences[fs]
	# # 	print
	#
	# x = numpy.linspace( 1, s.length, len( s.graded_likelihood ) )
	# plt.plot( x, map( lambda w: w[1] - f( w[0] ), zip( x, s.graded_likelihood )))
	# # plt.plot( x, f( x ), '--' )
	# for fs in s.frameshift_sequences:
	# 	F = s.frameshift_sequences[fs]
	# 	y = numpy.linspace( 1, F.length, len( F.graded_likelihood ))
	# 	plt.plot( y, map( lambda w: w[1] - f( w[0] ), zip( y, F.graded_likelihood )))
	#
	# # draw all the frameshift sites
	# most_likely_signals = s.most_likely_frameshift.signals
	# ymin, ymax = plt.ylim()
	# i = 0
	# for frame,position in s.most_likely_frameshift.path:
	# 	if position >= 0:
	# 		plt.axvline( x=position+1 ) # 0-based to 1-based
	# 		plt.annotate( most_likely_signals[i], xy=( position+1, ymax - 4 ), rotation=90 )
	# 		i += 1
	#
	# plt.grid()
	# plt.show()

if __name__ == "__main__":
	fn = sys.argv[1]
	main( fn )