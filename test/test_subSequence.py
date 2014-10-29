#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
from Sequence import *
from TransitionMatrix import *
from FrameshiftSequence import *

def main():
	# initialise the transition matrix
	TM = TransitionMatrix()
	TM.read( "transition_matrices/euplotid_transition_matrix.pic" )
	
	# start off with a sequence
	# comp1
	sequence = """ATCGGAAGAGGCAATTTCATAATAACATTTACACTAAGCTGAAGCCCAAGGAAATGAACT
CCCATTTTAGCGATTATAATGAGACCAGGTATGAAGTCACTACTGGGTTTAATAGCAAGA
AGATCAAACTTAAGGTTAGTACTAACACCGAAGGGGATGAACAGAGCAAGCCTATGAATA
AGAGGACTTTTAAGAACGATATACTTGAGAAGGAGGAGGATCCCACTGGGAAGAAAACTG
TGGTGTATAAGGATTTTATAGACTTTTGTAAAATTCATAAGGGGCAGATTTTTGGGTATA
GAACTATCATGCCTCTTGAGTTTTATATTATGTCAAAATAGAACCGATTATGAGAGGGAG
TTTATAAAAAGGGCTGAGAAGGAAGAGATTCAGAAGTTCTTTAATGAAGCATGAGTATCC
ATTGTGGCTAACTCAGCAATCGTTGAGACCTTTATCTTTGAAAAGGCTTTAATGTCTTTT
CTACCAGAACACCTGTCTCAGGCTTTCTTCAAGGATTTACTGAATTGCAAAGAGCATGAC
AGGCCTGCTAATATCGTTGAGAAAGGCAAGAACGATAGTATTATTGGAAGTATTATCAAG
ATGAATAAAGAAGATGATTTATGGGATAATACTAAAAACAAAATTATTGATAAGACGCTC
AAAGCGTCCTATATCGAGAGACACAAGGCTTTAG"""
	
	sequence = """TATAAATTATTGAGTACATAATGAAATGTTTCTGAGGGCAAGTTCGTTTAACCCCTTATG
GAGAACTAGCTTTCTAATGACTAACTACTTGAGACCAGCAACTTCACAAATTTCTCTGTT
TTCTTTAGTACAAAAATCGATGGTGCAGCTGATATAGACTGGGATGCTATCGATTCGAAA
TCTATTAAAAAGTTAACAAGTTATCATATGTAGAGTAAATTGTACAATCAGTTGTCTGAC
CAATAACAAAACCACTTGCCATAAAGCAGGCAAAGAAAACCAGGAGACTTGAAAGAGAAG
TCAATAGAGAACCAATATCAGAGGAAAATCAGAAAATTGTAGCTCATTGTGAAGTATAGC
AAATCACCTATCCCCAAGACCCTTATGGGATATTTGA"""

	sequence = """ATGAGAAGAGCAAGAAAGCTATTGAGGTCAGAAGGATCATAGTAACTCCTTCTAGAATGG
AGCTTGAGTTTCCATCTCAATCTGTTCCGAATAGGGCTATCAGGACTCACTTAAAATAAT
CTTGATGATTTTGTTTGTGTTTCTTTCCAAAATGAAGATCATCAGAAAGGAATGTATGCA
TACAAATCCGATGAAACTAATGATGTCTTGAATCACATCGAGCATACTTTAAAGGAAGGA
TTTCTCATCGGCAGTAAAAGATTTAAGTTCTTACATTACTCAAACTCGCAAATGAAAAGT
TACTCCTGCTGGTTTATGAATGAGGTTTCTCCCCACTT"""
	# comp12974
	sequence = """ATATAATTATTGCAAGTTTGAAATTGATGGAGAAATTAAGATCATTTGTTACTGTGAATC
TTTATTCAAAGATCATTGCAAGTTGTGTGCTCTTACTGACCATTCATACTTTTAACAGTG
AGCAATATCAAGGAGTCAAGGATTTTTGAAAATTGTTGGATATAGGAGATGAAGTATCAA
AGTGATTGACCTATTTTCTTGGAAGTTTAACATTTTTTGCAGCAATCTCTATATTGTTCA
GATGGTATAAAAGGGCAATGTATTTATGAATAATAGCATGAGGATGAAATCTGTTCGCAC
TTTTATGAGTTGGAATAACAAAAAATGATAATAACTCTGACCCTCCTCTCAGCCTGACTC
ATGCCCTCATCTCGCACATTCTACATCCAAACTCAGTTTGCCTTGAAAACGTGCTGACTT
GAACTGCAGTCCACTCCTGAATGTTCTATGTAGCTTACCAGTACATATAGGGTCTTGGGA
AGTAGAAGCAGGTGAAGAGGAGGGATTTTTAATA"""
	# comp16061
	#sequence = """TTTTCGCAGAAATACCTTTGCCAAACAACAAGATTTATGATGCATATTCAGCTCTATACC
#CAAAAGTGATGAACAAGAAGTTTGGGATACCCATAAATACCCCAAAACCAGTAGTAAAGC
#CAATAGTTGGCTTCCAGCTTAGATCTCAGCGTATAAAATCTCGGAATACTCCCCTTGTCA
#ATAACTTGACTCCTTCAAATATGAATGAGGAGAAAAAGTATATTCAGAATATACATGATA
#TCAGGAGTGAGTTGAGGAAATTTAATCATTCAAAGGACTCAACTGCTTTATTCCCCACTT
#ATCGCGCTAATGATACTGAACATGTTATGACGATTGAGAGGAACCAAGATAACTCACCTC
#TAAGAGATGACAGAATAAGGACCACTCTGAGCCATGATACTCAAACTATTTCTACAAAGA
#AGATTTTACGGACTTCTTCCAAGGCTGATAGCACTAAGTCTCGGTATTAGGCTCTCGAGA
#ACATATGAGTTCTAGCATACAAAATAAGCTGACATATAAGAAAGCATGGGAGGAAGAATC
#CAATGAAGAATCTAAGGGATGTGATGACTTTAACCAGGCTCAAAGATATTCAGAGAAATT
#CTTGCTTAAAAACCCTCAGCAGCCTCAAAGACCCCTGATGAAATTTCCAAGTGATGCTGT
#TCTTAACAAAATTGTAGAGAACTCCTCAAAGAACATTTTTATTTGTAACGAAGAGGATGA
#GAAAATAGAAGAAGAGAAAAATAAACAGCTTAATGATAAGCTGAGGCTTTCCGAAGTCCC
#ATTACAAGAAGAAAGTAAAGATCAGATCATCTCTGGAGAAAATCTTGACCAAATAGTCTC
#ATTTTCCTCAAGAAAACCACTAGCCAAAAAGATTGTTCCCAATTCAGAGTCTCCTTTACG
#AAAAAGTCGAACTAGAGAGTCCCTCAAAAAGAGGAAGGCTCTTCGAAAGAGGATCATCTT
#CACAAGCCAGGAAGGAAGTGAGAGGCTTGAATCTGATCAAAATCTCCTCCA"""
	# comp16092
	sequence = """ATGATCTTCACTTTCCTGGCCATTACACTGATAATCCGTGGGTTCACTTTAATACCGATT
ATAACTGATGAAGACTTTGAAGAGGAGCCTACCAAATGTCAAGAGGCTGCTTTTGAGTAC
ATGCCACTTACGACTTTTCTAATCGCAGTGATGGTGAACACTTATAGGTTTTATATAATC
CGGAGATCTCCGAAATTATGAACGAAGATAGCTTGACACATGATCTTATGGATCTATATA
TTATTCATCCCATTTATCAACATTCTATTCTGAGTGGACCACGGTATTAGCCTCTTAATC
AGAGTTGTGAAAATTACTGCCATTGTGATCTACTGCTTTATGAATTTTGCATGAGCTGTG
ATAAACCTGTACATAATGACCTGCCTCTTCCGGAAGATGAATATGCTTTTATTCACTTTG
ATATGCATATGTGTACTGAGAGTAATATCGATAACAGTTTTTACCTTTGTTAAAAATGAT
AATTTTGGAGCCTTGAATAATTTTGAGTTCATGCTTATACTGAGCTTCACAGAGCTTATC
CCTGTTATAATGATTATCACCTCATTTTTAAGATATATGAATAATCAACCTGGAGCAGAT
TTTATTGAAGAAAGTGCGAACTATGCTCCTCATAAGAAAGTATCCTTACAAAGGACTGAT
ACAGTGACTTCTCAAAATTCAGATGATGTTAATAAAACCAGCGACCAAACACAGGATCTA
ATAAAATCAATACTAGAGGTGGAAGAAAAGAGACACTCCAGAGCCGAAGAATTTGGAAAT
AAGAGGAACGGCCTTTCCAAAACAGTCGACGATCCTGCTGCCTTTGAGCAAGCATCTGCA
TATCGTCCTCCAAATGACAATAGCATTCCAGATAGCCCCAGGATGAGGATTATAGCTCGT
TCCAAGACGAAACAGTTTGAAGAGCAGAAAGAAGAAAGCAAGGAAGATCGCTTTCAAAAG
GTAGAATATGCCTTGAAATACATGATAAGTGACAATGCCATCAGTAGGTCCTCATCAGTG
ATGACAAAGAAATCTGTGGAGGATTCTGACAGCTCAAGAGTATCTCAAAGTCAGTCTCTT
TAAGACTGAGTTAATTGGTAAATATTTTTAAATT"""
	sequence = sequence.replace( "\n", "" )
	s = Sequence( sequence )
	s.truncate( verbose=True )
	no_of_leaves = s.count_leaves()
	if no_of_leaves > 1000:
		print >> sys.stderr, "Complex tree with %s leaves... omitting." % \
			no_of_leaves
		sys.exit( 0 )
	# set the transition matrix
	s.set_transition_matrix( TM )
	# get the stops in all frames
	# sanitise the list of stop sites
	s.build_tree()
	s.get_frameshift_signals()
	s.estimate_likelihood()
	s.estimate_frameshift_likelihood()
	# get the most likely sequence
	s.get_most_likely_frameshift()
	ML = s.most_likely_frameshift
	if ML is not None:
		# compute the GC content of the most likely sequence
		print ML.estimate_GC_content()
		# get the gradient of the most likely sequence, m_O
		#print ML.differential_graded_likelihood
		print "m_O\t=\t%.4f" % ML.estimate_gradient()
		# for each stop site in the most likely sequence
		print ML.path
		print ML.fragments
		print ML.signals
		no_of_codons = 5
		for i in xrange( len( ML.fragments )-1 ):
			#print ML.frameshifted_sequence[pos:pos+3]
			# get a window about the stop
			before_len = len( ML.fragments[i] )
			before = ML.fragments[i][max([ -1*no_of_codons*3, -before_len ]):]
			after_len = len( ML.fragments[i+1] )
			after = ML.fragments[i+1][:min([ no_of_codons*3, after_len ])]
			print before
			print after
			b = Sequence( before )
			b.set_transition_matrix( TM )
			b.estimate_likelihood()
			# compute gradient up to the stop, m_I
			m_I = b.estimate_gradient()
			print "m_I\t=\t%.4f" % m_I
			a = Sequence( before + after )
			a.set_transition_matrix( TM )
			a.estimate_likelihood()
			# compute the gradient across the stop, m
			m = a.estimate_gradient()
			print "m\t=\t%.4f" % m
			# compute rho = m/m_I
			rho = m/m_I
			print "rho\t=\t%.4f" % rho
			print

if __name__ == "__main__":
	main()