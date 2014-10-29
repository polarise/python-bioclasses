#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from GeneticCode import *
from Sequence import *

def main( gc ):	
	G = GeneticCode( gc )
	# G.build_CAI_table( "/Users/paulkorir/Dropbox/Euplotes/FrameshiftPredictionData/E.crassus_CDS.fasta" )
	G.build_CAI_table( "/home/paul/bioinf/Resources/H_sapiens/H_sapiens_Ens75_CDS.fasta" )
	#G.build_CAI_table( "a_fasta.fa" )
	G.write_CAI_table( "CAI_tables/homo_CAI_table.txt" )
	
	sys.exit( 0 )
	
	# G.write_CAI_table( "euplotid_CAI_table.txt" )
	#G.write_CAI_table( "test_CAI_table.txt" )
	G.read_CAI_table( "euplotid_CAI_table.txt" )
	#print G
	s = Sequence( "TAGAGATACACTGACTTACTTTCAAATACTATAAAACGGAATAGCCTAAGAATGAAATAAAGTAAAACATGACCATCAGGAGAAAGTTGAACAACTAGAGAGGGAGAATATTAAGCTTTATGCCCAATTAAAAAAGCTTGCAAAAAGTGAAAGAAATCTAATGAAGAAACTAGACGAAAGAGACCGGGAGATAACCAATCTAAAAGATACAAACATGAGGTTCAATTACAAACTCAATAGAGCACTCTATGCTAATGAAGAGCTGCAAAATAAAGTAACTGAATCTGACTACAAACTTCAACAAAAAAGAGATGAATTTATGAAAGACATAGAGCAAACTAACCAAATCC" )
	#s = Sequence( "TCAAACCGAGACTTACTAAAGTTGATCATCATAAGACTC" )
	s.set_genetic_code( G )
	s.estimate_CAI()
	# s.as_codons = True
	print s
	# print s.CAI_score
	s.truncate()
	print s
	s.build_tree()
	# print s.tree
	for i in xrange( 3 ):
		print s.binary_frame( i, "" )
	s.estimate_frameshift_CAI()
	for fs in s.frameshift_sequences:
		print s.frameshift_sequences[fs].repr_as_row()
	
if __name__ == "__main__":
	try:
		gc = sys.argv[1]
	except IndexError:
		print >> sys.stderr, "Error: missing genetic code..."
		sys.exit( 1 )
	main( gc )