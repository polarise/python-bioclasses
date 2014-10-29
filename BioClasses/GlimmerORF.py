# -*- encoding: utf-8 -*-
from __future__ import division
import sys

"""
>comp9463_c0_seq1
orf00001        1      900  +1     3.00
"""

def reverse_comp( sequence ):
	rev_sequence = list()
	for s in sequence:
		if s == "A":
			rev_sequence.append( "T" )
		elif s == "C":
			rev_sequence.append( "G" )
		elif s == "G":
			rev_sequence.append( "C" )
		elif s == "T":
			rev_sequence.append( "A" )
	
	return "".join( rev_sequence[::-1] )

class ORF( object ):
	def __init__( self, ID, start, end, frame, score ):
		self.ID = ID
		self.start = int( start )
		self.end = int( end )
		self.frame = frame
		self.score = float( score )
		if self.start <= self.end:
			self.strand = "+"
		else:
			self.strand = "-"
		self.sequence = None
	
	def __repr__( self ):
		return "\t".join( map( str, [ self.ID, self.start, self.end, \
			self.frame, self.score, self.strand ]))
	
	def print_row( self ):
		return "|".join( map( str, [ self.ID, self.start, self.end, \
			self.frame, self.score, self.strand ]))

################################################################################

class GlimmerORF( object ):
	def __init__( self, name ):
		self.name = name
		self.orfs = dict()
	
	def add_orf( self, orf ):
		orf_list = [ o for o in orf.split( " " ) if o != "" ]
		O = ORF( *orf_list )
		self.orfs[O.ID] = O
	
	def extract_orf( self, sequence ):
		# all coordinates are relative to the first ATG
		atg_pos = sequence.find( "ATG" )
		if atg_pos < 0:
			print >> sys.stderr, "Missing first ATG"
		else:
			sequence2 = sequence[atg_pos:]
			for o in self.orfs:
				if self.orfs[o].strand == "+":
					sequence3 = sequence2[self.orfs[o].start-1:self.orfs[o].end]
					if len( sequence3 ) > 6:
						self.orfs[o].sequence = sequence3
	
	def __repr__( self ):
		out_str = "ORFs in %s:\n" % self.name
		for o in self.orfs:
			if self.orfs[o].sequence is not None:
				out_str += self.orfs[o].__repr__() + "\n"
				out_str += self.orfs[o].sequence + "\n"
		return out_str
	
	def print_as_fasta( self ):
		for o in self.orfs:
			if self.orfs[o].sequence is not None:
				print ">%s" % self.orfs[o].print_row()
				print self.orfs[o].sequence
	
	def print_terminal_sequence( self ):
		for o in self.orfs:
			if self.orfs[o].sequence is not None:
				print "\t".join([ self.orfs[o].print_row(), self.orfs[o].sequence[:6], \
					self.orfs[o].sequence[-6:] ])