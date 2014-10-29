# -*- encoding: utf-8 -*-
from __future__ import division
import sys

class AminoAcid( object ):
	def __init__( self, name, lsymbol, ssymbol, codons ):
		self.name = name # full name e.g. 'Phelinalanine'
		self.lsymbol = lsymbol # long symbol e.g. 'Phe'
		self.ssymbol = ssymbol # short symbol e.g. 'F'
		self.codons = codons.split( "," ) # list
		self.codon_counts = dict()
		self.norm_counts = dict()
		self.wij = dict()
	
	def __repr__( self ):
		count_str = map( str, [ (c,self.codon_counts[c]) for c in self.codon_counts ])
		norm_str = map( str, [ (n,self.norm_counts[n]) for n in self.norm_counts ])
		wij_str = map( str, [ (w,self.wij[w]) for w in self.wij ] )
		return """\
Amino acid: %s (%s/%s)
Codons    : %s
Norm'd    : %s
Counts    : %s
w_ij      : %s""" % ( self.name, self.lsymbol, self.ssymbol, \
	",".join( self.codons ), ",".join( norm_str ), ",".join( count_str ), \
		",".join( wij_str ))
	
	def __str__( self ):
		return self.__repr__()
	
	def set_count( self, codon, codon_count ):
		self.codon_counts[codon] = codon_count
		
	def count_codon( self, codon ):
		if codon not in self.codon_counts:
			self.codon_counts[codon] = 1
		else:
			self.codon_counts[codon] += 1
	
	def normalise( self ):
		total_count = sum( self.codon_counts.values() )
		for c in self.codon_counts:
			self.norm_counts[c] = self.codon_counts[c]/total_count
		
		max_norm_count = max( self.norm_counts.values() )
		for n in self.norm_counts:
			self.wij[n] = self.norm_counts[n]/max_norm_count		
	
	def get_norm_count( self, codon ):
		return self.norm_counts[codon]