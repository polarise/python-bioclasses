# -*- encoding: utf-8 -*-
from __future__ import division
import sys

class UTR( object ):
	def __init__( self, record, terminus=None ):
		self.source = record.source
		self.seqname = record.seqname
		self.start = record.start
		self.end = record.end
		self.strand = record.strand
		self.terminus = terminus # '5' for 5' UTR or '3' for 3' UTR
		self.sequence = ""
		self.length = 0
		self.GC_percent = 0
		self.AT_percent = 0
	
	def __repr__( self ):
		string = ">%s:%s-%s:%s %s'UTR-exon|%sbp|%%GC=%s|%%AT=%s\n" % ( self.seqname, self.start, self.end, self.strand, self.terminus, self.length, round( self.GC_percent, 2 ), round( self.AT_percent, 2 ))
		string += "%s" % self.sequence
		return string
	
	def region_str( self, zero_based=False ):
		if zero_based:
			return "%s:%s-%s" % ( self.seqname, str( int( self.start ) + 1 ), str( int( self.end ) + 1 ))
		else:
			return "%s:%s-%s" % ( self.seqname, self.start, self.end )
		
	def length( self ):
		return int( self.end ) - int( self.start ) + 1
	
	def get_sequence( self, fastafile ):
		# TODO: raise exception if fastafile is not indexed
		region_str = self.region_str()
		sequence = fastafile.fetch( region=region_str )
		self.sequence = sequence
		self.length = len( sequence )
		self.GC_percent = ( sequence.upper().count( 'G' ) + sequence.upper().count( 'C' ))/len( sequence )*100
		self.AT_percent = ( sequence.upper().count( 'A' ) + sequence.upper().count( 'T' ))/len( sequence )*100
		return sequence