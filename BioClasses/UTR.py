# -*- encoding: utf-8 -*-
import sys

class UTR( object ):
	def __init__( self, record, terminus=None ):
		self.source = record.source
		self.seqname = record.seqname
		self.start = record.start
		self.end = record.end
		self.strand = record.strand
		self.terminus = terminus # '5' for 5' UTR or '3' for 3' UTR
	
	def region_str( self, zero_based=False ):
		if zero_based:
			return "%s:%s-%s" % ( self.seqname, str( int( self.start ) + 1 ), str( int( self.end ) + 1 ))
		else:
			return "%s:%s-%s" % ( self.seqname, self.start, self.end )
		
	def length( self ):
		return int( self.end ) - int( self.start ) + 1