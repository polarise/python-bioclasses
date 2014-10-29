# -*- encoding: utf-8 -*-
import sys
from Exon import *

class Transcript( object ):
	def __init__( self, record ):
		self.transcript_id = record.group_dict['transcript_id']
		self.source = record.source # protein_coding|...|miRNA| etc
		self.seqname = record.seqname
		self.start = record.start # 1-based
		self.end = record.end			# 1-based
		self.strand = record.strand
		self.exons = dict()
		self.CDS = list() # exons by the way
		self.UTR = dict() # exons by the way
		self.start_codon = None
		self.stop_codon = None
	
	#=============================================================================
	
	def is_complete( self ):
		if self.start_codon is not None and self.stop_codon is not None:
			return True
		else:
			return False
	
	#=============================================================================
	
	def region_str( self, zero_based=False ):
		if zero_based:
			return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.start ) - 1, int( self.end ) - 1 ]))
		elif not zero_based:
			return "%s:%s-%s" % ( self.seqname, self.start, self.end )
	
	#=============================================================================
	
	def cds_region_str( self, zero_based=False ):
		if self.strand == "+":
			if zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.start_codon ) - 1, int( self.stop_codon ) + 1 ]))
			elif not zero_based:
				return "%s:%s-%s" % ( self.seqname, self.start_codon , str( int( self.stop_codon ) + 2 ))
		elif self.strand == "-":
			if zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.stop_codon ) - 3, int( self.start_codon ) - 1 ]))
			elif not zero_based:
				return "%s:%s-%s" % ( self.seqname, str( int( self.stop_codon ) - 2 ), self.start_codon )
	
	#=============================================================================
	
	def get_cds_length( self ):
		length = 0
		for E in self.CDS:
			length += ( int( E.end ) - int( E.start ) + 1 )
		return length			
	
	#=============================================================================
	
	def get_full_length( self ):
		length = 0
		for exon_id,E in self.exons.iteritems():
			length += ( int( E.end ) - int( E.start ) + 1 )
		return length
	
	#=============================================================================
	
	def utr_region_str( self, zero_based=False ):
		if self.strand == "+":
			if zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.start ) - 1, int( self.start_codon ) - 2 ]))
			elif not zero_based:
				return "%s:%s-%s" % ( self.seqname, self.start, str( int( self.start_codon ) - 1 ))
		elif self.strand == "-":
			if zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.start_codon ), int( self.end ) - 1] ))
			elif not zero_based:
				return "%s:%s-%s" % ( self.seqname, str( int( self.start_codon ) + 1 ), self.end )
	
	#=============================================================================
	
	def process_exon( self, record ):
		self.exons[int( record.group_dict['exon_number'] )] = Exon( record )
#		self.exons[record.group_dict['exon_id']]  = Exon( record )
	
	#=============================================================================
	
	def process_CDS( self, record ):
		self.CDS.append( Exon( record ))
	
	#=============================================================================
		
	def process_UTR( self, record ):
		pass
	
	#=============================================================================	
	
	def process_start_codon( self, record ):
		self.start_codon = record.start
	
	#=============================================================================
		
	def process_stop_codon( self, record ):
		self.stop_codon = record.start
	
	#=============================================================================
	
	def __repr__( self ):
#		return "%s [%s]" % ( self.transcript_id, self.region_str())
		if self.strand == "+":
			return "%s\t[%s(%s):(%s)%s]" % ( self.transcript_id, self.start, self.start_codon, self.stop_codon, self.end )
		elif self.strand == "-":
			return "%s\t[%s(%s):(%s)%s]" % ( self.transcript_id, self.end, self.start_codon, self.stop_codon, self.start )
	
	#=============================================================================
	
	def __eq__( self, T ):
		"""
		we define equality based on having the same CDS
		"""
		if self.start_codon == T.start_codon and self.stop_codon == T.stop_codon:
			return True
		else:
			return False
		
