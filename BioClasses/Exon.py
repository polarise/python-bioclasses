# -*- encoding: utf-8 -*-
import sys

class Exon( object ):
	def __init__( self, record ):
		if record.source == "exon":
			self.exon_id = record.group_dict['exon_id']
		elif record.source == "CDS":
			self.exon_id = record.group_dict['protein_id']
		elif record.source == "start_codon" or record.source == "stop_codon":
			self.exon_id = record.group_dict['ccds_id']
		self.source = record.source
		self.exon_number = record.group_dict['exon_number']
		self.seqname = record.seqname
		self.start = record.start
		self.end = record.end
		self.strand = record.strand
	
	def region_str( self, zero_based=False ):
		if zero_based:
			return "%s:%s-%s" % ( self.seqname, str( int( self.start ) + 1 ), str( int( self.end ) + 1 ))
		else:
			return "%s:%s-%s" % ( self.seqname, self.start, self.end )
	
	def length( self ):
		return int( self.end ) - int( self.start )
		
