# -*- encoding: utf-8 -*-
import sys

class GTFRecord( object ):
	def __init__( self, seqname, source, feature, start, end, score, strand, frame, group ):
		self.seqname = seqname
		self.source = source
		self.feature = feature
		self.start = start
		self.end = end
		self.score = score
		self.strand = strand
		self.frame = frame # coding frame! ;-)
		self.group = group
		self.group_dict = self.process_group( group )
	
	def process_group( self, group ):
		group_dict = dict()
		group_list = group.split( " " )
		for i in xrange( 0, len( group_list[:-1] ), 2 ):
			key = group_list[i]
			value = group_list[i+1].strip( ";" ).strip( "\"" )
			group_dict[key] = value
		
		return group_dict
	
	def region_str( self ):
		return "%s:%s-%s" % ( self.seqname, self.start, self.end )
	
	def __repr__( self ):
		return "\t".join([ self.group_dict['gene_id'], self.group_dict['transcript_id'], self.start, self.end ])
		
		
