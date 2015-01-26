# -*- encoding: utf-8 -*-
import sys
import itertools
from Transcript import *
from Exon import *

class Gene( object ):
	def __init__( self, record ):
		self.gene_id = record.group_dict['gene_id']
		self.gene_name = record.group_dict['gene_name']
		self.seqname = record.seqname
		self.source = record.source
		self.start = record.start
		self.end = record.end
		self.strand = record.strand
		self.transcripts = dict()
	
	#=============================================================================
	
	def is_protein_coding( self ):
		if self.source == "protein_coding":
			return True
		else:
			return False
	
	#=============================================================================
	
	def region_str( self, zero_based=False ):
		if zero_based:
			return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.start ) + 1, int( self.end ) + 1 ]))
		elif not zero_based:
			return "%s:%s-%s" % ( self.seqname, self.start, self.end )		
	
	#=============================================================================
	
	def process_record( self, record ):
		if record.feature == "transcript":
			self.transcripts[record.group_dict['transcript_id']] = Transcript( record )
		elif record.feature == "exon":
			self.transcripts[record.group_dict['transcript_id']].process_exon( record )
		elif record.feature == "CDS":
			self.transcripts[record.group_dict['transcript_id']].process_CDS( record )
		elif record.feature == "UTR":
			self.transcripts[record.group_dict['transcript_id']].process_UTR( record )	
		elif record.feature == "start_codon":
					self.transcripts[record.group_dict['transcript_id']].process_start_codon( record )	
		elif record.feature == "stop_codon":
					self.transcripts[record.group_dict['transcript_id']].process_stop_codon( record )
	
	#=============================================================================
	
	def __repr__( self ):
#		return "%s [%s]" % ( self.gene_id, self.region_str())
		return self.gene_id + " [" + ":".join([ self.seqname, self.start, self.end, self.strand ]) + "]"
	
	#=============================================================================
	
	def get_protein_coding( self ):
		PCT = list()
		for transcript_id,T in self.transcripts.iteritems():
			if T.source == "protein_coding":
				PCT.append( T )
		
		return PCT
	
	#=============================================================================
	
	def get_complete_protein_coding( self ):
		PCT = self.get_protein_coding()
		complete_PCT = list()
		for P in PCT:
			if P.is_complete():
				complete_PCT.append( P )
		
		return complete_PCT
	
	def get_equal_cds( self, best=False ):
		"""
		if best=True then return the set of txs in the modal class of txs with equal CDS
		otherwise return all that share
		"""
		CPCT = self.get_complete_protein_coding() # complete protein coding transcript objects
		
		# we get the distribution of complete transcripts by (start_codon, stop_codon)
		start_stop_distr = dict()
		for P in CPCT:
			if P.is_complete():
				start_stop = P.start_codon,P.stop_codon
				if start_stop not in start_stop_distr:
					start_stop_distr[start_stop] = 1
				else:
					start_stop_distr[start_stop] += 1
			else:
				pass
#				print >> sys.stderr, "Incomplete transcript: %s" % P
		
		# pick the start_stop with the most transcripts
		equal_cds_transcripts = set()
		if len( start_stop_distr ) > 0:
			if best: # return the modal class
				start_stop_best = start_stop_distr.keys()[0] # initialise the best as the first
				start_stop_max = start_stop_distr[start_stop_distr.keys()[0]] # guess that first is max
			
				for start_stop,start_stop_count in start_stop_distr.iteritems():
					if start_stop_count > start_stop_max:
						start_stop_best = start_stop
					
				# now get the transcripts corresponding
				best_start, best_stop = start_stop
				for P in CPCT: # transcript objects
					if P.start_codon == best_start and P.stop_codon == best_stop:
						equal_cds_transcripts.add( P )
			else: # return all classes except those with only one transcript
				for P in CPCT:
					try:
						count = start_stop_distr[P.start_codon,P.stop_codon]
					except KeyError:
						count = 0
					if count > 1:
						equal_cds_transcripts.add( P )
		
		if len( equal_cds_transcripts ) > 1:
			return list( equal_cds_transcripts )
		else:
			return None
	
	#=============================================================================
	
	def get_overall_UTR_region( self, transcripts="all", as_region_str=True ):
		"""
		transcripts = all|equal_cds
		"""
		def process_region_str( region_str ):
			seqname, start_end = region_str.split( ":" )
			start, end = start_end.split( "-" )
			return seqname, start, end
		
		if transcripts == "all":
			candidate_transcripts = self.get_complete_protein_coding()
		elif transcripts == "equal_cds":
			candidate_transcripts = self.get_equal_cds()
		
		if candidate_transcripts is None:
			return None
		else:
			overall_utr_start = list()
			overall_utr_end = list()
			utr_starts = list()
			utr_ends = list()
			for T in candidate_transcripts:
				region_str = T.utr_region_str()
				seqname, start, end = process_region_str( region_str )
				start_i = int( start )
				end_i = int( end )
				utr_starts.append( start_i )
				utr_ends.append( end_i )
		
			try:
				overall_utr_start = min( utr_starts )
			except ValueError:
				overall_utr_start = None
				print >> sys.stderr, "Warning: empty utr_starts for %s" % self.gene_id
		
			try:
				overall_utr_end = max( utr_ends )
			except ValueError:
				overall_utr_end = None
				print >> sys.stderr, "Warning: empty utr_ends for %s" % self.gene_id
		
			# construct the final region string
			if as_region_str:
				return "%s:%s-%s" % ( self.seqname, overall_utr_start, overall_utr_end )
			elif not as_region_str:
				return self.seqname, str( overall_utr_start ), str( overall_utr_end )
	#=============================================================================
	
	def post_processing( self ):
		for transcript_id,T in self.transcripts.iteritems():
			# designate UTRs as either 5' or 3'
			T.designate_UTRs()
