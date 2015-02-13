# -*- encoding: utf-8 -*-
import sys
from Exon import *
from UTR import *
from TSS_Sequence import *

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
		self.UTR = list() # exons by the way
		self.start_codon = None
		self.stop_codon = None
		if self.strand == "+":
			self.TSS_Sequence = TSS_Sequence( self.seqname, self.start )
		elif self.strand == "-":
			self.TSS_Sequence = TSS_Sequence( self.seqname, self.end )
	
	#=============================================================================
	
	def is_complete( self ):
		if self.start_codon is not None and self.stop_codon is not None:
			return True
		else:
			return False
	
	#=============================================================================
	
	def is_protein_coding( self ):
		if self.source == "protein_coding":
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
		self.CDS.sort() # invoke the __cmp__() method of each exon
		if self.strand == "+":
			if self.start_codon is not None:
				cds_start = int( self.start_codon )
			elif len( self.CDS ) > 0:
				cds_start = int( self.CDS[0].start )
			else:
				cds_start = int( self.start )
                        
			if self.stop_codon is not None:
				cds_end = int ( self.stop_codon )
			elif len( self.CDS ) > 0:
				cds_end = int( self.CDS[-1].end )
			else:
				cds_end = int( self.end )
			
			if zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, cds_start - 1, cds_end + 1 ] )) 
				#return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.start_codon ) - 1, int( self.stop_codon ) + 1 ]))
			elif not zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, cds_start, cds_end + 2 ] ))	
				#return "%s:%s-%s" % ( self.seqname, self.start_codon , str( int( self.stop_codon ) + 2 ))
		elif self.strand == "-":
			if self.start_codon is not None:
				cds_start = int( self.start_codon )
			elif len( self.CDS ) > 0:
				cds_start = int( self.CDS[-1].end )
			else:
				cds_start = int( self.end )

			if self.stop_codon is not None:
				cds_end = int( self.stop_codon )
			elif len( self.CDS ) > 0:
				cds_end = int( self.CDS[-1].start )
			else:
				cds_end = int( self.start )

			if zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, cds_end - 3, cds_start - 1 ]))
				#return "%s:%s-%s" % tuple( map( str, [ self.seqname, int( self.stop_codon ) - 3, int( self.start_codon ) - 1 ]))
			elif not zero_based:
				return "%s:%s-%s" % tuple( map( str, [ self.seqname, cds_end - 2, cds_start ]))	
				#return "%s:%s-%s" % ( self.seqname, str( int( self.stop_codon ) - 2 ), self.start_codon )
	
	#=============================================================================
	
	def get_cds_unspliced_length( self ):
		length = 0
		if self.start_codon is not None and self.stop_codon is not None:
			if self.strand == "+":
				length = int( self.stop_codon ) + 2 - int( self.start_codon ) + 1
			elif self.strand == "-":
				length = int( self.start_codon ) + 2 - int( self.stop_codon ) + 1
		return length
	
	#=============================================================================
	def get_cds_length( self ):
		length = 0
		for E in self.CDS:
			length += E.length() #( int( E.end ) - int( E.start ) + 1 )
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
	def length_of_unspliced_5utr( self ):
		unspliced_length = 0
		if self.start_codon is not None:
			if self.strand == "+":
				unspliced_length = int( self.start_codon ) - int( self.start )
			elif self.strand == "-":
				unspliced_length = int( self.end ) - int( self.start_codon ) - 2
		return unspliced_length
	#=============================================================================
	def length_of_unspliced_3utr( self ):
		unspliced_length = 0
		if self.stop_codon is not None:
			if self.strand == "+":
				unspliced_length = int( self.end ) - int( self.stop_codon ) - 2
			elif self.strand == "-":
				unspliced_length = int( self.stop_codon ) - int( self.start )
		return unspliced_length
	#=============================================================================
	def length_of_spliced_5utr( self ):
		"""
		return the spliced length of the 5'UTR
		"""
		spliced_length = 0
		
		if self.start_codon is not None:
			if self.strand == "+":
				"""
				we are only interested in exonic regions up to the start codon
				there are therefore two types of exons:
				(i)		those that wholly lie on the left of the start codon
				(ii)	those that cover the start codon
				
				This algorithm relies on exon_number as the index.
				exon_number increases in the direction of transcription.
				But exon_number is NOT sorted!!!
				"""
				for exon_number in xrange( 1, len( self.exons )):
					E = self.exons[exon_number]
					# case (i)
					if int( E.start ) < int( self.start_codon ) and int( E.end ) < int( self.start_codon ):
						spliced_length += E.length()
					# case (ii)
					elif int( E.start ) < int( self.start_codon ) and int( E.end ) > int( self.start_codon ):
						spliced_length += int( self.start_codon ) - int( E.start )
					# break on all other cases
					else:
						break
					# spliced_length will be zero (0)
			elif self.strand == "-":
				for exon_number in xrange( 1, len( self.exons )):
					E = self.exons[exon_number]
					# case (i)
					if int( E.start ) > int( self.start_codon ) + 2: # the 3' end doesn't matter in minus strand genes
						spliced_length += E.length()
					# case (ii)
					elif int( E.start ) < int( self.start_codon ) + 2 and int( E.end ) > int( self.start_codon ) + 2:
						spliced_length += int( E.end ) - int( self.start_codon ) - 2
					# break on all other cases
					else:
						break
					# spliced_length will be zero (0)					
		
		return spliced_length
	
	#=============================================================================
	
	def length_of_spliced_3utr( self ):
		"""
		return th spliced length of the 3'UTR
		"""
		spliced_length = 0
		
		if self.stop_codon is not None:
			if self.strand == "+":
				for exon_number in xrange( len( self.exons ), 0, -1 ):
					E = self.exons[exon_number]
					# case (i)
					if int( E.start ) > int( self.stop_codon ) + 2:
						spliced_length += E.length()
					# case (ii)
					elif int( E.start ) < int( self.stop_codon ) + 2 and int( E.end ) > int( self.stop_codon ):
						spliced_length += int( E.end ) - int( self.stop_codon ) - 2
					else:
						break
			elif self.strand == "-":
				for exon_number in xrange( len( self.exons ), 0, -1 ):
					E = self.exons[exon_number]
					# case (i)
					if int( E.start ) < int( self.stop_codon ) and int( E.end ) < int( self.stop_codon ):
						spliced_length += E.length()
					# case (ii)
					elif int( E.start ) < int( self.stop_codon ) and int( E.end ) > int( self.stop_codon ):
						spliced_length += int( self.stop_codon ) - int( E.start )
					else:
						break
		
		return spliced_length
	
	#=============================================================================
	
	def process_exon( self, record ):
		self.exons[int( record.group_dict['exon_number'] )] = Exon( record )
#		self.exons[record.group_dict['exon_id']]  = Exon( record )
	
	#=============================================================================
	
	def process_CDS( self, record ):
		self.CDS.append( Exon( record ))
	
	#=============================================================================
		
	def process_UTR( self, record ):
		self.UTR.append( UTR( record ))
	
	#=============================================================================

	def process_5UTR( self, record ):
		self.UTR.append( UTR( record, terminus="5" ))
	
	#=============================================================================
	
	def process_3UTR( self, record ):
		self.UTR.append( UTR( record, terminus="3" ))
	
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
	
	#=============================================================================
	
	def designate_UTRs( self, zero_based=False ):
		# compare coordinates of each UTR exon to the cds start/end
		cds_region_str = self.cds_region_str( zero_based=zero_based )
		cds_seqname,cds_start_end = cds_region_str.split( ":" )
		cds_start_inferred,cds_end_inferred = map( int, cds_start_end.split( "-" ))
		for utr in self.UTR:
			if utr.terminus is not None:
				continue # don't waste time
			if self.strand == "+":
				if self.start_codon is not None:
					cds_start = int( self.start_codon )
				else:
					cds_start = cds_start_inferred
				
				if int( utr.end ) <= cds_start:
					utr.terminus = "5"
				
				if self.stop_codon is not None:
					cds_end = int( self.stop_codon ) + 2
				else:
					cds_end = cds_end_inferred
				
				if int( utr.start ) >= cds_end:
						utr.terminus = "3"
			elif self.strand == "-":
				if self.start_codon is not None:
					cds_start = int( self.start_codon )
				else:
					cds_start = cds_end_inferred
			
				if int( utr.start ) >= cds_start:
						utr.terminus = "5"
			
				if self.stop_codon is not None:
					cds_end = int( self.stop_codon ) - 2
				else:
					cds_end = cds_start_inferred

				if int( utr.end ) <= cds_end:
					utr.terminus = "3"
			
