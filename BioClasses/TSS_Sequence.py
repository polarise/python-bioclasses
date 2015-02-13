# -*- encoding: utf-8 -*-
import sys
import scipy

class TSS_Sequence( object ):
	def __init__( self, seqname, pos ):
		self.seqname = seqname
		self.pos = pos
		self.flank_right = 30
		self.flank_left = 30
		self.sequence = str()
		self.length = self.flank_right + self.flank_left + 1
		
		# pswm-specific attributes
		# these are set when the pswm method is invoked
		self.PSWM_width = None
		self.TOP_scores = list()
		self.TOP_candidate = None
		self.TOP_candidate_pos = None
		self.TOP_score = None
	
	def __repr__( self ):
		"""
		TSS on chr3 @ pos 02830982304
		Flanks: [30/30]
		TSS region: chr3:028320-20380283
		Sequence: CGATGCTAGCATGTAGCTGTACTG
		PSWM width: 11
		TOP score(s): [-35,-35-3,-3,-53,-35,-3]
		TOP candidate: -3
		TOP position: 3
		TOP sequence: ATCGACTAGC"""
		if self.sequence == "":
			TOP_sequence = ""
		else:
			TOP_sequence = self.sequence[self.TOP_candidate_pos:self.TOP_candidate_pos + self.PSWM_width]
			
		string = """\
TSS on %s @ pos %s
-----------------------------
Flanks:        [%s/%s]
TSS region:    %s
Sequence:      %s
PSWM width:    %s
TOP score(s):  %s
TOP candidate: %s
TOP position:  %s
TOP score:     %s
TOP sequence:  %s
""" % ( self.seqname, self.pos, self.flank_left, self.flank_right, self.get_tss_region_str(), self.sequence, self.PSWM_width, str( map( lambda e: round( e, 2 ), self.TOP_scores )), self.TOP_candidate, self.TOP_candidate_pos, self.TOP_score, TOP_sequence )
		return string
	
	def set_flank_right( self, flank ):
		self.flank_right = flank
	
	def set_flank_left( self, flank ):
		self.flank_left = flank
	
	def get_tss_region_str( self, zero_based=False ):
		if zero_based:
			return "%s:%s-%s" % ( self.seqname, int( self.pos ) - 1 - self.flank_left, int( self.pos ) - 1 + self.flank_right )
		elif not zero_based:
			return "%s:%s-%s" % ( self.seqname, int( self.pos ) - self.flank_left, int( self.pos ) + self.flank_right )
	
	def get_sequence( self, fastafile ): # indexed fasta
		# TODO: raise exception if fastafile is not indexed
		tss_region = self.get_tss_region_str()
		sequence = fastafile.fetch( region=tss_region )
		self.sequence = sequence
		return sequence
	
	def compute_TOP_score( self, pswm_obj ):
		# set the width of the PSWM_width attribute
		self.PSWM_width = pswm_obj.cols
		
		# internal function
		def convert_sequence_to_matrix( sequence ):
			ref = {'A': [ 1, 0, 0, 0 ],
						'C': [ 0, 1, 0, 0 ],
						'G': [ 0, 0, 1, 0 ],
						'T': [ 0, 0, 0, 1 ]}
			matrix = list()
			for s in sequence:
				matrix.append( ref[s] )
			matrix = scipy.array( matrix )
			return matrix.T

		# is the sequence long enough
		TOP_scores = list()
		if pswm_obj.cols > self.length: # can be resolved later
			#raise ValueError( "Sequence length shorter than PSWM col width. Use shorter sequence." )
			TOP_scores.append( -99 )		
		else:
			for i in xrange( 0, self.length - pswm_obj.cols + 1, 1 ):
				subseq = self.sequence[i:pswm_obj.cols + i]
				subseq_matrix = convert_sequence_to_matrix( subseq )
				score = sum( sum( pswm_obj * subseq_matrix ))
				TOP_scores.append( score )
		
		self.TOP_scores = TOP_scores
		
		# get the highest TOP score and region
		i = self.TOP_scores.index( max( self.TOP_scores ))
		TOP_candidate = self.sequence[i:self.PSWM_width + i]
		self.TOP_candidate = TOP_candidate
		self.TOP_candidate_pos = i
		self.TOP_score = TOP_scores[i]
		
		return TOP_scores, TOP_candidate, i, TOP_scores[i]