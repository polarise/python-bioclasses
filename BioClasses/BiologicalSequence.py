# -*- encoding: utf-8 -*-
from Sequence import *

class BiologicalSequence( Sequence ):
	def __init__( self, sequence, bases='ACGT', starts=[ 'ATG' ], stops=[ 'TAA', 'TAG' ] ):
		super( BiologicalSequence, self ).__init__( sequence=sequence, bases=bases, starts=starts, stops=stops )
		self.length = len( sequence )
	
		# must have a valid binary codon matrix
		self.set_binary_codon_matrix( weighted_start=True )
		
		self.frameshifted_sequence = None
		self.path = None
		self.fragments = list()
		self.frameshift_signals = list()
		
	#*****************************************************************************
	
	def info( self, comment="" ):
		return "Unknown biological sequence of length %s bases." % self.length
		
	#*****************************************************************************
	
	def frameshift_from_path( self, path ):
		"""
		"""
		# first get all frame of self.sequence
		sequence_in_frames = dict()
		for i in xrange( 3 ):
			sequence_in_frames[i] = self.sequence[i:]
		
		#for f in sequence_in_frames:
			#print f, ":", sequence_in_frames[f]
		#print "    " + "         |"*(( self.length )//10 )
		#print
		
		self.frameshifted_sequence = ""
		self.fragments = list()
		self.frameshift_signals = list()
		i = 0
		f_i = 0
		for f,j in path:
			self.frameshifted_sequence += self.sequence[i+(f-f_i):j]
			self.fragments.append( self.sequence[i+(f-f_i):j] )
			self.frameshift_signals.append( self.sequence[j-3:j+3] )
			i = j
			f_i = f
			# we could factor in the last trivial frameshift...
		self.frameshifted_sequence += self.sequence[-1]
		self.fragments[-1] += self.sequence[-1]
			# or (preferably) allow the last fragment to run until the end
			#self.frameshifted_sequence += self.sequence[j:]
			#self.fragments[-1] += self.sequence[-1]
		
		self.path = path
			
		return self.frameshifted_sequence, self.fragments, self.frameshift_signals[:-1]
	
	def colour_frameshifted_sequence( self, frame=0, sep=" " ):
		"""
		Method to return in colour for frame frame
		"""
		# ensure that you have a valid frame; no need for errorcheck
		frame = frame % 3
		
		coloured_sequence = ""
	
		# front
		codon = self.frameshifted_sequence[0:frame]
		if codon in self.starts:
			codon = termcolor.colored( codon, 'yellow', 'on_yellow', attrs=['bold'] )
		elif codon in self.stops:
			codon = termcolor.colored( codon, 'white', 'on_red', attrs=['bold'] )
		elif codon in self.non_stops:
			codon = termcolor.colored( codon, 'blue', 'on_white', attrs=['bold'] )
		else:
			codon = termcolor.colored( codon, 'red', 'on_green', attrs=['bold'] )
	
			#		if frame % 3 != 0:
		coloured_sequence += codon + sep
	
		# body
		i = frame
		while i < len( self.frameshifted_sequence ):
			codon = self.frameshifted_sequence[i:i+3]
			if codon in self.starts:
				codon = termcolor.colored( codon, 'yellow', 'on_green', attrs=['bold'] )
			elif codon in self.stops:
				codon = termcolor.colored( codon, 'white', 'on_red', attrs=['bold'] )
			elif codon in self.non_stops:
				codon = termcolor.colored( codon, 'blue', 'on_white', attrs=['bold'] )
			else:
				codon = termcolor.colored( codon, 'red', 'on_green', attrs=['bold'] )
			coloured_sequence += codon + sep
			i += 3
		
		return coloured_sequence