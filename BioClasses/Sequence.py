# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import itertools
import random
import math
import copy
import termcolor
import scipy.signal
from Node import *
from Branch import *
from Tree import *
from FrameshiftSequence import *
from GeneticCode import *
from TransitionMatrix import *
import numpy
import matplotlib.pyplot as plt
from LeafCounter import *
from Utils import *

class Sequence( object ):
	"""
	The base class
	"""
	def __init__( self, sequence=None, name=None, length=None, bases='ACGT', starts=[ 'ATG' ], stops=[ 'TAA', 'TAG' ] ):
		"""
		Initialise a Sequence object
		"""
		# basic stuff
		self.bases = bases
		self.starts = starts
		self.stops = stops
		
		itercodons = itertools.product( self.bases, repeat=3 )
		self.non_stops = [ "".join( codon ) for codon in itercodons ]
		for stop in self.stops:
			self.non_stops.remove( stop )
		
		self.name = name
		if sequence is None:
			self.sequence = sequence # frame 0
			self.sequence_ = sequence[1:] # frame 1
			self.sequence__ = sequence[2:] # frame 2
			self.length = length
		else:
			self.sequence = sequence.upper()
			self.sequence_ = sequence[1:].upper() # frame 1
			self.sequence__ = sequence[2:].upper() # frame 2
			self.length = len( sequence )
		
		self.CAI_score = None
		self.likelihood = None # frame 0
		self.likelihood_ = None # frame 1
		self.likelihood__ = None	# frame 2
		self.graded_likelihood = None # frame 0
		self.graded_likelihood_ = None # frame 1
		self.graded_likelihood__ = None # frame 2
		self.differential_graded_likelihood = None # frame 0
		self.differential_graded_likelihood_ = None # frame 1
		self.differential_graded_likelihood__ = None # frame 2
		self.gradient = None
		
		self.binary_codon_matrix = None
		
		self.is_frameshift = False
		self.as_codons = False
		
		self.start_pos = None
		self.start_positions = dict()
		self.stop_positions = dict()
		self.start_sequence = list()
		self.stop_sequence = list()
		self.filtered_stop_sequence = list()
		self.unique_stop_sequence = list()
		
		self.tree = None
		self.branches = list()		
		self.paths = list() # all tree paths
		self.no_paths = None
		self.frame_paths = dict() # all paths per frame
		self.sorted_frame_paths = dict() # sorted by length
		
		# dictionary of FrameshiftSequence objects
		# keys are path tuples
		self.frameshift_sequences = dict()
		self.most_likely_frameshift = None
		self.least_likely_frameshift = None
		self.frameshift_signals = list()
		
		# translation
		self.genetic_code = None
		self.transition_matrix = None
		
	#*****************************************************************************
	
	def count_leaves( self ):
		"""
		Method to count the number of leaves in the tree
		"""
		self.get_stop_sequence()
		self.sanitise_stop_sequence()
		nodes = [ Node( *d ) for d in self.unique_stop_sequence ]
		# initiate a LeafCounter object
		L = LeafCounter()
		for n in nodes[:-3]:
			L.add_node( n )
			
		return L.leaf_count()
	
	#*****************************************************************************
	
	def truncate( self, start_from="ATG", effect_truncation=True, verbose=True ):
		"""
		
		"""
		self.start_pos = None
		self.start_pos = self.sequence.find( start_from )
		
		if effect_truncation: 
			if self.start_pos < 0:
				if verbose:
					print >> sys.stderr, "Warning: unable to find %s signal... using whole sequence." % start_from
				return -1
			else:
				if verbose:
					print >> sys.stderr, "Found start (%s) from position %d... truncating sequence" % ( start_from, self.start_pos )
				self.sequence = self.sequence[self.start_pos:] # frame 0
				self.sequence_ = self.sequence[1:]	# frame 1
				self.sequence__ = self.sequence[2:]	# frame 2
				self.length = len( self.sequence )
				self.start_pos = 0
				return 0
		else: # What if there is not 'start'? - fixed PK 20140806
			if verbose:
				if self.start_pos < 0:
					print >> sys.stderr, "Warning: unable to find %s signal... using whole sequence." % start_from
					return None
				else:
					print >> sys.stderr, "Not effecting truncation but return position of first start (ATG)..."
					return self.start_pos
	
	#*****************************************************************************
	
	def colour_frame( self, frame, sep=" " ):
		"""
		Method to return in colour for frame frame
		"""
		# ensure that you have a valid frame; no need for errorcheck
		frame = frame % 3
		
		coloured_sequence = ""
	
		# front
		codon = self.sequence[0:frame]
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
		while i < len( self.sequence ):
			codon = self.sequence[i:i+3]
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
	
	#*****************************************************************************
	
	def __repr__( self ):
		if self.as_codons:
			sequence = list()
			for i in xrange( 0, len( self.sequence ), 3 ):
				sequence.append( self.sequence[i:i+3] )
			return " ".join( sequence )
		else:
			return self.sequence
	
	#*****************************************************************************
	
	def repr_as_row( self, sep="\t" ):
		return sep.join([ "...".join([ self.sequence[:20], 
						self.sequence[-20:] ]), str( self.length ), \
							"0", str( self.CAI_score ), str( self.CAI_score ), \
								"None", "None" ])
	
	#*****************************************************************************
	
	def repr_frameshift_sites( self, sep="\t", include_nulls=True ):
		if self.most_likely_frameshift is None:
			self.get_most_likely_frameshift()
		if self.most_likely_frameshift is not None:
			self.most_likely_frameshift.find_frameshift_sites()
			ML = self.most_likely_frameshift
			ML.find_frameshift_sites()
			if len( ML.frameshift_sites ) == 0 and include_nulls:
				print sep.join( map( str, [ self.name, ML.length, ML.likelihood, \
					None, None, None, None, None, None, None ] ))
			else:
				for fss in ML.frameshift_sites:
					FS = ML.frameshift_sites[fss]
					print sep.join( map( str, [ self.name, ML.length, ML.likelihood, \
						FS.signal, FS.initial_frame, FS.final_frame, FS.designation, \
							FS.position, FS.position_score, FS.radians_vector ] ))			
	
	#*****************************************************************************
		
	def info( self, comment="" ):
		raise NotImplementedError
	
	#*****************************************************************************
	
	def binary_frame( self, frame, sep=" " ):
		"""
		Method that returns a binary codon sequence: 0 - stop; 1 - non-stop;
		red: stop; blue: start
		"""
		binary_sequence = ""
			
		# body
		i = frame
		while i < self.length + frame:
			if i < 0:
				bit = termcolor.colored( "0", 'red', 'on_green', attrs=['bold'] )			
			else:
				codon = self.sequence[i:i+3]
				if codon in self.starts:
					bit = termcolor.colored( "1", 'white', 'on_blue', attrs=['bold'] )
				elif codon in self.stops:
					bit = termcolor.colored( "0", 'white', 'on_red', attrs=['bold'] )
				elif codon in self.non_stops:
					bit = termcolor.colored( "1", 'blue', 'on_white', attrs=['bold'] )
				else:
					bit = termcolor.colored( "0", 'red', 'on_green', attrs=['bold'] )
			binary_sequence += bit + sep
			i += 3
		
		return binary_sequence
	
	#*****************************************************************************
		
	def codon_frame_vector( self, frame, SS=10, ES=0, NS=1, IS=-1, weighted_start=False ):
		"""
		SS - start codon score
		ES - stop (end) codon score
		NS - non-stop codon score
		IS - incomplete codon score
		"""
		frame_vector = list()
		
		# body
		i = frame
		while i < self.length + frame:
			if i < 0:
				ind = IS
			else:
				codon = self.sequence[i:i+3]
				if codon in self.starts:
					if weighted_start:
						ind = max( NS, SS - len( frame_vector )) # arbitrary value
					else:
						ind = SS
				elif codon in self.stops:
					ind = ES
				elif codon in self.non_stops:
					ind = NS
				else:
					ind = IS
			frame_vector.append( ind )
			i += 3
		
		return frame_vector
	
	#*****************************************************************************
	
	def set_binary_codon_matrix( self, min_frame=-2, max_frame=2, SS=10, ES=0, NS=1, IS=-1, weighted_start=False ):
		self.binary_codon_matrix = dict()
		for i in xrange( min_frame - 2, max_frame + 3, 1 ):
			row = self.codon_frame_vector( i, SS, ES, NS, IS, weighted_start )
			self.binary_codon_matrix[i] = row
			
	#*****************************************************************************
	
	def show_binary_codon_matrix( self ):
		keys = self.binary_codon_matrix.keys()
		keys.sort( reverse=True )
		
		max_frame = max( map( abs, keys ))
		rjust_int = int( math.floor( math.log( max_frame, 10 ) + 1 ))
		
		output_str = "BinaryCodonMatrix of length %d:\n" % len( self.binary_codon_matrix[0] )
		for k in keys:
			if k > 0:
				output_str += "+%s: %s\n" % ( str( k ).rjust( rjust_int ), " ".join( map( lambda x: x.rjust( 2 ), map( str, self.binary_codon_matrix[k] ))))
			elif k == 0:
				output_str += " %s: %s\n" % ( str( k ).rjust( rjust_int ), " ".join( map( lambda x: x.rjust( 2 ), map( str, self.binary_codon_matrix[k] ))))
			elif k < 0:
				output_str += "%s: %s\n" % ( str( k ).rjust( rjust_int ), " ".join( map( lambda x: x.rjust( 2 ), map( str, self.binary_codon_matrix[k] ))))
		
		return output_str.strip( "\n" )
		
	#*****************************************************************************
	
	def get_start_positions( self ):
		"""
		Method to return the position of start codons in all frames
		"""
		start_pos = dict()
		
		for frame in xrange( 3 ):
			i = frame
			while i < self.length:
				codon = self.sequence[i:i+3]
				if codon in self.starts:
					start_pos[i] = frame
				i += 3
		
		return start_pos
	
	#*****************************************************************************
	
	def get_stop_positions( self ):
		"""
		Method to return the position of stop codons in all frames
		"""		
		stop_pos = dict()
		
		for frame in xrange( 3 ):
			i = frame
			while i < self.length:
				codon = self.sequence[i:i+3]
				if codon in self.stops:
					stop_pos[i] = frame
				i += 3

		return stop_pos
		
	#*****************************************************************************
	
	def get_start_sequence( self ):
		self.start_positions = list()
		self.start_positions = self.get_start_positions()
		positions = self.start_positions.keys()
		positions.sort()
		
		for p in positions:
			self.start_sequence.append( ( self.start_positions[p], p ))
		
		return self.start_sequence
	
	#*****************************************************************************
	
	def get_stop_sequence( self ):
		self.stop_positions = list()
		self.stop_positions = self.get_stop_positions()
		positions = self.stop_positions.keys()
		positions.sort()
		
		for p in positions:
			self.stop_sequence.append( ( self.stop_positions[p], p ))
		
		return self.stop_sequence
			
	#*****************************************************************************
	
	def sanitise_stop_sequence( self ):
		"""
		Method to remove duplicates and append terminal frames
		"""
		stop_sequence = copy.copy( self.stop_sequence )
		
		if len( self.stop_sequence ) > 0:
			self.unique_stop_sequence.append( stop_sequence[0] )
			for i in xrange( 1, len( stop_sequence ) ):
				if stop_sequence[i][0] == self.unique_stop_sequence[-1][0]:
					continue
				else:
					self.unique_stop_sequence.append( stop_sequence[i] )
		else:
			self.unique_stop_sequence = list()

		# print unique_stop_sequence
		# if len( self.unique_stop_sequence ) == 0:
		self.unique_stop_sequence += zip( range( 3 ), [ -1 ]*3 )
		# else:
		# 	if self.unique_stop_sequence[-1][0] == 0:
		# 		self.unique_stop_sequence += [ (1,-1), (2,-1) ]
		# 	elif self.unique_stop_sequence[-1][0] == 1:
		# 		self.unique_stop_sequence += [ (0,-1), (2,-1) ]
		# 	elif self.unique_stop_sequence[-1][0] == 2:
		# 		self.unique_stop_sequence += [ (0,-1), (1,-1) ]
		
		return self.unique_stop_sequence
		
	#*****************************************************************************
	
	def create_branches( self ):
		"""
		Method to create the individual branches from which the tree is built
		"""
		branches = list()
		positions = list()
		unique_stop_sequence = copy.copy( self.unique_stop_sequence )
		while len( unique_stop_sequence ) > 3:
			branch = list()
			position = list()
			i = 0
			while len( branch ) < 3:
				if unique_stop_sequence[i][0] not in branch:
					branch.append( unique_stop_sequence[i][0] )
					position.append( unique_stop_sequence[i][1] )
				i += 1
			branches.append( zip( branch, position ))
			unique_stop_sequence.pop( 0 )
		
		# create a list of Branch objects
		for b in branches:
			self.branches.append( Branch( Node( *b[0] ), Node( *b[1] ), Node( *b[2] )))

	#*****************************************************************************
	
	def build_tree( self, verbose=False ):
		"""
		Method to build a tree associated with the sequence
		"""
		# initialise
		self.start_sequence = list()
		self.stop_sequence = list()
		self.unique_stop_sequence = list()
		
		# initialise the tree and branches
		self.tree = Tree()
		self.branches = list()
		
		# initialise other attributes
		self.paths = list() # all tree paths
		self.no_paths = None
		self.frame_paths = dict() # all paths per frame
		self.sorted_frame_paths = dict() # sorted by length
		
		self.frameshift_sequences = dict()
		
		# get the raw start sequence
		self.get_start_sequence()
		if verbose: print >> sys.stderr, "Generated start sequence of length %d..." % len( self.start_sequence )
		
		# get the raw stop sequence
		self.get_stop_sequence()
		if verbose: print >> sys.stderr, "Generated stop sequence of length %d..." % len( self.stop_sequence )
		
		# sanitise the raw stop sequence
		self.sanitise_stop_sequence()
		if verbose: print >> sys.stderr, "Generated unique stop sequence of length %d..." % len( self.unique_stop_sequence )
		
		# create a list of branches
		self.create_branches()
		if verbose: print >> sys.stderr, "Creating branches..."
		
		# graft the branches to the tree
		if verbose: print >> sys.stderr, "Grafting branches to tree..."
		for B in self.branches:
			self.tree.graft( B, verbose )
		
		# get the paths
		self.paths = self.tree.get_paths( simplify=True )
		self.no_paths = len( self.paths )
		if verbose: print >> sys.stderr, "Found %d paths in tree." % self.no_paths
		
		# get paths per frame
		for frame in xrange( 3 ):
			self.frame_paths[frame] = self.tree.get_frame_paths( frame )
			
		# sort frame paths by length
		sorted_frame_paths = dict()
		for frame in xrange( 3 ):
			for path in self.frame_paths[frame]:
				path_len = len( path )
				if path_len not in sorted_frame_paths:
					sorted_frame_paths[path_len] = [ path ]
				else:
					sorted_frame_paths[path_len] += [ path ]
		
		all_path_len = sorted_frame_paths.keys()
		all_path_len.sort()
		
		for a in all_path_len:
			if a not in self.sorted_frame_paths:
				self.sorted_frame_paths[a] = sorted_frame_paths[a]
			else:
				self.sorted_frame_paths[a] += sorted_frame_paths[a]
		
		# get all frameshift sequences
		for frame in self.frame_paths:
			for path in self.frame_paths[frame]:
				self.frameshift_sequences[tuple(path)] = FrameshiftSequence( self.sequence, path )
	
	#*****************************************************************************
	
	def set_genetic_code( self, genetic_code ):
		self.genetic_code = genetic_code
	
	#*****************************************************************************
	
	def set_transition_matrix( self, transition_matrix ):
		self.transition_matrix = transition_matrix
	
	#*****************************************************************************
		
	def translate( self, genetic_code=None ):
		if genetic_code is None and self.genetic_code is None:
			raise ValueError( """Translation cannot proceed without a valid \
				genetic code.""" )
	
	#*****************************************************************************
	
	def estimate_CAI( self ):
		if self.genetic_code is None:
			raise ValueError( "Unable to compute CAI without genetic code. First run 's.set_genetic_code( G )'." )
		else:
			score = 0
			i = 0
			length = 0
			while i <= len( self.sequence ) - 3:
				codon = self.sequence[i:i+3]
				if len( codon ) < 3:
					print >> sys.stderr, "Warning: sequence end is not a codon: %s" % codon
					break
				score += math.log( self.genetic_code.get_wij( codon ))
				length += 1
				i += 3
			self.CAI_score = math.exp( score/length )
		
			return self.CAI_score
	
	#*****************************************************************************
	
	def estimate_frameshift_CAI( self ):
		if self.genetic_code is None:
			raise ValueError( "Unable to compute CAI without genetic code. First run 's.set_genetic_code( G )'." )
		else:
			for fs in self.frameshift_sequences:
				F = self.frameshift_sequences[fs]
				score = 0
				i = 0
				length = 0
				while i <= len( F.frameshifted_sequence ) - 3:
					codon = F.frameshifted_sequence[i:i+3]
					if len( codon ) < 3:
						print >> sys.stderr, "Warning: sequence end is not a codon: %s" % codon
						break
					score += math.log( self.genetic_code.get_wij( codon ))
					length += 1
					i += 3
				F.CAI = math.exp( score/length )
	
	#*****************************************************************************
	
	def estimate_likelihood( self, loglik=True ):
		if self.transition_matrix is None:
			raise ValueError( "Missing transition matrix. \
				First run 's.set_transition_matrix( TM )'." )
		else:
			# frame 0
			self.likelihood = self.transition_matrix.likelihood( self.sequence, \
				loglik=loglik )
			self.graded_likelihood = self.transition_matrix.graded_likelihood( \
				self.sequence, loglik=loglik, initial=True )
			self.differential_graded_likelihood = \
				self.transition_matrix.differential_graded_likelihood( self.sequence, \
					loglik=loglik, initial=True )
			# frame 1
			self.likelihood_ = self.transition_matrix.likelihood( self.sequence_, \
				loglik=loglik, initial=True )
			self.graded_likelihood_ = self.transition_matrix.graded_likelihood( \
				self.sequence_, loglik=loglik, initial=True )
			self.differential_graded_likelihood_ = \
				self.transition_matrix.differential_graded_likelihood( self.sequence_, \
					loglik=loglik, initial=True )
			# frame 2
			self.likelihood__ = self.transition_matrix.likelihood( self.sequence__, \
				loglik=loglik, initial=True )
			self.graded_likelihood__ = self.transition_matrix.graded_likelihood( \
				self.sequence__, loglik=loglik, initial=True )
			self.differential_graded_likelihood__ = \
				self.transition_matrix.differential_graded_likelihood( self.sequence__,\
					loglik=loglik, initial=True )
			
			# only return results for frame 0
			return self.likelihood, self.graded_likelihood, \
				self.differential_graded_likelihood
			
	#*****************************************************************************
	
	def estimate_frameshift_likelihood( self, loglik=True ):
		if self.transition_matrix is None:
			raise ValueError( "Missing transition matrix. \
				First run 's.set_transition_matrix( TM )'." )
		else:
			for fs in self.frameshift_sequences:
				F = self.frameshift_sequences[fs]
				F.likelihood = self.transition_matrix.likelihood( \
					F.frameshifted_sequence, loglik=loglik, initial=True )
				F.graded_likelihood = self.transition_matrix.graded_likelihood( \
					F.frameshifted_sequence, loglik=loglik, initial=True )
				F.differential_graded_likelihood = \
					self.transition_matrix.differential_graded_likelihood( \
						F.frameshifted_sequence, loglik=loglik, initial=True )
				#F.estimate_partial_gradients()
	
	#*****************************************************************************
	
	def estimate_gradient( self ):
		self.gradient = self.differential_graded_likelihood[-1]/self.length
		return self.gradient
	
	#*****************************************************************************
	
	def get_most_likely_frameshift( self ):
		max_likelihood = None
		min_likelihood = None
		self.most_likely_frameshift = None
		self.least_likely_frameshift = None
		for fs in self.frameshift_sequences:
			F = self.frameshift_sequences[fs]
			if max_likelihood is None:
				self.most_likely_frameshift = F
				self.least_likely_frameshift = F
				max_likelihood = F.likelihood
				min_likelihood = F.likelihood
			else:
				if F.likelihood > max_likelihood:
					max_likelihood = F.likelihood
					self.most_likely_frameshift = F
				if F.likelihood < min_likelihood:
					min_likelihood = F.likelihood
					self.least_likely_frameshift = F
				
		return self.most_likely_frameshift
	
	#*****************************************************************************
		
	def get_frameshift_signals( self ):
		"""
		"""
		self.frameshift_signals = list()
		# first get all frame of sequence
		sequence_in_frames = dict()
		for i in xrange( 3 ):
			sequence_in_frames[i] = self.sequence[i:]
	
		i = 0
		f_i = 0
		for f,j in self.unique_stop_sequence[:-3]:
			self.frameshift_signals.append( self.sequence[j-3:j+3] )
			i = j
			f_i = f
			
		return self.frameshift_signals[:-1]
			
	#*****************************************************************************
	
	def get_indexes( self ):
		"""
		Method to determine likely frameshift sequences
		for each frameshift sequence
			for each fragment
				compute the differential_graded_likelihood
				compute the slope
				create a vector of slopes for all three frames
				compare the vector to ( 1, -1, -1 ), ( -1, 1, -1 ), and ( -1, -1, 1 )
					to determine the right frame
					we're looking for frameshifts in which adjacent frames match in
					the zeroth frame
		"""
		# local functions
		def my_dot0( v1 ):
			"""
			frame 0
			"""
			v2 = ( 1, -1, -1 )
			#v2 = ( 1, -0.5, -0.5 )
			try:
				result = numpy.dot( v1, v2 )/numpy.linalg.norm( v1 )/numpy.linalg.norm( v2 )
			except TypeError:
				result = None
			return result

		def my_dot1( v1 ):
			"""
			frame 1
			"""
			v2 = ( -1, 1, -1 )
			#v2 = ( -0.5, 1, -0.5 )
			try:
				result = numpy.dot( v1, v2 )/numpy.linalg.norm( v1 )/numpy.linalg.norm( v2 )
			except TypeError:
				result = None
			return result

		def my_dot2( v1 ):
			"""
			frame 2
			"""
			v2 = ( -1, -1, 1 )
			#v2 = ( -0.5, -0.5, 1 )
			try:
				result = numpy.dot( v1, v2 )/numpy.linalg.norm( v1 )/numpy.linalg.norm( v2 )
			except TypeError:
				result = None
			return result

		def my_arccos( v ):
			if v is not None:
				return numpy.arccos( v )
			else:
				return None
		
		def index_of_mins( radians ):
			indexes = list()
			for r in radians:
				# allowing ( None, None, None ) to return 0
				if r[0] == None and r[1] == None and r[2] == None:
					index = None
				else:
					index = r.index( min( r ))
				indexes.append( index )
			
			return indexes
		
		radians = list()
		
		for fs in self.frameshift_sequences:
			F = self.frameshift_sequences[fs]
			# differential gradient likelihood in each frame
			dgl_list0 = map( self.transition_matrix.differential_graded_likelihood, \
				F.fragments )
			dgl_list1 = map( lambda x: \
				self.transition_matrix.differential_graded_likelihood( x[1:] ), \
					F.fragments )
			dgl_list2 = map( lambda x: \
				self.transition_matrix.differential_graded_likelihood( x[2:] ), \
					F.fragments )
			
			# slopes
			vector_list = zip( map( self.transition_matrix.likelihood_slope, \
				dgl_list0 ), map( self.transition_matrix.likelihood_slope, \
					dgl_list1 ), map( self.transition_matrix.likelihood_slope, \
						dgl_list2 ))
			
			# arccos values
			#print vector_list
			arccos0 = map( my_dot0, vector_list )
			arccos1 = map( my_dot1, vector_list )
			arccos2 = map( my_dot2, vector_list )
			
			# angles in radians
			rad0 = map( my_arccos, arccos0 )
			rad1 = map( my_arccos, arccos1 )
			rad2 = map( my_arccos, arccos2 )
			
			# radians
			radians = zip( rad0, rad1, rad2 )
			F.radians = radians
			radian_sums = list()
			for r in radians:
				the_sum = 0
				try:
					the_sum = sum( r )
				except TypeError:
					the_sum = None
				radian_sums.append( the_sum )
			F.radian_sums = radian_sums
			F.indexes = index_of_mins( radians )
		
	#*****************************************************************************
	def plot_differential_graded_likelihood( self, outfile=None, show_first_start=True, show_starts=False, show_signals=True, show_path_str=True, show_name=True, show_ML=False ):
		"""
		Method to plot a sequence and its likelihood tributaries
		"""
		fig = plt.figure()
		
		left, width = 0.1, 0.8
		
		# left, bottom, width, height
		rect1 = [ left, 0.5, width, 0.45 ]
		rect2 = [ left, 0.36, width, 0.14 ]
		rect3 = [ left, 0.22, width, 0.14 ]
		rect4 = [ left, 0.08, width, 0.14 ]
		
		
		ax1 = fig.add_axes( rect1 )
		ax2 = fig.add_axes( rect2, sharex=ax1 )
		ax3 = fig.add_axes( rect3, sharex=ax1 )
		ax4 = fig.add_axes( rect4, sharex=ax1 )
		
		# frame 0
		x = numpy.linspace( 1, len( self.differential_graded_likelihood ), \
			len( self.differential_graded_likelihood ))
		ax1.plot( x*3, self.differential_graded_likelihood, ":", color='r', \
			linewidth=1.5 )
		ax1.annotate( "No shift (fr 0)", xy=( self.length + 4, \
			self.differential_graded_likelihood[-1] ), size='x-small', \
				horizontalalignment='left' )
		
		# frame 1
		x1 = numpy.linspace( 1, len( self.differential_graded_likelihood_ ), \
			len( self.differential_graded_likelihood_ ))
		ax1.plot( x1*3 + 1, self.differential_graded_likelihood_, ":", color='g', \
			linewidth=1.5 )
		ax1.annotate( "No shift (fr 1)", xy=( self.length + 4, \
			self.differential_graded_likelihood_[-1] ), size='x-small', \
				horizontalalignment='left' )
		
		# frame 2
		x2 = numpy.linspace( 1, len( self.differential_graded_likelihood__ ), \
			len( self.differential_graded_likelihood__ ))
		ax1.plot( x2*3 + 2, self.differential_graded_likelihood__, ":", color='b', \
			linewidth=1.5 )
		ax1.annotate( "No shift (fr 2)", xy=( self.length + 4, \
			self.differential_graded_likelihood__[-1] ), size='x-small', \
				horizontalalignment='left' )
		
		"""
		print >> sys.stderr, len( self.differential_graded_likelihood ), len( self.differential_graded_likelihood_ ), len( self.differential_graded_likelihood__ )
		the_length = min([ len( self.differential_graded_likelihood ), len( self.differential_graded_likelihood_ ), len( self.differential_graded_likelihood__ )])
		val = list()
		for i in xrange( the_length ):
			# val.append( math.sqrt( self.differential_graded_likelihood[i]**2 + self.differential_graded_likelihood_[i]**2 + self.differential_graded_likelihood__[i]**2 ))
			val.append( self.differential_graded_likelihood[i] + self.differential_graded_likelihood_[i] + self.differential_graded_likelihood__[i] )
		xk = numpy.linspace( 1, the_length, the_length )
		plt.plot( xk*3, val, linewidth=2 )
		"""
		
		ax1.set_xlim(( 0, self.length + 40 ))
		
		# x- and y-labels
		ax1.set_ylabel( r"$\lambda(\phi(t))$" )
		
		# the frameshift sequences
		up = True
		for path in self.paths:
			F = self.frameshift_sequences
			# only plot the name of the most_likely_frameshift
			# note: path[1:]
			# why? because starting with the first frame is pointless
			x = numpy.linspace( 1, \
				len( F[tuple( path[1:] )].differential_graded_likelihood ), \
					len( F[tuple( path )].graded_likelihood ))
			ax1.plot( x*3, F[tuple( path )].differential_graded_likelihood )
			# write the shift sequence at the end
			# add the indicator sequence
			if show_path_str:
				if up:
					ax1.annotate( F[tuple( path[1:] )].path_str, xy=( self.length + 4, \
						F[tuple( path )].differential_graded_likelihood[-1] + 0.25 ), \
							size='xx-small', horizontalalignment='left' )
					up = False
				else:
					ax1.annotate( F[tuple( path[1:] )].path_str, xy=( self.length + 4, \
						F[tuple( path )].differential_graded_likelihood[-1] - 0.25 ), \
							size='xx-small', horizontalalignment='left' )
					up = True
			
		# the frameshift sites
		# the frameshift signal
		ymin,ymax = plt.ylim()
		xmin,xmax = plt.xlim()
		if show_signals:
			for i in xrange( len( self.frameshift_signals )):
				# vertical dashed frameshift signal markers
				ax1.axvline( self.unique_stop_sequence[i][1], color=( .5, .5, .5 ), \
					linestyle="dashed" )
				# the frame
				if i == 0:
					ax1.annotate( self.unique_stop_sequence[i][0], \
						xy=( self.unique_stop_sequence[i][1]/2, ymin + 3 ), size='x-small', \
							color=( 0.5, 0.5, 0.5 ), horizontalalignment='center' )
				else:
					ax1.annotate( self.unique_stop_sequence[i][0], \
						xy=( self.unique_stop_sequence[i-1][1] + \
							( self.unique_stop_sequence[i][1] - \
								self.unique_stop_sequence[i-1][1] )/2, ymin + 3 ), \
									size='x-small', color=( 0.5, 0.5, 0.5 ), \
										horizontalalignment='center' )
				# frameshift signal
				ax1.annotate( self.frameshift_signals[i], \
					xy=( self.unique_stop_sequence[i][1], ymax ), rotation=90, \
						size='x-small', horizontalalignment='right', \
							verticalalignment='right' )
		
		# terminal region
		ax1.axvline( self.length - 1, color='r', linestyle="dashed" )
	
		# mark the position of the first start (ATG)
		if show_first_start:
			#if self.start_pos >= 0:
			ax1.axvline( self.start_pos, color='r', linestyle='dashed' )
			ax1.annotate( "ATG", xy=( self.start_pos, ymax ), rotation=90, \
				size='x-small', color='r', horizontalalignment='right', \
					verticalalignment='right' )
		
		# mark the positions of all starts
		if show_starts:
			if len( self.start_sequence ) > 0:
				for fr,pos in self.start_sequence:
					if fr == 0:
						ax1.axvline( pos, color='r', linestyle='dashed' )
						ax1.annotate( "ATG[0]", xy=( pos, ymax ), rotation=90,\
							size='x-small', color='r', horizontalalignment='right',\
								verticalalignment='right' )
					elif fr == 1:
						ax1.axvline( pos, color='g', linestyle='dashed' )
						ax1.annotate( "ATG[1]", xy=( pos, ymax ), rotation=90,\
							size='x-small', color='g', horizontalalignment='right',\
								verticalalignment='right' )
					elif fr == 2:
						ax1.axvline( pos, color='b', linestyle='dashed' )
						ax1.annotate( "ATG[2]", xy=( pos, ymax ), rotation=90,\
							size='x-small', color='b', horizontalalignment='right',\
								verticalalignment='right' )
		
		# the sequence name (number of paths)
		if show_name:
			ax1.annotate( "%s (%s paths)" % ( self.name, len( self.paths )), \
				xy=( xmin + 0.05*( xmax - xmin ), ymax - 0.05*( ymax - ymin ) ), \
					size='large', horizontalalignment='left', verticalalignment='top', \
						bbox=dict( boxstyle="square", ec=( 1, .5, .5 ), fc=( 1, 1, 1 )))
		
		if show_ML:
			ML = self.most_likely_frameshift
			ax1.annotate( "MLFS: %s %s" % ( ML.path, ML.indexes ), \
				xy=( xmin + 0.05*( xmax - xmin ), ymin + 0.15*( ymax - ymin ) ), \
					size='small', horizontalalignment='left', verticalalignment='top', \
						bbox=dict( boxstyle="square", ec=( 1, .5, .5 ), fc=( 1, 1, 1 )))
							 
		# add a grid
		ax1.grid()
		
		"""
		# the filtered gradient plot
		xml = numpy.linspace( 1, len( self.most_likely_frameshift.differential_graded_likelihood ), \
			len( self.most_likely_frameshift.differential_graded_likelihood ))
		#print >> sys.stderr, len( xml ), self.most_likely_frameshift.length, len( self.most_likely_frameshift.differential_graded_likelihood )
		ax2.plot( xml*3, savgol( self.most_likely_frameshift.differential_graded_likelihood, \
			window_size=5, deriv=1 ))
		"""
		
		# frame 0
		dgl_prime = scipy.signal.savgol_filter( self.differential_graded_likelihood, \
			window_length=7, polyorder=3, deriv=1 )
		x = numpy.linspace( 1, len( dgl_prime ), len( dgl_prime ))
		ax2.plot( x*3, dgl_prime, color='r', linewidth=1.5 )
		if show_signals:
			for i in xrange( len( self.most_likely_frameshift.path )):
				# vertical dashed frameshift signal markers
				ax2.axvline( self.most_likely_frameshift.path[i][1], color='r', \
					linestyle="dashed" )
		ax2.set_ylabel( r"$m_0$" )
		ax2.grid()
		
		# frame 1
		dgl_prime_ = scipy.signal.savgol_filter( self.differential_graded_likelihood_, \
			window_length=7, polyorder=3, deriv=1 )
		x1 = numpy.linspace( 1, len( dgl_prime_ ), len( dgl_prime_ ))
		ax3.plot( x1*3 + 1, dgl_prime_, color='g', linewidth=1.5 )
		if show_signals:
			for i in xrange( len( self.most_likely_frameshift.path )):
				# vertical dashed frameshift signal markers
				ax3.axvline( self.most_likely_frameshift.path[i][1], color='r', \
					linestyle="dashed" )
		ax3.set_ylabel( r"$m_1$" )
		ax3.grid()
		
		# frame 2
		dgl_prime__ = scipy.signal.savgol_filter( self.differential_graded_likelihood__, \
			window_length=7, polyorder=3, deriv=1 )
		x2 = numpy.linspace( 1, len( dgl_prime__ ), len( dgl_prime__ ))
		ax4.plot( x2*3 + 2, dgl_prime__, color='b', linewidth=1.5 )
		if show_signals:
			for i in xrange( len( self.most_likely_frameshift.path )):
				# vertical dashed frameshift signal markers
				ax4.axvline( self.most_likely_frameshift.path[i][1], color='r', \
					linestyle="dashed" )
		ax4.set_ylabel( r"$m_2$" )
		ax4.set_xlabel( "Sequence position, $t$ (bp)" )
		ax4.grid()
		
		"""
		if show_signals:
			for i in xrange( len( self.most_likely_frameshift.path )):
				# vertical dashed frameshift signal markers
				ax2.axvline( self.most_likely_frameshift.path[i][1], color='r', \
					linestyle="dashed" )
		"""
		
		# space between plots
		#fig.subplots_adjust( hspace=0.1 )
		
		if outfile is not None:
			outfile.savefig()
			plt.close()
		else:
			plt.show()
		