# -*- encoding: utf-8 -*-
import sys
from Node import *	
from Branch import *

class Paths( object ):
	"""
	A class that constructs paths from branches
	"""
	def __init__( self, frame_sequence ):
		"""
		initiate an empty Paths object
		"""
		self.frame_sequence = frame_sequence
		self.unique_frame_sequence = self.sanitise()
		self.paths = list()
		self.branches = list()
		self.all_path_sequences = list()
		
	def sanitise( self ):
		"""
		Method to remove duplicates and append terminal frames
		"""
		unique_frame_sequence = [ self.frame_sequence[0] ]
		for i in xrange( 1, len( self.frame_sequence ) ):
			if self.frame_sequence[i][0] == unique_frame_sequence[-1][0]:
				continue
			else:
				unique_frame_sequence.append( self.frame_sequence[i] )

		# print unique_frame_sequence
		unique_frame_sequence += zip( range( 3 ), [ -1 ]*3 )
		
		return unique_frame_sequence
	
	def create_branches( self ):
		"""
		Method to create the individual branches from which the tree is built
		"""
		branches = list()
		positions = list()
		unique_frame_sequence = self.unique_frame_sequence
		while len( unique_frame_sequence ) > 3:
			branch = list()
			position = list()
			i = 0
			while len( branch ) < 3:
				if unique_frame_sequence[i][0] not in branch:
					branch.append( unique_frame_sequence[i][0] )
					position.append( unique_frame_sequence[i][1] )
				i += 1
			branches.append( zip( branch, position ))
			unique_frame_sequence.pop( 0 )
		
		# create a list of Branch objects
		for b in branches:
			self.branches.append( Branch( Node( *b[0] ), Node( *b[1] ), Node( *b[2] )))
		
	def extend( self, branch ):
		"""
		Consider the branch and generate all possible paths in path
		"""
		new_paths = list() # a new list to put new paths
		if len( self.paths ) == 0: # if we're starting from scratch...
			for d in branch.descendants: # for each descendant in the branch...
				new_paths.append([ d, branch.root ] ) # create a list ending in the branch root
		else: 
			for d in branch.descendants: # for each descendant in the branch...
				for p in self.paths: # for each path we already had...
					if p[0] == branch.root: # if the leaf resembles the branch root...
						new_paths.append( [d] + p ) # extend this path; put in our new list
					else:
						if p not in new_paths: # otherwise, put only one copy in our new list
							new_paths.append( p )
		
		self.paths = new_paths # sub the paths with the updated list of paths
	
	def build_paths( self ):
		"""
		Method that builds all possible paths from leaves to the root using the 
		extend method (above)
		"""
		for B in self.branches:
			self.extend( B )
	
	def get_path_sequences( self, leaf ):
		"""
		Return a list of all paths having leaf as its terminus
		"""
		path_sequences = list()
		for p in self.paths:
			if p[0] == leaf:
				path_sequences.append( map( lambda x: ( x.name, x.value ), p[::-1] ))
		
		return path_sequences
	
	def get_all_path_sequences( self ):
		"""
		Return all paths present from the branches
		"""
		self.all_path_sequences = list()
		for p in self.paths:
			self.all_path_sequences.append( map( lambda x: ( x.name, x.value ), p[::-1] ))
		
		return self.all_path_sequences
	
	def get_frame_paths( self, frame ):
		frame %= 3 # validate frame
		
		frame_paths = list()
		for p in self.all_path_sequences:
			if p[0][0] == frame:
				frame_paths.append( p )
			elif p[1][0] == frame:
				frame_paths.append( p[1:] )
		
		return frame_paths