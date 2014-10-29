from __future__ import division
import sys
import copy
from Node import *

class Tree( object ):
	def __init__( self ):
		self.root = None
		self.head = list()
		self.leaves = list()
		self.paths = list()
		
	def graft( self, branch, verbose ):
		"""
		Method to graft a branch to a growing tree
		"""
		if self.root is None: # a new tree
			if verbose: print >> sys.stderr, "Placing root %s to tree..." % branch
			self.root = branch.root
			self.root.left_leaf = copy.deepcopy( branch.descendants[0] )
			self.root.right_leaf = copy.deepcopy( branch.descendants[1] )
			self.head.append( self.root )
		else: # a growing tree
			if verbose: print >> sys.stderr, "Grafting %s to tree with %d leaves..." % ( branch, len( self.leaves ))
			for h in self.head:
				# left_leaf
				if h.left_leaf == branch.root:
					h.left_leaf.left_leaf, h.left_leaf.right_leaf = copy.deepcopy( branch.descendants )
					self.head.append( h.left_leaf )
				# right_leaf
				elif h.right_leaf == branch.root:
					h.right_leaf.left_leaf, h.right_leaf.right_leaf = copy.deepcopy( branch.descendants )
					self.head.append( h.right_leaf )
			
			# prune head: remove those without descendants
			if verbose: print >> sys.stderr, "Prunning head..."
			for h in self.head:
				if h.left_leaf is None and h.right_leaf is None:
					self.head.remove( h )
		
		# refresh leaves
		if verbose: print >> sys.stderr, "Refreshing leaves..."
		self.leaves = list()
		for l in self.head:
			if l.left_leaf.left_leaf is None and l.left_leaf.right_leaf is None:
				self.leaves.append( l.left_leaf )
			if l.right_leaf.right_leaf is None and l.right_leaf.left_leaf is None:
				self.leaves.append( l.right_leaf )
		if verbose: print >> sys.stderr
		
	def __repr__( self ):
		# tree_str = "Tree with the %d composite nodes and %d leaves.\n" % ( len( self.nodes ), len( self.leaves ))
		tree_str = "Tree with %d leaves.\n" % len( self.leaves )
		tree_str += "Root:\n\t%s\n" % str( self.root )
		# tree_str += "Nodes:\n\t%s\n" % ", ".join( map( str, self.nodes.values()))
		leave_str = ", ".join( map( str, self.leaves ))
		#tree_str += "Leaves:\n\t%s\n" % ", ".join( map( str, self.leaves ))
		tree_str += "Leaves:\n\t%s\n" % leave_str
		return tree_str
		
	def get_paths( self, simplify=False ):
		stack = list()
		count_leaves = len( self.leaves )
		current_node = self.root
		while count_leaves > 0:
			#print >> sys.stderr, stack
			if current_node.left_leaf is not None and current_node.right_leaf is not None:
				if not current_node.left_leaf.flagged and not current_node.right_leaf.flagged:
					stack.append( current_node )
					current_node = current_node.left_leaf
				elif current_node.left_leaf.flagged and not current_node.right_leaf.flagged:
					stack.append( current_node )
					current_node = current_node.right_leaf
				elif not current_node.left_leaf.flagged and current_node.right_leaf.flagged:
					print >> sys.stderr, "How did I get here?"
				elif current_node.left_leaf.flagged and current_node.right_leaf.flagged:
					current_node.flagged = True
					current_node = stack[-1]
					stack.pop()
			elif current_node.left_leaf is None and current_node.right_leaf is None:
				this_path = [ current_node ] + stack[::-1]
				self.paths.append( this_path[::-1] )
				current_node.flagged = True
				current_node = stack[-1]
				stack.pop()
				count_leaves -= 1
		
		if simplify: # return tuples
			return [ map( lambda x: ( x.name, x.value ), path ) for path in self.paths ]
		else: # return Node objects
			return self.paths
	
	def get_frame_paths( self, frame ):
		frame %= 3 # validate frame
		
		frame_paths = list()
		for p in self.paths:
			if p[0].name == frame:
				frame_paths.append( p )
			elif p[1].name == frame:
				frame_paths.append( p[1:] )
		
		return [ map( lambda x: ( x.name, x.value ), path ) for path in frame_paths ]