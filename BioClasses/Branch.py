# -*- encoding: utf-8 -*-
from __future__ import division
import sys

class Branch( object ):
	def __init__( self, root_node, node1, node2 ):
		self.root = root_node
		self.descendants = [ node1, node2 ]
	
	def __repr__( self ):
		branch_str = "Branch root:\n\t%s\n" % str( self.root )
		branch_str += "Descendants:\n\t%s, %s\n" % ( str( self.descendants[0] ), str( self.descendants[1] ))
		return branch_str	
