# -*- encoding: utf-8 -*-
from __future__ import division
import sys

class LeafCounter( object ):
	def __init__( self ):
		self.state_count = 0
		self.zero_count = 0
		self.one_count = 0
		self.two_count = 0
		self.branch_count = 0
	
	def add_node( self, node ):
		if self.state_count == 0:
			if node.name == 0:
				self.zero_count = 0
				self.one_count = 1
				self.two_count = 1
			elif node.name == 1:
				self.zero_count = 1
				self.one_count = 0
				self.two_count = 1
			elif node.name == 2:
				self.zero_count = 1
				self.one_count = 1
				self.two_count = 0
			self.branch_count = 0
		else:
			if node.name == 0:
				self.branch_count = self.zero_count
				self.zero_count = 0
				self.one_count += self.branch_count
				self.two_count += self.branch_count
			elif node.name == 1:
				self.branch_count = self.one_count
				self.zero_count += self.branch_count
				self.one_count = 0
				self.two_count += self.branch_count
			elif node.name == 2:
				self.branch_count = self.two_count
				self.zero_count += self.branch_count
				self.one_count += self.branch_count
				self.two_count = 0
		self.state_count += 1
	
	def leaf_count( self ):
		return self.zero_count + self.one_count + self.two_count
			