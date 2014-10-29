# -*- encoding: utf-8 -*-
from __future__ import division
import sys

class Read( object ):
	def __init__( self, name, trim_by=5, min_length=20 ):
		self.name = name
		self.bases = None
		self.name_ = None
		self.quals = None
		self.trim_by = trim_by
		self.min_length = min_length
	
	def remove_polyA( self, automatic=True ):
		if automatic:
			no_A = 0
			rev_seq = self.bases[::-1]
			i = 0
			is_A = True
			while is_A and i <= self.trim_by:
				l = rev_seq[i]
				if l != 'A':
					is_A = False
				i += 1
			if i > 3:
				self.bases = self.bases[:( len( self.bases ) - i + 1)]
				self.quals = self.quals[:( len( self.quals ) - i + 1)]			
		else:
			if self.bases[-self.trim_by:] == "A"*self.trim_by:
				self.bases = self.bases[:( len( self.bases ) - self.trim_by )]
				self.quals = self.quals[:( len( self.quals ) - self.trim_by )]
	
	def __repr__( self ):
		if len( self.bases ) >= self.min_length:
			return "\n".join([ self.name, self.bases, self.name_, self.quals ])