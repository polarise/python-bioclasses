# -*- encoding: utf-8 -*-
from __future__ import division
import sys

class FrameshiftSite( object ):
	def __init__( self, initial_node, final_node, signal, length, position_score, \
		radians_vector, desig=None ):
		self.initial_frame = initial_node[0]
		self.final_frame = final_node[0]
		self.position = initial_node[1]
		if desig == None:
			desig = self.final_frame - self.initial_frame
			if desig == 1:
				self.designation = "+1"
			elif desig == 2:
				self.designation = "-1"
			elif desig == -1:
				self.designation = "-1"
			elif desig == -2:
				self.designation = "+1"
			else:
				self.designation = "-"
		else:
			self.designation = desig
		self.signal = signal
		self.distance_from_5prime = length - self.position
		self.position_score_f = position_score
		self.position_score = "%.4f" % position_score
		self.radians_vector_f = radians_vector
		self.radians_vector = "%.4f\t%.4f\t%.4f" % radians_vector
	
	def posscore2proportion( self, length ):
		if self.position < length/2:
			return self.position_score_f/2
		else:
			return 1 - self.position_score_f/2
		
	def __eq__( self, other ):
		if self.position == other.position and self.signal == other.signal:
			return True
		else:
			return False
	
	def __repr__( self ):
		return "\t".join( map( str, [ self.signal, self.distance_from_5prime, self.designation ] ))
		