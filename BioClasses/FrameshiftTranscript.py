# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import scipy
from FrameshiftSite import *

class FrameshiftTranscript( object ):
	def __init__( self, name, length ):
		self.name = name
		self.length = length
		self.frameshift_sites = dict()
	
	def add_frameshift_site( self, position, signal, radians_vector=( 2*scipy.pi/3, 2*scipy.pi/3, 2*scipy.pi/3 ) ):
		def frameshift_position_score( x, L ):
			"""
			triangular function
			P( frameshift ) is maximum in the middle and decreases to the edges
			"""
			if x < L/2:
				return x/(L/2)
			else:
				return ( L - x )/(L/2)
		
		position_score = frameshift_position_score( position, self.length )
		
		self.frameshift_sites[position] = FrameshiftSite( ( 0, position ), \
			( 0, 0 ), signal, self.length, position_score, radians_vector )
	
	def __repr__( self ):
		output_str = "Transcript: %s of length %s\n" % ( self.name, self.length )
		i = 1
		for pos,FS in self.frameshift_sites.iteritems():
			output_str += "Frameshift #%s: %s at %s (pos-score = %s).\n" % ( i, \
				FS.signal, FS.position, FS.position_score )
			i += 1
		
		return output_str
	
	def filtered_print( self, p0=0, p1=1, theta0=scipy.pi ):
		output_str = "Transcript: %s of length %s\n" % ( self.name, self.length )
		i = 1
		for pos,FS in self.frameshift_sites.iteritems():
			if p0 <= FS.posscore2proportion( self.length ) <= p1 and FS.radians_vector_f[0] <= theta0:
				output_str += "Frameshift #%s: %s at %s (pos-score = %s).\n" % ( i, \
					FS.signal, FS.position, FS.position_score )
			i += 1
		
		return output_str
	
	def frameshifts( self, p0=0, p1=1, theta0=scipy.pi ):
		for fss_i,fss in self.frameshift_sites.iteritems():
			if p0 <= fss.posscore2proportion( self.length ) <= p1 and fss.radians_vector_f[0] <= theta0:
				yield self.name, fss
	
	def has_frameshift( self, p0=0, p1=1, theta0=scipy.pi ):
		"""
		beware!
		"""
		frameshift_count = 0
		for fss_i,fss in self.frameshift_sites.iteritems():
			if p0 <= fss.posscore2proportion( self.length ) <= p1 and fss.radians_vector_f[0] <= theta0:
				frameshift_count += 1
		if frameshift_count > 0:
			return True
		else:
			return False
	
	def has_exact_frameshift( self, other, p0=0, p1=1, theta0=scipy.pi, tol=3 ):
		"""
		beware!
		"""
		self_fsss = self.frameshift_sites.values()
		other_fsss = other.frameshift_sites.values()
		present = False
		for fss in self_fsss:
			for oss in other_fsss:
				if p0 <= fss.posscore2proportion( self.length ) <= p1 and fss.radians_vector_f[0] <= theta0 and -tol <=  fss.distance_from_5prime - oss.distance_from_5prime <= tol and fss.signal == oss.signal:
					present = True
		
		return present
		
	def rough_equality( self, other ):
		if len( self.frameshift_sites ) > 0 and len( other.frameshift_sites ) > 0:
			return True
		else:
			return False
		
	def is_equal( self, other, p0, p1, theta0 ):
		# each FSTObject has one or more FSSObjects
		# we look for equality on FSSObjects by comparing positions and signals
		equal = False
		number_equal = 0
		frameshift_sites_self = [ fss for fsss in self.frameshift_sites.values() if p0 <= fss.posscore2proportion( self.length ) <= p1 and fss.radians_vector <= theta0 ]
		frameshift_sites_other = other.frameshift_sites.values()
		for fsss in frameshift_sites_self:
			for fsso in frameshift_sites_other:
				if fsss == fsso:
					equal = True
					number_equal += 1
		
		return equal, number_equal
		
	
