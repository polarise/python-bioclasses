# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import scipy

class PSWM( object ):
	def __init__( self, data ): # data is a numpy array
		data = data + 1
		self.data = scipy.log10( data/data.sum( axis=0 )*4 )
		self.rows, self.cols = data.shape
	
	def __mul__( self, array ):
		if array.shape != self.data.shape:
			prod = -scipy.ones( self.data.shape )
		else:
			prod = array*self.data
		return prod		
	
	def __repr__( self ):
		return str( self.data )