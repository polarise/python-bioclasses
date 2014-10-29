# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import multiprocessing

class WorkerProcess( multiprocessing.Process ):
	def __init__( self, i, q1, q2, callback, var, lock ):
		"""
		i		- process_id
		q1	- input queue
		q2	- output queue
		"""
		multiprocessing.Process.__init__( self ) # inherit from ~
		self.i = i
		self.q1 = q1
		self.q2 = q2
		self.callback = callback # an externally-defined function that takes myobject as input
		self.var = var
		self.lock = lock
	
	def run( self ):
		"""
		overwrite the run method of the multiprocessing.Process class
		"""
		while True:
			myobject = self.q1.get()
			if myobject == "END": # the poison pill; assume that 'END' is not a valid myobject value
				self.q2.put( "END" ) # mark that this process is done
				self.q1.task_done() # unblock q1
				self.q2.task_done() # unblock q2
				break # die
			else:
				result = self.callback( myobject, self.var, self.lock )
				self.q2.put( result )
				self.q1.task_done() # unblock q1
				self.q2.task_done() # unblock q2