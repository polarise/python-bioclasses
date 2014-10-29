# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import multiprocessing
from WorkerProcess import *

class ParallelTask( object ):
	def __init__( self, input_container, callback, num_of_processors=1, var=None, lock=None ):
		self.input_container = input_container
		self.callback = callback
		self.var = var
		self.lock = lock
		if num_of_processors == 0:
			self.num_of_processors = multiprocessing.cpu_count()
		else:
			self.num_of_processors = num_of_processors
	
	def run( self ):
		# create the input and output queues
		input_queue = multiprocessing.JoinableQueue()
		output_queue = multiprocessing.JoinableQueue()
		
		for myobject in self.input_container:
			input_queue.put( myobject )
		
		# put the poison pill
		for i in xrange( self.num_of_processors ):
			input_queue.put( "END" )
		
		# create the WorkerProcess objects
		procs = [ WorkerProcess( i, input_queue, output_queue, self.callback, self.var, self.lock ) \
			for i in xrange( self.num_of_processors )]
		
		# start the processes
		for p in procs:
			p.daemon = True
			p.start()
		
		# block until the input queue joins
		input_queue.join()
		
		end = list()
		output_container = list()
		while len( end ) < self.num_of_processors:
			result = output_queue.get()
			if result != "END":
				output_container.append( result )
			else:
				end.append( result )
		
		# block until the output queue joins
		output_queue.join()
		
		return output_container