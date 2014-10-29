# -*- encoding: utf-8 -*-
from __future__ import division
import sys
import math
import itertools
import cPickle
from Bio import SeqIO

class TransitionMatrix( object ):
	def __init__( self ):
		"""
		Constructor
		"""
		self.transition_counts = dict()
		self.initial_counts = dict()
		self.transition_probabilities = dict()
		self.initial_probabilities = dict()
	
	#*****************************************************************************
	
	def __repr__( self ):
		data = ""
		keys = self.transition_probabilities.keys()
		keys.sort()
		data += " " + " ".join( keys ) + "\n"
		for codon in keys:
			row = [ codon ]
			for next_codon in keys:
				row.append( "%.3f" % self.transition_probabilities[codon][next_codon] )
			data += " ".join( row ) + "\n"
		
		return data
	
	#*****************************************************************************
	
	def row_sums( self ):
		keys = self.transition_probabilities.keys()
		keys.sort()
		
		for k in keys:
			print k, sum( self.transition_probabilities[k].values() )
	
	#*****************************************************************************
	
	def col_sums( self ):
		keys = self.transition_probabilities.keys()
		keys.sort()
		
		for k in keys:
			colsum = 0
			for l in keys:
				colsum += self.transition_probabilities[l][k]
			print k, colsum
	
	#*****************************************************************************
	
	def all_sum( self ):
		keys = self.transition_probabilities.keys()
		keys.sort()
		
		sum = 0
		for k in keys:
			for l in keys:
				sum += self.transition_probabilities[k][l]
		
		print sum
		
	#*****************************************************************************
		
	def build( self, fastafile ):
		# count the transitions
		self.transition_counts = dict()
		self.initial_counts = dict()
		for seq_record in SeqIO.parse( fastafile, "fasta" ):
			sequence = str( seq_record.seq )
			if sequence[:3] == "Seq":
				print >> sys.stderr, "Warning: for %s sequence is unavailable..." % \
					seq_record.id
				continue
			i = 0
			while i <= len( sequence ) - 6:
				codon = sequence[i:i+3]
				next_codon = sequence[i+3:i+6]
				if codon.find( "N" ) >= 0 or next_codon.find( "N" ) >= 0 or \
					codon.find( "R" ) >= 0 or next_codon.find( "R" ) >= 0: # we have 'N' or 'R' (?)
					print >> sys.stderr, "Warning: found 'N'... truncating CDS."
					break # break because they will mess the reading frame!!! :-(
				# initial counts
				if i == 0:
					if codon not in self.initial_counts:
						self.initial_counts[codon] = 1
					else:
						self.initial_counts[codon] += 1
				if codon not in self.transition_counts:
					self.transition_counts[codon] = dict()
				if next_codon not in self.transition_counts[codon]:
					self.transition_counts[codon][next_codon] = 1
				else:
					self.transition_counts[codon][next_codon] += 1
				i += 3
	
		# the are some empty transitions e.g. Stop->Ci
		# we fill these with low but non-zero counts
		# add a pseudocount of 1 to everything
		itercodons = itertools.product( "ACGT", repeat=3 )
		all_codons = [ "".join( codon ) for codon in itercodons ]
	
		for C1 in all_codons:
			try:
				a_codon = self.transition_counts[C1]			# C1 is present
				for C2 in all_codons:
					try:
						self.transition_counts[C1][C2] += 1		# C1 and C2 are present
					except KeyError:
						self.transition_counts[C1][C2] = 1		# C1 present; C2 absent
			except KeyError:
				self.transition_counts[C1] = dict()				# neither C1 nor C2 are present
				for C2 in all_codons:
					self.transition_counts[C1][C2] = 1
		
		# do the same for initial counts
		for C1 in all_codons:
			try:
				self.initial_counts[C1] += 1
			except KeyError:
				self.initial_counts[C1] = 1
	
		# normalise the counts into probabilities
		self.transition_probabilities = dict()
		for codon in self.transition_counts:
			self.transition_probabilities[codon] = dict()
			total = sum( self.transition_counts[codon].values() )
			for next_codon in self.transition_counts[codon]:
				self.transition_probabilities[codon][next_codon] = \
					self.transition_counts[codon][next_codon]/total
		
		self.initial_probabilities = dict()
		total = sum( self.initial_counts.values())
		for codon in self.initial_counts:
			self.initial_probabilities[codon] = self.initial_counts[codon]/total
	
	#*****************************************************************************
	
	def probability( self, C1, C2, loglik=True, logbase=math.exp( 1 ) ): # probability of C2 given C1 i.e. C1->C2
		if loglik:
			return math.log( self.transition_probabilities[C1][C2], logbase )
		else:
			return self.transition_probabilities[C1][C2]
	
	#*****************************************************************************
	
	def likelihood( self, sequence, loglik=True, initial=False ): # sequence likelihood
		if loglik:
			loglikelihood = 0
			if len( sequence ) >= 3 and initial:
				loglikelihood += math.log( self.initial_probabilities[sequence[:3]] )
			i = 3
			while i <= len( sequence ) - 6:
				codon = sequence[i:i+3]
				next_codon = sequence[i+3:i+6]
				loglikelihood += self.probability( codon, next_codon, loglik )
				i += 3
			return loglikelihood
		else: # potential problem naming a variable like the function!!!
			likelihood = 1
			if len( sequence ) >= 3 and initial:
				likelihood *= self.initial_probabilities[sequence[:3]]
			i = 3
			while i <= len( sequence ) - 6:
				codon = sequence[i:i+3]
				next_codon = sequence[i+3:i+6]
				likelihood *= self.probability( codon, next_codon, loglik )
				i += 3
			return likelihood
	
	#*****************************************************************************
	
	def graded_likelihood( self, sequence, loglik=True, initial=False ): 
		# graded - cumulative likelihood across sequence
		if loglik:
			graded_loglikelihood = list()
			loglikelihood = 0
			if len( sequence ) >= 3 and initial:
				loglikelihood = math.log( self.initial_probabilities[sequence[:3]] )
				graded_loglikelihood.append( loglikelihood )
			i = 3
			while i <= len( sequence ) - 6:
				codon = sequence[i:i+3]
				next_codon = sequence[i+3:i+6]
				loglikelihood += self.probability( codon, next_codon, loglik )
				graded_loglikelihood.append( loglikelihood )
				i += 3
			return graded_loglikelihood
		else: # potential problem naming a variable like the function!!!
			graded_likelihood = list()
			likelihood = 1
			if len( sequence ) >= 3 and initial:
				likelihood = self.initial_probabilities[sequence[:3]]
				graded_likelihood.append( likelihood )
			i = 3
			while i <= len( sequence ) - 6:
				codon = sequence[i:i+3]
				next_codon = sequence[i+3:i+6]
				likelihood += self.probability( codon, next_codon, loglik )
				graded_likelihood.append( likelihood )
				i += 3
			return graded_likelihood
	
	#*****************************************************************************
	
	def differential_graded_likelihood( self, sequence, loglik=True, initial=False ): 
		# minus uniform random sequence
		if loglik:
			diff_graded_loglikelihood = list()
			graded_loglikelihood = self.graded_likelihood( sequence, loglik, \
				initial=initial )
			diff_graded_loglikelihood = [ graded_loglikelihood[i] + \
				( i + 1 )*math.log( 64 ) for i in xrange( len( graded_loglikelihood ))]
			return diff_graded_loglikelihood
		else: # potential problem naming a variable like the function!!!
			diff_graded_likelihood = list()
			graded_likelihood = self.graded_likelihood( sequence, initial=initial )
			diff_graded_likelihood = [ 64**( i + 1 )*graded_likelihood[i]]
			return diff_graded_likelihood
		
	#*****************************************************************************
	
	def likelihood_slope( self, dgl ):
		"""
		dgl = a differential gradient likelihood list
		"""
		if len( dgl ) == 0:
			return None
		elif len( dgl ) == 1:
			return 0
		else:
			return ( dgl[-1] - dgl[0] )/len( dgl )
	
	#*****************************************************************************
	
	def write( self, outfile ):
		data = self.transition_counts, self.transition_probabilities, \
			self.initial_counts, self.initial_probabilities
		with open( outfile, 'w' ) as f:
			cPickle.dump( data, f, cPickle.HIGHEST_PROTOCOL )
	
	#*****************************************************************************
		
	def read( self, infile ):
		with open( infile ) as f:
			data = cPickle.load( f )
			self.transition_counts, self.transition_probabilities, \
				self.initial_counts, self.initial_probabilities = data
