#!/usr/bin/env python
from __future__ import division
import sys
from Bio import SeqIO
from GlimmerORF import *

def main( glimmer_fn, fasta_fn ):
	# make the list of GlimmerORF objects
	glimmer_objs = dict()
	with open( glimmer_fn ) as f:
		name_state = True
		c = 0
		for row in f:
			if c >= 2000:
				break
			l = row.strip()
			if row[0] == ">":
				name = l[1:]
			else:
				G = GlimmerORF( name )
				if G.name not in glimmer_objs:	
					G.add_orf( l )
					glimmer_objs[G.name] = G
				else:
					glimmer_objs[G.name].add_orf( l )
			c += 0
	
	for seq_record in SeqIO.parse( fasta_fn, "fasta" ):
		seq_name = seq_record.id
		sequence = str( seq_record.seq )
		if seq_name in glimmer_objs:
			glimmer_objs[seq_name].extract_orf( sequence )
	
	for g in glimmer_objs:
		glimmer_objs[g].print_terminal_sequence()

if __name__ == "__main__":
	try:
		glimmer_fn = sys.argv[1]
		fasta_fn = sys.argv[2]
	except IndexError:
		print >> sys.stderr, "./scripy.py <glimmer-fn> <fasta-fn>"
		raise ValueError( "Missing files!" )
	main( glimmer_fn, fasta_fn )