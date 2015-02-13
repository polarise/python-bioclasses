#!/usr/bin/env python
import sys
import BioClasses
import cPickle
import pysam
import scipy

def main():
	# load data from the GTF pic
	with open( "/home/paul/Resources/H_sapiens/test.hg19.chr.pic" ) as f:
		genes = cPickle.load( f )
		
	# load the indexed FASTA file
	fastafile = pysam.Fastafile( "/home/paul/Resources/H_sapiens/hg19.ens.chr.fa" )
	
	#  Nucl. Acids Res. (2008) 36 (11): 3707-3715. doi: 10.1093/nar/gkn248 
	matrix = """5	2	0	0	0	0	0	0	1	2	2
33	6	6	25	48	2	25	2	12	21	29
7	16	8	0	0	0	0	0	0	4	11
3	24	34	23	0	46	23	46	35	21	6"""

	M = list()
	for row in matrix.split( "\n" ):
		M.append( map( int, row.split( "\t" )))

	MM = scipy.array( M )
	P = BioClasses.PSWM( MM )
	
	print MM, MM.shape
	print
	print P
		
	# iterate over a few Transcript objects
	c = 0
	for gene_id,G in genes.iteritems():
		for transcript_id,T in G.transcripts.iteritems():
			if c > 100:
				break
			
			# arbitrary flanks
			print T
			T.TSS_Sequence.set_flank_right( 50 )
			T.TSS_Sequence.set_flank_left( 30 )
			
			# get the sequence before you compute TOP scores
			T.TSS_Sequence.get_sequence( fastafile )
			
			# compute TOP scores and get the highest
			T.TSS_Sequence.compute_TOP_score( P )
			
			print T.TSS_Sequence
			c += 1
	# 

if __name__ == "__main__":
	main()