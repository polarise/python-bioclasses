# -*- encoding: utf-8 -*-
from __future__ import division
from AminoAcid import *
from Bio import SeqIO

class GeneticCode( object ):
	def __init__( self, gc_file ):
		self.gc_file = gc_file
		self.codons = dict()
		self.amino_acids = dict()
		self.has_CAI = False
		with open( self.gc_file ) as f:
			for row in f:
				l = row.strip( "\n" ).split( "\t" )
				self.amino_acids[l[0]] = AminoAcid( *l )
				for c in l[3].split( "," ):
					self.codons[c] = l[0]
	
	def __repr__( self ):
		gc_str = ""
		for a in self.amino_acids:
			gc_str += str( self.amino_acids[a] ) + "\n\n"
		return gc_str
	
	def build_CAI_table( self, cds_file, DNA=True ):
		for seq_record in SeqIO.parse( cds_file, "fasta" ):
			sequence = str( seq_record.seq )
			if sequence[:3] == "Seq":
				print >> sys.stderr, "Sequence unavailable for %s..." % seq_record.id
				continue
			if DNA:
				sequence = sequence.replace( "T", "U" ) # convert DNA to RNA
			else:
				sequence = sequence.replace( "U", "T" ) # convert RNA to DNA
			i = 0
			while i <= len( sequence ) - 3:
				codon = sequence[i:i+3]
				if codon.find( "N" ) >= 0 or codon.find( "R" ) >= 0:
					print >> sys.stderr, "Warning: founding wrong base in codon..."
					break
				aa = self.codons[codon]
				self.amino_acids[aa].count_codon( codon )
				i += 3
		
		for g in self.amino_acids:
			self.amino_acids[g].normalise()
		
		self.has_CAI = True
	
	def write_CAI_table( self, CAI_file ):
		aa_names = self.amino_acids.keys()
		aa_names.sort()
		with open( CAI_file, 'w' ) as f:
			print >> f, "\t".join([ "Amino acid", "Lsym", "Ssym", "Codon", "#Codons", "RSCU",\
			"w_ij" ])
			for aa in aa_names:
				A = self.amino_acids[aa]
				for codon in A.codons:
					print >> f, "\t".join( map( str, [ A.name, A.lsymbol, A.ssymbol, codon, \
						A.codon_counts[codon], A.norm_counts[codon], A.wij[codon] ]))
	
	def read_CAI_table( self, CAI_file ):
		with open( CAI_file ) as f:
			for row in f:
				l = row.strip( "\n" ).split( "\t" )
				if l[0] == "Amino acid":
					continue
				# Alanine	Ala	A	GCU	37810	0.409598093381	1.0
				name, lsym, ssym, codon, count, RSCU, wij = l
				self.amino_acids[name].set_count( codon, int( count ))
		
		for aa in self.amino_acids:
			self.amino_acids[aa].normalise()
		
		self.has_CAI = True
	
	def get_wij( self, codon ):
		codon = codon.replace( "T", "U" )
		aa = self.codons[codon]
		if self.has_CAI:
			return self.amino_acids[aa].wij[codon]
		else:
			return None

	
	