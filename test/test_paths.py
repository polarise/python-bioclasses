#!/usr/bin/env python
import sys
from Node import *	
from Branch import *

class Paths( object ):
	def __init__( self ):
		self.paths = list()
	
	def extend( self, branch ):
		new_paths = list()
		if len( self.paths ) == 0:
			for d in branch.descendants:
				new_paths.append([ d, branch.root ] )
		else:
			for d in branch.descendants:
				for p in self.paths:
					if p[0] == branch.root:
						new_paths.append( [d] + p )
					else:
						if p not in new_paths:
							new_paths.append( p )
		
		self.paths = new_paths
	
	def get_path_sequences( self, leaf ):
		path_sequences = list()
		for p in self.paths:
			if p[0] == leaf:
				path_sequences.append( map( lambda x: x.name, p[::-1] ))
		
		return path_sequences
	
	def get_all_path_sequences( self ):
		all_path_sequences = list()
		for p in self.paths:
			all_path_sequences.append( map( lambda x: x.name, p[::-1] ))
		
		return all_path_sequences

def extend( paths, branch ):
	new_paths = list()
	if len( paths ) == 0:
		for b in branch[1:]:
			new_paths.append([ b, branch[0] ])
	else:
		for b in branch[1:]:
			for p in paths:
				if p[0] == branch[0]:
					new_paths.append( [b] + p )
				else:
					if p not in new_paths:
						new_paths.append( p )
	
	return new_paths

paths = list()
b1 = [(1,0), (0,1), (2,2)]
b11 = Branch( Node( 1, 0 ), Node( 0, 1 ), Node( 2, 2 ))
paths = extend( paths, b1 )
print b1
print paths, len( paths )
print
b2 = [(0,1), (2,2), (1,4)]
b21 = Branch( Node( 0, 1 ), Node( 2, 2 ), Node( 1, 4 ))
paths = extend( paths, b2 )
print b2
print paths, len( paths )
print
b3 = [(2,2), (1,4), (0,5)]
b31 = Branch( Node( 2, 2 ), Node( 1, 4 ), Node( 0, 5 ))
paths = extend( paths, b3 ) 
print b3
print paths, len( paths )
print
b4 = [(1,4), (0,5), (2,-1)]
b41 = Branch( Node( 1, 4 ), Node( 0, 5 ), Node( 2, -1 ))
paths = extend( paths, b4 )
print b4
print paths, len( paths )
print
b5 = [(0,5), (1,6), (2,-1)]
b51 = Branch( Node( 0, 5 ), Node( 1, 6 ), Node( 2, -1 ))
paths = extend( paths, b5 )
print b5
print paths, len( paths )
print
b6 = [(1,6), (0,-1), (2,-1)]
b61 = Branch( Node( 1, 6 ), Node( 0, -1 ), Node( 2, -1 ))
paths = extend( paths, b6 )
print b6
print paths, len( paths )
print

p = Paths()
p.extend( b11 )
p.extend( b21 )
p.extend( b31 )
p.extend( b41 )
p.extend( b51 )
p.extend( b61 )
all_paths = p.get_all_path_sequences()
print all_paths
print len( all_paths )
print p.paths, len( p.paths )
sys.exit()


path_str = list()
for p in paths:
	one_path = ""
	for q in p:
		one_path = str( q[0] ) + one_path
	path_str.append( one_path )

path_str.sort()

for p in path_str:
	print p
	
	

