#!/usr/bin/env python
from __future__ import division
import sys
from Node import *
from Branch import *

def main():
	b1 = Branch( Node( 1, 0 ), Node( 0, 1 ), Node( 2, 2 ))
	print b1.leaves, len( b1.leaves )

	b2 = Branch( Node( 0, 1 ), Node( 2, 2 ), Node( 1, 4 ))
	b1.graft_branch( b2 )
	print b1.leaves, len( b1.leaves )
	
	b3 = Branch( Node( 2, 2 ), Node( 1, 4 ), Node( 0, 5 ))
	b1.graft_branch( b3 )
	print b1.leaves, len( b1.leaves )
	
	b4 = Branch( Node( 1, 4 ), Node( 0, 5 ), Node( 2, -1 ))
	b1.graft_branch( b4 )
	print b1.leaves, len( b1.leaves )
	
	b5 = Branch( Node( 0, 5 ), Node( 1, 6 ), Node( 2, -1 ))
	b1.graft_branch( b5 )
	print b1.leaves, len( b1.leaves )
	
	b6 = Branch( Node( 1, 6 ), Node( 0, -1 ), Node( 2, -1 ))
	b1.graft_branch( b6 )
	print b1.leaves, len( b1.leaves )
	
	for leaf in b1.leaves:
#		print leaf, leaf.parent, leaf.parent.parent, leaf.parent.parent.parent, leaf.parent.parent.parent.parent
		print leaf, b1.path_to_leaf( leaf )

if __name__ == "__main__":
	main()
