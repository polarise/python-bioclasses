#!/usr/bin/env python
from __future__ import division
import sys
from Node import *
from Branch import *
from Tree import *

def main():
	branches = list()
	# n1 = Node( 1, 0 ) # the root node
	# n1.add_descendant([ Node( 0, 1 ), Node( 2, 2 )])
	b1 = Branch( Node( 1, 0 ), Node( 0, 1 ), Node( 2, 2 ))
	branches.append( b1 )
	
	t = Tree()
	t.graft_branch( b1 )
	# t.graft( branches[n1.identify()] )
	
	
	print t
	print
	
	for tn in t.leaves:
		print tn, tn.parent, tn.parent.parent
	
	print
	
	for leaf in t.leaves:
		print leaf, t.path_to_leaf( leaf )
	
	print
	# n2 = Node( 0, 1 )
	# n2.add_descendant([ Node( 2, 2 ), Node( 1, 4 ) ])
	b2 = Branch( Node( 0, 1 ), Node( 2, 2 ), Node( 1, 4 ))
	branches.append( b2 )
	t.graft_branch( b2 )
	
	# t.graft( branches[n2.identify()] )
	
	print t
	print
	
	for tn in t.leaves:
		print tn, tn.parent, tn.parent.parent
	
	print
	
	for leaf in t.leaves:
		print leaf, t.path_to_leaf( leaf )
	print
	
	# n3 = Node( 2, 2 )
	# n3.add_descendant([ Node( 1, 4 ), Node( 0, 5 )])
	# t.graft( n3 )
	print "begin here:"
	b3 = Branch( Node( 2, 2 ), Node( 1, 4 ), Node( 0, 5 ))
	branches.append( b3 )
	t.graft_branch( b3 )
	print "end here:"
	
	print t
	print
	
	for tn in t.leaves:
		print tn, tn.parent, tn.parent.parent
	
	print
	
	for leaf in t.leaves:
		print leaf, t.path_to_leaf( leaf )
		
	# sys.exit()
	
	
	# n4 = Node( 1, 4 )
	# n4.add_descendant([ Node( 0, 5 ), Node( 2, -1 )])
	# t.graft( n4 )
	b4 = Branch( Node( 1, 4 ), Node( 0, 5 ), Node( 2, -1 ))
	branches.append( b4 )
	t.graft_branch( b4 )
	
	print t
	print

	# n5 = Node( 0, 5 )
	# n5.add_descendant([ Node( 1, 6 ), Node( 2, -1 )])
	# t.graft( n5 )
	b5 = Branch( Node( 0, 5 ), Node( 1, 6 ), Node( 2, -1 ))
	branches.append( b5 )
	t.graft_branch( b5 )
	
	
	print t
	print
	
	# n6 = Node( 1, 6 )
	# n6.add_descendant([ Node( 0, -1 ), Node( 2, -1 )])
	# t.graft( n6 )
	b6 = Branch( Node( 1, 6 ), Node( 0, -1 ), Node( 2, -1 ))
	branches.append( b6 )
	t.graft_branch( b6 )
	
	print t
	print
	
	for tn in t.leaves:
		print tn, tn.parent, tn.parent.parent, tn.parent.parent.parent, tn.parent.parent.parent.parent
	
	print
	
	# for n in t.leaves:
	# 	while n != t.root:
	# 		print n, n.name
	# 		n = n.parent
	# 	print n, n.name
	# 	print
	
	for leaf in t.leaves:
		print leaf, t.path_to_leaf( leaf )

	
if __name__ == "__main__":
	main()