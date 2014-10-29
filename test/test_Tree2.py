#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
from Node import *
from Branch import *
from Tree import *

def main():
	t = Tree()
	b1 = Branch( Node( 1, 0 ), Node( 0, 1 ), Node( 2, 2 ))
	t.graft( b1 )
	print t
	print
	b2 = Branch( Node( 0, 1 ), Node( 2, 2 ), Node( 1, 4 ))
	t.graft( b2 )
	print t
	print
	b3 = Branch( Node( 2, 2 ), Node( 1, 4 ), Node( 0, 5 ))
	t.graft( b3 )
	print t
	print
	b4 = Branch( Node( 1, 4 ), Node( 0, 5 ), Node( 2, -1 ))
	t.graft( b4 )
	print t
	print
	b5 = Branch( Node( 0, 5 ), Node( 1, 6 ), Node( 2, -1 ))
	t.graft( b5 )
	print t
	print
	b6 = Branch( Node( 1, 6 ), Node( 0, -1 ), Node( 2, -1 ))
	t.graft( b6 )
	print t
	print
	
	paths = t.get_paths()
	print paths
	simple_paths = [ map( lambda x: ( x.name, x.value ), path ) for path in paths ]
	print len( paths )
	for p in simple_paths:
		print p

if __name__ == "__main__":
	main()