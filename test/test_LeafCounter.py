#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import division
import sys
from Node import *
from LeafCounter import *

def main():
	L = LeafCounter()
	
	stop_seq = [ 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 ]
	nodes = [ Node( *d ) for d in zip( stop_seq, range( len( stop_seq ))) ]
	
	for n in nodes:
		L.add_node( n )
		print L.leaf_count()

if __name__ == "__main__":
	main()	