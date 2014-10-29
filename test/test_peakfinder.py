#!/usr/bin/env python
from __future__ import division
import sys
import random
import termcolor

peak_pvalues = list()
for i in xrange( 225 ):
	if random.random() < 0.5:
		peak_pvalues.append( 0 )
	else:
		peak_pvalues.append( 1 )

# look for peaks
found_peaks = dict()
peak = list()
i = 0
n = 4	
peak_pos = -1
for p in peak_pvalues:
	if p == 0:
		if peak_pos < 0:
			peak_pos = i
		peak.append( p )	
	elif p == 1:
		if len( peak ) >= n:
			found_peaks[peak_pos] = len( peak )
		peak = list()
		peak_pos = -1
	i += 1

sequence = "".join( map( str, peak_pvalues ))
print found_peaks
print sequence
keys = found_peaks.keys()
keys.sort()
for k in keys:
	v = found_peaks[k]
	print " "*k + sequence[k:k+v]
print sequence


