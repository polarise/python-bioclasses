# -*- encoding: utf-8 -*-
import sys
import numpy as np
import math
import scipy as sp
import time

def PrintStatic( line, stderr=True ):
	"""
	Print to sdtderr without a newline at the end
	"""
	if stderr:
		sys.stderr.write( "\r%s".ljust( 20 ) % line )
		sys.stderr.flush()
	else:
		sys.stdout.write( "\r%s".ljust( 20 ) % line )
		sys.stdout.flush()
	
def msg( mesg, newline=True, show_time=True ):
	"""
	Neat function to write error messages to the screen
	"""
	if show_time:
		if newline:
			sys.stderr.write( "[%s] %s\n" % ( time.ctime( time.time()), mesg ))
		else:
			sys.stderr.write( "[%s] %s" % ( time.ctime( time.time()), mesg ))
	else:
		if newline:
			sys.stderr.write( mesg + "\n" )
		else:
			sys.stderr.write( mesg )
	sys.stderr.flush()

def savgol(x, window_size=3, order=2, deriv=0, rate=1):
	"""
	Savitzky-Golay filter
	Author: Thomas Haslwanter Version: 1.0 Date: 25-July-2012
	"""
	# Check the input
	try:
		window_size = np.abs(np.int(window_size))
		order = np.abs(np.int(order))
	except ValueError:
		raise ValueError("window_size and order have to be of type int")
	if window_size > len(x):
		raise TypeError("Not enough data points!")
	if window_size % 2 != 1 or window_size < 1:
		raise TypeError("window_size size must be a positive odd number")
	if window_size < order + 1:
		raise TypeError("window_size is too small for the polynomials order")
	if order <= deriv:
		raise TypeError("The 'deriv' of the polynomial is too high.")
	
	# Calculate some required parameters
	order_range = range(order+1)
	half_window = (window_size -1) // 2
	num_data = len(x)

	# Construct Vandermonde matrix, its inverse, and the Savitzky-Golay coefficients   
	a = [[ii**jj for jj in order_range] for ii in range(-half_window, half_window+1)]
	pa = np.linalg.pinv(a)
	sg_coeff = pa[deriv] * rate**deriv * math.	factorial(deriv)

	# Get the coefficients for the fits at the beginning and at the end of the data
	coefs = np.array(order_range)**np.sign(deriv)
	coef_mat = np.zeros((order+1, order+1))
	row = 0
	for ii in range(deriv,order+1):
		coef = coefs[ii]
		for jj in range(1,deriv):
			coef *= (coefs[ii]-jj)
			coef_mat[row,row+deriv]=coef
			row += 1
	coef_mat *= rate**deriv

	# Add the first and last point half_window times
	firstvals = np.ones(half_window) * x[0] 
	lastvals  = np.ones(half_window) * x[-1]
	x_calc = np.concatenate((firstvals, x, lastvals))

	y = np.convolve( sg_coeff[::-1], x_calc, mode='full')

	# chop away intermediate data
	y = y[window_size-1:window_size+num_data-1]

	# filtering for the first and last few datapoints
	y[0:half_window] = sp.dot(sp.dot(sp.dot(a[0:half_window], coef_mat), \
							np.mat(pa)), x[0:window_size])
	y[len(y)-half_window:len(y)] = sp.dot(sp.dot(sp.dot(\
		a[half_window+1:window_size],	coef_mat), pa), x[len(x)-window_size:len(x)])

	return y
