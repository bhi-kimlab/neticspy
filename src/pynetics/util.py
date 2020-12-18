import numpy as np

from scipy.special import gammainc

def row_normalize(adj):
	return adj / adj.sum(axis=1).reshape(-1, 1)

def fdr(pval, qval=0.05, cV=1):
	# Accept arguments
	if cV == 0:
		cV = sum([1. / i for i in range(1, pval.size+1)])

	# Check if pval is a vector
	if pval.size != len(pval):
		raise ValueError('p-values should be a row or column vector, not an array.')

	# Check if pvals are within the interval
	if (min(pval) < 0) or (max(pval) > 1):
		raise ValueError('p-values out of range (0-1).')

	# Check if qval is withing the interval
	if (qval < 0) or (qval > 1):
		raise ValueError('q-value out of range (0-1).')

	# Sort p-values
	oidx = np.argsort(pval)
	pval.sort()

	# Number of observations
	V = pval.size

	# Order (indices), in the same size as the pvalues
	idx = np.reshape(np.arange(V), pval.size) + 1

	# Line to be used as cutoff
	thrline = idx * qval / (V * cV)

	try:
		# Find the largest pval, still under the line
		thr = max(pval[pval <= thrline])
		# Deal with the case when all the points under the line are equal to zero, and other points are above the line
		if thr == 0:
			thr = max(thrline[pval <= thrline])
	except ValueError:
		# Case when it does not cross
		thr = 0
	
	# Returns the result
	pthr = thr

	# p-corrected
	pcor = pval * V * cV / idx

	# Sort back to the original order and output
	oidxR = np.argsort(oidx)
	pcor = pcor[oidxR]

	# Loop over each sorted original p-value
	padj = np.zeros(pval.size)
	prev = 1
	for i in range(V-1, -1, -1):
		# The p-adjusted for the current p-value is the smallest slope among all the slopes of each of the p-values larger than the current one
		padj[i] = min(prev, pval[i] * V * cV / (i + 1))
		prev = padj[i]
	padj = padj[oidxR]
	return pthr, pcor, padj

def pchisq(x, a):
	# The chisquare distribution function
	# F = pchisq(x, Degrees Of Freedom)
	if a <= 0:
		raise ValueError('Degrees Of Freedom is wrong')
	F = pgamma(x/2, a*0.5)
	return F

def pgamma(x, a):
	# The gamma distribution function
	# F = pgamma(x, a)
	if a <= 0:
		raise ValueError('Parameter a is wrong')
	F = gammainc(a, x)
	I0 = (x < 0).nonzero() # tuple
	F[I0] = np.zeros(I0[0])
	return F