version = '0.0.2'

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import beta
import argparse
import warnings

warnings.filterwarnings(action='ignore', category=RuntimeWarning)

def parse_args():
	parser = argparse.ArgumentParser(description='Python implementation of NetICS')

	parser.add_argument('--abberant', type=str, default=None, help='Input file that contains the genetically aberrant genes of each sample.\n It contain two columns that map every gene (1st column) to the samples that it it genetically aberrant (2nd column).')
	parser.add_argument('--adj', type=str, default=None, help='Adjacency matrix of the directed interaction network.')
	parser.add_argument('--beta', type=float, default=0.4, help='Restart probability for the insulated diffusion. For (Wu et al., 2010) network use 0.4.')
	parser.add_argument('--rank', type=str, default='SUM', help="'MEDIAN' uses the median of the sample-specific ranks.\n'RRA' uses the Robust Rank Aggregation method to integrate sample-specific ranked lists.\n'SUM' uses the sum of the sample-specific ranks.")
	parser.add_argument('--network', type=str, default=None, help='Input file that contains the list of the genes that are present in the network.\nThey should be in the same order as in the rows of the adjacency matrix adj. An example file is given that contains the gene names of the network described in (Wu et al., 2010).')
	parser.add_argument('--deg', type=str, default=None, help='Tab delimited file with two columns. First column contains the genes for which differential expression\nbetween the tumor and normal samples at the RNA level was measured. Second column contains the p-values\nof these measurements. This file can be the result of a tool for differential expression analysis such as DESeq2.\nEach gene in this file should have only one entry.')
	parser.add_argument('--dep', type=str, default=None, help='Tab delimited file with two columns. First column contain the proteins for which differential expression between\nthe tumor and normal samples at the protein level was measured. Second column contains the p-values of these\nmeasurements. Each gene in this file should have only one entry.')
	parser.add_argument('--output', type=str, default='NetICSpy_result.tsv', help='Output file of NetICSpy.')

	return parser.parse_args()

def main():
	args = parse_args()
	print(args)
	netics_fun(args.abberant, args.adj, args.beta, args.rank, args.network, args.deg, args.dep, args.output)

if __name__ == '__main__':
	main()

def netics_fun(filenameMu=None, filenameAdj=None, restart_prob=0.4, rank_method_str='SUM', filenameNet=None, filenameRNA=None, filenamePR=None, output=None):
	# Accept arguments and check input
	if None in [filenameMu, filenameAdj, filenameNet]:
		raise ValueError('Missing input arguments.')
	elif rank_method_str not in ['SUM', 'MEDIAN', 'RRA']:
		print('Wrong input: %s' % rank_method_str)
		print('Input should be MEDIAN, RRA or SUM')
		return None
	else:
		mutation_data = read_mutations(filenameMu)
		fid = open(filenameAdj, 'r')
		adj = np.loadtxt(fid, delimiter='\t')
		fid.close()
		if filenameRNA == None:
			choose_mut_diff = 1
		else:
			choose_mut_diff = 3

	# Read network genes
	fidnet = open(filenameNet, 'r')
	network_genes = np.loadtxt(fidnet, delimiter='\n', dtype=str)
	fidnet.close()

	# Copy data from inpurt mutation_data
	samples = {}
	samples['DNA'] = []
	samples['DNA_network'] = []
	for i in range(len(mutation_data)):
		samples['DNA'].append(np.char.upper(mutation_data[i]))
		samples['DNA_network'].append(np.char.upper(np.intersect1d(network_genes, samples['DNA'][i])))
	del mutation_data

	# Read differentially expressed genes
	diff_expr_genes = []
	if choose_mut_diff != 1:
		diff_expr_genes = read_diff_expr(filenameRNA, filenamePR)
		if (diff_expr_genes is None) or (len(diff_expr_genes) == 0):
			choose_mut_diff == 1
	else:
		print('No input files was provided for differentially expressed genes.')
	samples['RNA'] = []
	samples['RNA_network'] = []
	for i in range(len(samples['DNA'])):
		samples['RNA'].append(diff_expr_genes)
		samples['RNA_network'].append(np.char.upper(np.intersect1d(network_genes, samples['RNA'][i]).astype(str)))

	# Compute diffused matrix
	print('Computing diffused matrix...')
	F = insulated_diff(norm_adj(adj), restart_prob)
	F_opp = insulated_diff(norm_adj(adj.conj().transpose()), restart_prob)
	print('Running NetICS')
	ranked_list_genes, scores = prioritization(samples, F, F_opp, network_genes, choose_mut_diff, rank_method_str)

	if output is not None:
		res = pd.DataFrame([ranked_list_genes, scores], index=['Genes', 'Scores']).T
		res.to_csv(output, sep='\t', header=True, index=False)

	return ranked_list_genes, scores

def read_mutations(filename):
	fid = open(filename, 'r')
	g = np.loadtxt(fid, delimiter='\t', dtype=str)
	fid.close()
	unique_samples = np.unique(g[:,1])
	mutation_data = []
	for i in range(len(unique_samples)):
		mutation_data.append(g[:,0][np.isin(g[:,1], unique_samples[i])])
	return mutation_data # len=81

def read_diff_expr(filenameDE, filenamePR):
	diff_expr_genes = []
	RNA_names = []
	rppa_names = []

	if filenamePR != None:
		rppa_names, rppa_pval = read_file(filenamePR)
	else:
		print('No input file was provided for differentially expressed genes at the proteome level.')

	if filenameDE != None:
		RNA_names, RNA_pval = read_file(filenameDE)
	else:
		print('No input file was provided for differentially expressed genes at the RNA level.')

	# Check for duplicate genes in the input files
	if (len(np.unique(RNA_names)) != len(RNA_names)) or (len(np.unique(rppa_names)) != len(rppa_names)):
		print('Input files for differentially expressed genes should contain only one entry per gene.')
		print('Input files for differentially expressed genes were ignored.')
		return None
	# Check if the input files for differential expression are empty
	if (len(rppa_names) == 0) and (len(RNA_names) == 0):
		print('Only genetically aberrant genes are used for diffusion.')
		return None
	# Only diffferential expressed genes at the RNA level are provided
	elif (len(rppa_names) == 0) and (len(RNA_names) != 0):
		_, _, pval_all = fdr(RNA_pval, 0.05)
		diff_expr_genes = RNA_names[pval_all < 0.05]
	# Only differential expressed genes at the protein level are provided
	elif (len(RNA_names) == 0) and (len(rppa_names) != 0):
		_, _, pval_all = fdr(rppa_pval, 0.05)
		diff_expr_genes = rppa_names[pval_all < 0.05]
	else:
		names_intersection = np.intersect1d(RNA_names, rppa_names)
		names_setdiff_rna = np.setdiff1d(RNA_names, rppa_names)
		names_setdiff_pr = np.setdiff1d(rppa_names, RNA_names)

		rna_pval_names_intersection = RNA_pval[np.isin(RNA_names, names_intersection)].astype(np.float)
		rppa_pval_names_intersection = rppa_pval[np.isin(rppa_names, names_intersection)].astype(np.float)

		rna_pval_names_setdiff_rna = RNA_pval[np.isin(RNA_names, names_setdiff_rna)]
		pr_pval_names_setdiff_pr = rppa_pval[np.isin(rppa_names, names_setdiff_pr)]

		all_names = np.concatenate((names_intersection, names_setdiff_rna, names_setdiff_pr))
		# Fishers method
		pval_all = (1 - pchisq(-2 * (np.log(rna_pval_names_intersection) + np.log(rppa_pval_names_intersection)), 4))
		pval_all = np.concatenate((pval_all, rna_pval_names_setdiff_rna, pr_pval_names_setdiff_pr))
		_, _, pval_all = fdr(pval_all, 0.05)
		diff_expr_genes = all_names[pval_all < 0.05]
	return diff_expr_genes

def read_file(filenamePR):
	fid = open(filenamePR, 'r')
	g = pd.read_csv(fid, delimiter='\t', header=None).values
	fid.close()
	names = g[:,0]
	pval = g[:,1]
	return names, pval

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

def norm_adj(adj):
	degree = adj.sum(axis=1)
	W = np.zeros((len(adj[:,0]), len(adj[0,:])))
	for i in range(len(adj[:,0])):
		for j in range(len(adj[0,:])):
			if adj[i,j] != 0:
				W[i,j] = 1 / degree[i]
	return W

def insulated_diff(W, b):
	temp = np.identity(len(W)) - (1-b) * W
	F = b * np.linalg.inv(temp)
	return F

def prioritization(samples, F, F_opp, network_genes, choose_up_down_flag, rank_method_str):
	weights_all = np.zeros((len(samples['DNA_network']), len(network_genes)))
	positions_sorted = np.zeros((len(network_genes), len(samples['DNA_network'])))
	final_weights_median = np.zeros(len(network_genes))
	final_weights_sum = np.zeros(len(network_genes))
	final_weights_rho = np.zeros(len(network_genes))

	for i in range(len(samples['DNA_network'])): # Number of samples
		weights_all[i,:] = diffuse_all(choose_up_down_flag, samples['DNA_network'][i], samples['RNA_network'][i], network_genes, F, F_opp)
		to_sort_sorted = sortrows(np.array([weights_all[i,:], np.arange(len(weights_all[i,:]))]).conj().transpose())
		positions_sorted[:,i] = np.flip(to_sort_sorted[:,1].conj().transpose())
	for j in range(len(network_genes)): # Number of genes
		rank = np.zeros((len(samples['DNA_network']), 1))
		for i in range(len(samples['DNA_network'])): # Number of samples
			#rank[i] = (positions_sorted[:,i] == j).nonzero()
			rank[i] = (positions_sorted[:,i] == j).nonzero()[0] + 1
		final_weights_median[j] = np.median(rank)
		final_weights_sum[j] = rank.sum(axis=0)
		final_weights_rho[j] = rhoScores(rank / max(rank))

	# Choose the aggregation method
	if rank_method_str == 'MEDIAN':
		scores = final_weights_median.astype(np.int64)
	elif rank_method_str == 'SUM':
		scores = final_weights_sum.astype(np.int64)
	else:
		scores = final_weights_rho

	# Sort scores
	to_sort_sorted_rho = sortrows(np.array([scores, np.arange(len(scores))]).conj().transpose())
	ranked_list_genes = network_genes[to_sort_sorted_rho[:,1].astype(np.int64)]
	ranked_scores = scores[to_sort_sorted_rho[:,1].astype(np.int64)]
	return ranked_list_genes, ranked_scores

def diffuse_all(choose_up_down_flag, sample_DNA, sample_RNA, network_genes, F, F_opp):
	# Diffuse mutation labels
	if choose_up_down_flag == 1 or choose_up_down_flag == 3:
		mutation_weights = diffuse(np.intersect1d(sample_DNA, network_genes), network_genes, F)
	# Diffuse DE labels
	if choose_up_down_flag == 2 or choose_up_down_flag == 3:
		DE_weights = diffuse(np.intersect1d(sample_RNA, network_genes), network_genes, F_opp)
	# Combines scores by multiplication (use both mutations and DE)
	if choose_up_down_flag == 3:
		ret_weights_all = mutation_weights * DE_weights
	elif choose_up_down_flag == 1:
		ret_weights_all = mutation_weights
	else:
		ret_weights_all = DE_weights
	return ret_weights_all

def diffuse(inter_mut_sampl_network, network_genes, F):
	positions = (np.isin(network_genes, inter_mut_sampl_network)).nonzero()[0] # len=6
	if len(positions) != 0:
		mutation_weights = (1 / len(positions)) * np.matmul(np.ones((1, len(positions))), F[positions,:])
	else:
		mutation_weights = np.inf * np.matmul(np.ones((1, len(positions))), F[positions,:])
	return mutation_weights

def sortrows(M):
	if len(M.shape) != 2:
		raise ValueError('M must be 2d numpy.array')
	M_columns = tuple(M[:,c] for c in range(M.shape[1]-1, -1, -1))
	return M[np.lexsort(M_columns), :]

def rhoScores(r, topCutoff=None):
	r = r.reshape((1,-1))

	rho = np.empty(1)
	rho[:] = np.nan

	r1 = r[0,:]
	x = betaScores(r1)[0]
	# Correct using Bonferroni method.
	rho[0] = correctBetaPvalues(min(x), sum(~np.isnan(x)))
	return rho

def betaScores(r):
	n = sum(~np.isnan(r))
	p = np.empty((1,len(r)))
	p[:] = np.nan
	# Sort the values
	r.sort()
	# Get the order statistics and calculates p-values for each of the order statistics. These are based on their expected distribution under the null hypothesis of uniform distribution.
	p[:,:n] = beta.cdf(r[0:n], np.arange(1, n+1), np.arange(n, 0, -1))
	return p

def correctBetaPvalues(p, k):
	pval = beta.cdf(p, 1, k)
	return pval
