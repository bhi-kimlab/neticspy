__version__ = '0.0.0'

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.special import gammainc
from scipy.stats import beta

def netics_fun(filenameMu=None, filenameAdj=None, restart_prob=0.4, rank_method_str='SUM', filenameNet=None, filenameRNA=None, filenamePR=None):
	ranked_list_genes = {}

	# Check input
	if None in [filenameMu, filenameAdj, filenameNet]:
		raise ValueError('\nError: Missing input arguments.')
	elif rank_method_str not in ['SUM', 'MEDIAN', 'RRA']:
		print('Wrong input: %s\nInput should be MEDIAN, RRA, SUM' % rank_method_str)
		return None
	else:
		mutation_data = read_mutations(filenameMu)
		fidadj = open(filenameAdj, 'r')
		adj = np.loadtxt(fidadj, delimiter='\t')
		fidadj.close()
		choose_mut_diff = 1

	if filenameRNA != None:
		choose_mut_diff = 3

	# read network genes
	fidnet = open(filenameNet, 'r')
	network_genes = np.loadtxt(fidnet, delimiter='\n', dtype=str)
	fidnet.close()

	# copy data from input mutation_data
	samples = {}
	samples['DNA'] = []
	samples['DNA_network'] = []
	for i in range(len(mutation_data)):
		samples['DNA'].append(np.char.upper(mutation_data[i]))
		samples['DNA_network'].append(np.intersect1d(network_genes, samples['DNA'][i]))
	del mutation_data

	# read differentially expressed genes
	diff_expr_genes = []
	if choose_mut_diff != 1:
		diff_expr_genes = read_diff_expr(filenameRNA, filenamePR)
		if len(diff_expr_genes) == 0:
			choose_mut_diff = 1
			diff_expr_genes = ['']
	else:
		print('No Input files were provided for differentially expressed genes.')

	samples['RNA'] = []
	samples['RNA_network'] = []
	for i in range(len(samples['DNA'])):
		samples['RNA'].append(diff_expr_genes)
		samples['RNA_network'].append(np.char.upper(np.intersect1d(network_genes, samples['RNA'][i])))

	# compute diffused matrix
	print('Computing diffused matrix...')
	F = insulated_diff(norm_adj(adj), restart_prob)
	F_opp = insulated_diff(norm_adj(adj.conj().transpose()), restart_prob)
	print('Running NetICS...')
	ranked_list_genes, scores = prioritization(samples, F, F_opp, network_genes, choose_mut_diff, rank_method_str)
	return ranked_list_genes, scores

def read_mutations(filename):
	fid = open(filename, 'r')
	g = np.loadtxt(fid, delimiter='\t', dtype=str)
	fid.close()
	unique_samples = np.unique(g[:,1])
	mutation_data = []
	for i in range(len(unique_samples)):
		mutation_data.append(g[g[:,1]==unique_samples[i],0].reshape(-1,1))
	return mutation_data

def read_diff_expr(filenameDE=None, filenamePR=None):
	diff_expr_genes = []
	RNA_names = []
	rppa_names = []

	if filenamePR is not None:
		rppa_names, rppa_pval = read_file(filenamePR)
	else:
		print('No Input file was provided for differentially expressed genes at the proteome level.')

	if filenameDE is not None:
		RNA_names, RNA_pval = read_file(filenameDE)
	else:
		print('No Input file was provided for differentially expressed genes at the RNA level.')

	# check for duplicate genes in the input files
	if (len(np.unique(RNA_names)) != len(RNA_names)) or (len(np.unique(rppa_names)) != len(rppa_names)):
		print('Input files for differentially expressed genes should contain only one entry per gene.\nInput files for differentially expressed genes were ignored.')
		return None

	# check if the input files for differential expression are empty
	if rppa_names == [] and RNA_names == []:
		print('Only genetically aberrant genes are used for diffusion.')
		return None
	# only differential expressed genes at the RNA level are provided
	elif rppa_names == [] and RNA_names != []:
		pval_all = multipletests(RNA_pval, alpha=0.05, method='fdr_bh')[1]
		diff_expr_genes = RNA_names[pval_all < 0.05]
	# only differential expressed genes at the protein level are provided
	elif RNA_names == [] and rppa_names != []:
		pval_all = multipletests(rppa_pval, alpha=0.05, method='fdr_bh')[1]
		diff_expr_genes = rppa_names[pval_all < 0.05]
	else:
		names_intersection = np.intersect1d(RNA_names, rppa_names)
		names_setdiff_rna = np.setdiff1d(RNA_names, rppa_names)
		names_setdiff_pr = np.setdiff1d(rppa_names, RNA_names)

		rna_pval_names_intersection = []
		for i in range(len(RNA_names)):
			if RNA_names[i] in names_intersection:
				rna_pval_names_intersection.append(RNA_pval[i])
		rppa_pval_names_intersection = []
		for i in range(len(rppa_names)):
			if rppa_names[i] in names_intersection:
				rppa_pval_names_intersection.append(rppa_pval[i])

		rna_pval_names_setdiff_rna = []
		for i in range(len(RNA_names)):
			if RNA_names[i] in names_setdiff_rna:
				rna_pval_names_setdiff_rna.append(RNA_pval[i])
		pr_pval_names_setdiff_pr = []
		for i in range(len(rppa_names)):
			if rppa_names[i] in names_setdiff_pr:
				pr_pval_names_setdiff_pr.append(rppa_pval[i])
		all_names = np.concatenate((names_intersection, names_setdiff_rna, names_setdiff_pr))
		# Fishers method
		pval_all = (1 - pchisq(-2 * (np.log(rna_pval_names_intersection) + np.log(rppa_pval_names_intersection)), 4))
		pval_all = np.concatenate((pval_all, rna_pval_names_setdiff_rna, pr_pval_names_setdiff_pr))

		pval_all = multipletests(pval_all, alpha=0.05, method='fdr_bh')[1]
		diff_expr_genes = all_names[pval_all < 0.05]
	return diff_expr_genes

def read_file(filenamePR):
	fid = open(filenamePR, 'r')
	g =  pd.read_csv(fid, header=None, sep='\t').values
	fid.close()
	names = g[:,0]
	pval = g[:,1]
	return names, pval

def pchisq(x, a):
	F = pgamma(x/2, a*0.5)
	return F

def pgamma(x, a):
	F = gammainc(a, x)
	I0 = (x < 0).nonzero()
	F[I0] = np.zeros(len(I0))
	return F

def norm_adj(adj):
	degree =  adj.sum(axis=1)
	W = np.zeros((len(adj[:,0]),len(adj[0,:])))
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

	for i in range(len(samples['DNA_network'])):
		weights_all[i,:] = diffuse_all(choose_up_down_flag, samples['DNA_network'][i], samples['RNA_network'][i], network_genes, F, F_opp)
		to_sort_sorted = sortrows(np.array([weights_all[i,:], np.arange(len(weights_all[i,:]))]).conj().transpose())
		positions_sorted[:,i] = np.flip(to_sort_sorted[:,1].conj().transpose()).conj().transpose()
	for j in range(len(network_genes)):
		rank = np.zeros((len(samples['DNA_network']), 1))
		for i in range(len(samples['DNA_network'])):
			rank[i] = (positions_sorted[:,i] == j).nonzero()
		final_weights_median[j] = np.median(rank, axis=0)
		final_weights_sum[j] = rank.sum(axis=0)
		final_weights_rho[j] = rhoScores(rank/max(rank))

	# choose the aggregation method
	if rank_method_str == 'MEDIAN':
		scores = final_weights_median
	elif rank_method_str == 'SUM':
		scores = final_weights_sum
	else:
		scores = final_weights_rho

	# sort scores
	to_sort_sorted_rho = sortrows(np.array([scores, np.arange(len(scores))]).conj().transpose()).astype(np.int64)
	ranked_list_genes = network_genes[to_sort_sorted_rho[:,1]]
	ranked_scores = scores[to_sort_sorted_rho[:,1]].astype(np.int64)
	return ranked_list_genes, ranked_scores

def diffuse_all(choose_up_down_flag, sample_DNA, sample_RNA, network_genes, F, F_opp):
	# diffuse mutation labels
	if choose_up_down_flag == 1 or choose_up_down_flag == 3:
		mutation_weights = diffuse(np.intersect1d(sample_DNA, network_genes), network_genes, F)
	# diffuse DE labels
	if choose_up_down_flag == 2 or choose_up_down_flag == 3:
		DE_weights = diffuse(np.intersect1d(samples_RNA, network_genes), network_genes, F_opp)
	# combines scores by multiplication (use both mutations and DE)
	if choose_up_down_flag == 3:
		ret_weights_all = mutation_weights * DE_weights
	elif choose_up_down_flag == 1:
		ret_weights_all = mutation_weights
	else:
		ret_weights_all = DE_weights
	return ret_weights_all

def diffuse(inter_mut_sampl_network, network_genes, F):
	positions = np.isin(network_genes, inter_mut_sampl_network).ravel().nonzero()
	if len(positions[0]) != 0:
		mutation_weights = (1 / len(positions[0])) * np.matmul(np.ones((1, len(positions[0]))), F[positions,])
	else:
		mutation_weights = np.nan_to_num(np.inf * np.matmul(np.ones((1, len(positions[0]))), F[positions,]))
	return mutation_weights

def sortrows(M):
	M_columns = tuple(M[:,c] for c in range(M.shape[1]-1, -1, -1))
	return M[np.lexsort(M_columns), :]

def rhoScores(r):
	r = r.transpose()
	rho = np.array([np.NaN])
	r1 = r[0,:]
	x = betaScores(r1)
	rho[0] = correctBetaPvalues(min(x), sum(~np.isnan(x)))
	return rho

def betaScores(r):
	n = sum(~np.isnan(r))
	p = np.empty((r.shape[0]))
	p[:] = np.NaN
	r.sort()
	p[0:n] = beta.cdf(r[0:n], np.arange(1,n+1), np.arange(n,0,-1))
	return p

def correctBetaPvalues(p, k):
	pval = beta.cdf(p, 1, k)
	return pval
